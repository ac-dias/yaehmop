/*******************************************************

Copyright (C) 2026 Alexandre C. Dias
All rights reserved

This file is part of yaehmop (fork ac-dias/yaehmop, branch adias-dev).

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

********************************************************************/

/*******************************************************
 * Spin-orbit coupling for the tightbinding program (adias-dev).
 *
 * On-site approximation, as in DFTB+ and in the yaehmop-bandfit (soc.py):
 *
 *     H_SOC = sum over atoms A and shells l >= 1 of  lambda_{A,l} L.S
 *
 * acting only between orbitals of the same shell of the same atom.  The
 * spinor basis has 2n functions: the n orbitals with spin up, then the same
 * n orbitals with spin down.  With S = sigma/2,
 *
 *     up-up = +(lambda/2) Lz            dn-dn = -(lambda/2) Lz
 *     up-dn =  (lambda/2) (Lx - i Ly)   dn-up =  (lambda/2) (Lx + i Ly)
 *
 *     H(k) -> I2 (x) H(k) + H_SOC,      S(k) -> I2 (x) S(k)
 *
 * (orbitals of one atom don't overlap each other, so H_SOC needs no overlap
 * correction).  p levels split into j = 3/2 (+lambda/2) and j = 1/2 (-lambda),
 * d levels into j = 5/2 (+lambda) and j = 3/2 (-3 lambda/2).
 *
 * The real orbitals are the YAeHMOP ones (px py pz; dx2-y2 dz2 dxy dxz dyz;
 * fz3 fxz2 fyz2 fxyz fz(x2-y2) fx(x2-3y2) fy(3x2-y2)), built from complex
 * spherical harmonics with the Condon-Shortley phase.
 *
 * Complex hermitian matrices use the packed layout of the program: element
 * (i,j), i<j, of the physical matrix is mat[i*n+j] - I*mat[j*n+i] (this is
 * what the LAPACK copy in kpoints.c/bands.c and the readers of .HAM/.OV in
 * the yaehmop-bandfit assume), the diagonal is real.
 *
 * The spinor states (2n per k point) go through the usual average
 * properties code with num_spin_comp = 2: the spin-up part of state J is
 * kept in avg_prop_info.orbs, the spin-down part in orbs_dn, and the charge
 * matrix chg_mat[J][mu] (2n x n) adds the two spin parts of orbital mu.
 * Since S_spinor = I2 (x) S does not couple the spins, Mulliken charges,
 * projected DOS and COOPs of a spinor state are exactly the sums of those of
 * its spin-up and spin-down parts.
 *
 * Input file (--soc FILE), one line per element and shell, ';' or '#' start
 * comments:
 *      Mo  d  0.080
 *      S   p  0.050
 *******************************************************/
#include "bind.h"

extern atom_type custom_atoms[];   /* fileio.c: parameters given in the input file */
#ifndef MAX_CUSTOM_ATOMS
#define MAX_CUSTOM_ATOMS 40         /* as in fileio.c */
#endif

#define SOC_MAX_ENTRIES 512

typedef struct {
  char symb[ATOM_SYMB_LEN];
  int l;
  real lambda;
  int used;
} soc_entry_type;

typedef struct { double re, im; } zc_type;

static int soc_dim = 0;            /* 2*num_orbs */
static real *socR = 0, *socI = 0;  /* dense H_SOC (2n x 2n) */
static real *Hs = 0, *Ss = 0;      /* packed spinor H(k), S(k) */
static complex *zH = 0, *zS = 0, *zwork = 0;
static real *soc_vals = 0, *soc_rwork = 0;
static int soc_ham_file = -1, soc_ov_file = -1;
/* eigenvectors (state J, orbital mu): spin-up and spin-down parts, n2 x n */
static real *vupR = 0, *vupI = 0, *vdnR = 0, *vdnI = 0;
static real *soc_chg = 0;          /* charge matrix, n2 x n */
static real *soc_occ = 0;          /* occupations of the n2 states (one k point) */
static real *SdR = 0, *SdI = 0;    /* dense S(k), n x n */
static real *soc_tmp = 0;          /* work, n2 x n2 */
static real *soc_spin = 0;         /* spin-up minus spin-down weight of each AO, n2 x n */
/* SOC + spin polarization: second spin channel */
static int soc_two_tables = 0;
static atom_type *dn_atoms = 0;    /* parameters of every atom, spin-down table */
static atom_type *up_prev = 0;     /* spin-up parameters at the last soc_follow_down */
static int soc_diff_zetas = 0;     /* the two tables have different zetas */
static real const_dn = -1.0;       /* K of the spin-down channel */
static hermetian_matrix_type HdnR = {0,0}, HdnK = {0,0};
static hermetian_matrix_type SdnR = {0,0}, SdnK = {0,0};  /* spin-down S (different zetas) */
static real *SdR2 = 0, *SdI2 = 0;  /* dense spin-down S(k), n x n */

int soc_active(void) { return soc_file_name != 0; }
int soc_num_states(void) { return soc_dim; }

/********
  Lx, Ly, Lz (hbar = 1) in the real YAeHMOP orbitals of shell l.
  real_a = sum_m U[a][m] Y_{l,m}, m = -l..l, and L_real = conj(U) L U^T.
********/
static void real_orbital_U(int l, zc_type U[7][7])
{
  /* (m, kind) of each real orbital, kind 0: m=0, 1: cos-like, 2: sin-like */
  static const int p_m[3] = {1,1,0}, p_k[3] = {1,2,0};
  static const int d_m[5] = {2,0,2,1,1}, d_k[5] = {1,0,2,1,2};
  static const int f_m[7] = {0,1,1,2,2,3,3}, f_k[7] = {0,1,2,2,1,1,2};
  const int *mm, *kk;
  int a, m, sgn, n = 2*l+1;
  real r2 = 1.0/sqrt(2.0);

  if( l == 1 ){ mm = p_m; kk = p_k; }
  else if( l == 2 ){ mm = d_m; kk = d_k; }
  else { mm = f_m; kk = f_k; }
  for(a=0;a<7;a++) for(m=0;m<7;m++){ U[a][m].re = 0.0; U[a][m].im = 0.0; }
  for(a=0;a<n;a++){
    m = mm[a];
    sgn = (m % 2) ? -1 : 1;
    if( kk[a] == 0 ){
      U[a][l].re = 1.0;
    } else if( kk[a] == 1 ){   /* (1/sqrt2) [ (-1)^m Y_{l,m} + Y_{l,-m} ] */
      U[a][l+m].re = sgn*r2;
      U[a][l-m].re = r2;
    } else {                   /* (i/sqrt2) [ Y_{l,-m} - (-1)^m Y_{l,m} ] */
      U[a][l-m].im = r2;
      U[a][l+m].im = -sgn*r2;
    }
  }
}

static void angular_momentum(int l, zc_type Lx[7][7], zc_type Ly[7][7], zc_type Lz[7][7])
{
  zc_type U[7][7], Lm[3][7][7];
  int n = 2*l+1, a, b, i, j, c;
  real cp;

  real_orbital_U(l,U);
  /* Lx, Ly, Lz in |l,m>, m = -l..l (index i = m+l) */
  for(c=0;c<3;c++) for(i=0;i<7;i++) for(j=0;j<7;j++){
    Lm[c][i][j].re = 0.0; Lm[c][i][j].im = 0.0;
  }
  for(i=0;i<n;i++) Lm[2][i][i].re = (real)(i-l);
  for(i=0;i<n-1;i++){
    real mm = (real)(i-l);
    cp = sqrt(l*(l+1) - mm*(mm+1));      /* <m+1| L+ |m> */
    /* Lx = (L+ + L-)/2 ; Ly = (L+ - L-)/(2i) */
    Lm[0][i+1][i].re += cp/2;  Lm[0][i][i+1].re += cp/2;
    Lm[1][i+1][i].im += -cp/2; Lm[1][i][i+1].im += cp/2;
  }
  /* L_real[a][b] = sum_ij conj(U[a][i]) L[i][j] U[b][j] */
  for(c=0;c<3;c++){
    zc_type (*out)[7] = (c == 0) ? Lx : (c == 1) ? Ly : Lz;
    for(a=0;a<n;a++) for(b=0;b<n;b++){
      double sr = 0.0, si = 0.0;
      for(i=0;i<n;i++) for(j=0;j<n;j++){
        double ur = U[a][i].re, ui = -U[a][i].im;          /* conj(U[a][i]) */
        double lr = Lm[c][i][j].re, li = Lm[c][i][j].im;
        double vr = U[b][j].re, vi = U[b][j].im;
        double tr = ur*lr - ui*li, ti = ur*li + ui*lr;
        sr += tr*vr - ti*vi;
        si += tr*vi + ti*vr;
      }
      out[a][b].re = sr; out[a][b].im = si;
    }
  }
}

static int same_symbol(const char *a, const char *b)
{
  while( *a == ' ' ) a++;
  while( *b == ' ' ) b++;
  while( *a && *b && *a != ' ' && *b != ' ' ){
    char ca = *a, cb = *b;
    if( ca >= 'a' && ca <= 'z' ) ca += 'A'-'a';
    if( cb >= 'a' && cb <= 'z' ) cb += 'A'-'a';
    if( ca != cb ) return 0;
    a++; b++;
  }
  return (!*a || *a == ' ') && (!*b || *b == ' ');
}

static int read_soc_file(const char *name, soc_entry_type *entries)
{
  FILE *f = fopen(name,"r");
  char line[MAX_STR_LEN], symb[80], shell[80], err[MAX_STR_LEN+100];
  int n = 0, lineno = 0;
  double lam;

  if( !f ){
    sprintf(err,"Can't open the SOC file %s",name);
    fatal(err);
  }
  while( fgets(line,MAX_STR_LEN,f) ){
    char *p = line;
    lineno++;
    while( *p == ' ' || *p == '\t' ) p++;
    if( !*p || *p == '\n' || *p == ';' || *p == '#' || *p == '!' ) continue;
    if( sscanf(p,"%79s %79s %lf",symb,shell,&lam) != 3 ){
      sprintf(err,"SOC file %s, line %d: expected 'Element shell lambda'",name,lineno);
      fatal(err);
    }
    if( n >= SOC_MAX_ENTRIES ) fatal("Too many lines in the SOC file.");
    if( strlen(symb) >= ATOM_SYMB_LEN ) fatal("Element symbol too long in the SOC file.");
    strcpy(entries[n].symb,symb);
    switch( shell[0] ){
    case 'p': case 'P': entries[n].l = 1; break;
    case 'd': case 'D': entries[n].l = 2; break;
    case 'f': case 'F': entries[n].l = 3; break;
    case 's': case 'S':
      fprintf(stderr,"SOC file: s shell of %s ignored (L = 0)\n",symb);
      continue;
    default:
      sprintf(err,"SOC file %s, line %d: shell must be p, d or f",name,lineno);
      fatal(err);
    }
    entries[n].lambda = lam;
    entries[n].used = 0;
    n++;
  }
  fclose(f);
  return n;
}

/****************************************************************************
 *
 *                   Procedure soc_setup_two_tables
 *
 *   SOC on top of a spin-polarized calculation: H_spinor = diag(H_up(k),
 *   H_down(k)) + H_SOC.  The main parameter table (already in cell->atoms)
 *   is the spin-up one; this reads the spin-down table (soc_parm_down) into
 *   dn_atoms.  Both tables must have the same shells; Hii and zetas may
 *   differ.  With different zetas the spin-down functions have their own
 *   overlaps, S_spinor = diag(S_up(k), S_down(k)) (functions of different
 *   spin are orthogonal whatever their radial parts), and the on-site
 *   lambda of the SOC file is used as it is (it is a parameter: the radial
 *   overlap of the up and down functions is absorbed in it).  Special atoms
 *   ('*', parameters in the input file) are the same in both channels.
 *
 ****************************************************************************/
static int zetas_differ(const atom_type *u, const atom_type *d)
{
  return fabs(u->exp_s - d->exp_s) > 1e-8 || fabs(u->exp_p - d->exp_p) > 1e-8 ||
    fabs(u->exp_p2 - d->exp_p2) > 1e-8 || fabs(u->coeff_p1 - d->coeff_p1) > 1e-8 ||
    fabs(u->coeff_p2 - d->coeff_p2) > 1e-8 ||
    fabs(u->exp_d - d->exp_d) > 1e-8 || fabs(u->exp_d2 - d->exp_d2) > 1e-8 ||
    fabs(u->coeff_d1 - d->coeff_d1) > 1e-8 || fabs(u->coeff_d2 - d->coeff_d2) > 1e-8 ||
    fabs(u->exp_f - d->exp_f) > 1e-8 || fabs(u->exp_f2 - d->exp_f2) > 1e-8 ||
    fabs(u->coeff_f1 - d->coeff_f1) > 1e-8 || fabs(u->coeff_f2 - d->coeff_f2) > 1e-8;
}

static void soc_setup_two_tables(cell_type *cell, detail_type *details, int num_orbs)
{
  atom_type *dn;
  int i, na = cell->num_atoms;
  char err[MAX_STR_LEN], *special;
  real input_const = details->the_const;

  if( details->Execution_Mode == THIN )
    fatal("SOC with two spin tables is not implemented in THIN mode.");
  if( details->do_muller_it )
    fatal("The Muller iteration is not implemented for SOC with two spin tables (it replaces\n"
          "\tthe parameters of the tables by spin-independent ones); use Charge Iteration.");

  dn = (atom_type *)calloc(na,sizeof(atom_type));
  special = (char *)calloc(na,sizeof(char));
  if( !dn || !special ) fatal("Can't allocate the spin-down parameters.");
  memcpy(dn,cell->atoms,na*sizeof(atom_type));
  /* special atoms keep the parameters of the input in both channels
     (fill_atomic_parms skips the atoms whose symbol starts with '&') */
  for(i=0;i<na;i++){
    int c;
    if( dn[i].symb[0] == '*' ){ special[i] = 1; continue; }
    /* atoms with parameters in the input file ('*' in the geometry) carry the
       symbol of their custom entry after read_inputfile */
    for(c=0;c<MAX_CUSTOM_ATOMS && custom_atoms[c].symb[0];c++){
      char a[ATOM_SYMB_LEN], b[ATOM_SYMB_LEN];
      strncpy(a,dn[i].symb,ATOM_SYMB_LEN); strncpy(b,custom_atoms[c].symb,ATOM_SYMB_LEN);
      a[ATOM_SYMB_LEN-1] = b[ATOM_SYMB_LEN-1] = 0;
      upcase(a); upcase(b);
      if( !strcmp(a,b) ){ special[i] = 1; break; }
    }
  }
  for(i=0;i<na;i++) if( special[i] ){ special[i] = dn[i].symb[0]; dn[i].symb[0] = '&'; }
  fill_atomic_parms(dn,na,(FILE *)0,soc_parm_down);
  for(i=0;i<na;i++) if( special[i] ) dn[i].symb[0] = special[i];
  for(i=0;i<na;i++){
    atom_type *u = &(cell->atoms[i]), *d = &(dn[i]);
    if( u->symb[0] == '&' || special[i] ) continue;
    if( u->ns != d->ns || u->np != d->np || u->nd != d->nd || u->nf != d->nf ){
      sprintf(err,"SOC with two spin tables: %s has different shells in the spin-down table.\n"
              "\tBoth channels must have the same orbitals (Hii and zetas may differ).",u->symb);
      fatal(err);
    }
    if( zetas_differ(u,d) ) soc_diff_zetas = 1;
  }
  free(special);
  dn_atoms = dn;
  up_prev = (atom_type *)calloc(na,sizeof(atom_type));
  if( !up_prev ) fatal("Can't allocate the spin-up parameters.");
  memcpy(up_prev,cell->atoms,na*sizeof(atom_type));

  if( soc_diff_zetas ){
    /* what still needs one set of functions per atom */
    if( details->num_FMO_frags || details->num_FCO_frags )
      fatal("FMO/FCO is not implemented for SOC with two spin tables with different zetas.");
    if( details->num_MOs_to_print )
      fatal("MO printing for viewkel is not implemented for SOC with two spin tables with\n"
            "\tdifferent zetas (the .MO file has one set of functions per atom).");
    if( details->Execution_Mode == MOLECULAR && details->eval_electrostat )
      fatal("The electrostatic energy is not implemented for SOC with two spin tables with\n"
            "\tdifferent zetas.");
  }

  if( soc_const_up > 0.0 ) details->the_const = soc_const_up;
  const_dn = soc_const_down > 0.0 ? soc_const_down : input_const;

  HdnR.dim = HdnK.dim = num_orbs;
  HdnR.mat = (real *)calloc(num_orbs*num_orbs,sizeof(real));
  HdnK.mat = (real *)calloc(num_orbs*num_orbs,sizeof(real));
  if( !HdnR.mat || !HdnK.mat ) fatal("Can't allocate the spin-down hamiltonian.");
  if( soc_diff_zetas ){
    SdnK.dim = num_orbs;
    SdnK.mat = (real *)calloc(num_orbs*num_orbs,sizeof(real));
    SdR2 = (real *)calloc(num_orbs*num_orbs,sizeof(real));
    SdI2 = (real *)calloc(num_orbs*num_orbs,sizeof(real));
    if( !SdnK.mat || !SdR2 || !SdI2 ) fatal("Can't allocate the spin-down overlap.");
  }
  soc_two_tables = 1;

  fprintf(output_file,"\n; SOC on top of two spin channels: H = diag(H_up, H_down) + H_SOC\n");
  fprintf(output_file,"#SOC_SPIN_TABLES: up (main table)  down %s\n",soc_parm_down);
  fprintf(output_file,"#SOC_SPIN_ZETAS: %s\n",soc_diff_zetas ? "different" : "same");
  if( soc_diff_zetas )
    fprintf(output_file,"; different zetas: S = diag(S_up, S_down), each channel with its own functions\n");
  fprintf(output_file,"#THE_CONST_UP: %lf\n#THE_CONST_DOWN: %lf\n",details->the_const,const_dn);
  fprintf(output_file,"; atom   Hii up (s p d f)   /   Hii down (s p d f)\n");
  for(i=0;i<na;i++){
    atom_type *u = &(cell->atoms[i]), *d = &(dn[i]);
    if( u->symb[0] == '&' ) continue;
    fprintf(output_file,"; %s%-3d %8.4lf %8.4lf %8.4lf %8.4lf / %8.4lf %8.4lf %8.4lf %8.4lf\n",
            u->symb,i+1,u->coul_s,u->coul_p,u->coul_d,u->coul_f,
            d->coul_s,d->coul_p,d->coul_d,d->coul_f);
    if( zetas_differ(u,d) )
      fprintf(output_file,";      zetas up (s p d)  %7.4lf %7.4lf %7.4lf  down %7.4lf %7.4lf %7.4lf\n",
              u->exp_s,u->exp_p,u->exp_d,d->exp_s,d->exp_p,d->exp_d);
  }
}

/* SOC with two tables and different zetas: the R space overlaps of the
   spin-down functions (COOP and average overlap populations), else 0 */
real *soc_down_overlap_R(void)
{
  return (soc_two_tables && soc_diff_zetas) ? SdnR.mat : (real *)0;
}
/* the overlaps must be kept in R space (memory.c) */
int soc_needs_R_overlaps(void)
{
  return soc_active() && soc_parm_down && soc_diff_zetas;
}

/* renormalize a contraction c1 chi(z1) + c2 chi(z2) of principal number n */
static void renorm_contraction(int n, real z1, real z2, real *c1, real *c2)
{
  real S12, norm;
  if( *c2 == 0.0 || z1 <= 0.0 || z2 <= 0.0 ) return;
  S12 = pow(4.0*z1*z2/((z1+z2)*(z1+z2)),(real)n+.5);
  norm = sqrt((*c1)*(*c1) + (*c2)*(*c2) + 2.0*S12*(*c1)*(*c2));
  if( norm > 0.0 ){ *c1 /= norm; *c2 /= norm; }
}

/****************************************************************************
 *
 *                   Procedure soc_follow_down
 *
 *   SOC with two spin tables: after charge iteration or zeta variation
 *   changed the spin-up parameters, give the spin-down ones the same
 *   changes (Hii and exponents), so the difference between the two tables
 *   (the exchange splitting) is kept.  Contractions are renormalized.
 *
 ****************************************************************************/
void soc_follow_down(cell_type *cell)
{
  int i;
  if( !soc_two_tables ) return;
  for(i=0;i<cell->num_atoms;i++){
    atom_type *u = &(cell->atoms[i]), *p = &(up_prev[i]), *d = &(dn_atoms[i]);
    if( u->symb[0] == '&' ) continue;
    d->coul_s += u->coul_s - p->coul_s;  d->coul_p += u->coul_p - p->coul_p;
    d->coul_d += u->coul_d - p->coul_d;  d->coul_f += u->coul_f - p->coul_f;
    d->exp_s += u->exp_s - p->exp_s;  d->exp_p += u->exp_p - p->exp_p;
    d->exp_d += u->exp_d - p->exp_d;  d->exp_f += u->exp_f - p->exp_f;
    if( d->coeff_p2 != 0.0 ){
      d->exp_p2 += u->exp_p2 - p->exp_p2;
      renorm_contraction(d->np,d->exp_p,d->exp_p2,&d->coeff_p1,&d->coeff_p2);
    }
    if( d->coeff_d2 != 0.0 ){
      d->exp_d2 += u->exp_d2 - p->exp_d2;
      renorm_contraction(d->nd,d->exp_d,d->exp_d2,&d->coeff_d1,&d->coeff_d2);
    }
    if( d->coeff_f2 != 0.0 ){
      d->exp_f2 += u->exp_f2 - p->exp_f2;
      renorm_contraction(d->nf,d->exp_f,d->exp_f2,&d->coeff_f1,&d->coeff_f2);
    }
    *p = *u;
  }
  fprintf(output_file,"; SOC spin-down table after the same changes: atom  Hii (s p d f)  zetas (s p d)\n");
  for(i=0;i<cell->num_atoms;i++){
    atom_type *d = &(dn_atoms[i]);
    if( cell->atoms[i].symb[0] == '&' ) continue;
    fprintf(output_file,";   %s%-3d %9.4lf %9.4lf %9.4lf %9.4lf   %7.4lf %7.4lf %7.4lf\n",
            cell->atoms[i].symb,i+1,d->coul_s,d->coul_p,d->coul_d,d->coul_f,
            d->exp_s,d->exp_p,d->exp_d);
  }
}

/****************************************************************************
 *
 *                   Procedure soc_build_down_hamiltonian
 *
 *   SOC with two spin tables: the R space hamiltonian of the spin-down
 *   channel (same call as for spin up, with the spin-down Hii and K).
 *   Called from driver.c right after the spin-up one is built.
 *
 ****************************************************************************/
/* exchange the EHT parameters (Hii, exponents, contraction coefficients)
   of two atoms; the rest (positions, symbols, ...) stays where it is */
static void swap_parms(atom_type *a, atom_type *b)
{
#define SWP(f) { real t_ = a->f; a->f = b->f; b->f = t_; }
  SWP(coul_s) SWP(coul_p) SWP(coul_d) SWP(coul_f)
  SWP(exp_s) SWP(exp_p) SWP(exp_d) SWP(exp_f) SWP(exp_d2) SWP(exp_f2)
  SWP(exp_p2) SWP(coeff_p1) SWP(coeff_p2) SWP(coeff_d1) SWP(coeff_d2) SWP(coeff_f1) SWP(coeff_f2)
#undef SWP
}

void soc_build_down_hamiltonian(cell_type *cell, detail_type *details,
                                hermetian_matrix_type overlapR, int num_orbs,
                                int *orbital_lookup_table, int tot_overlaps)
{
  int i;
  real save_const;
  if( !soc_two_tables ) return;
  /* the spin-down channel: its own parameters (and overlaps, with different
     zetas) through the usual routines; only the parameters are exchanged, so
     positions and the rest of the atoms stay current (Walsh steps, ...) */
  save_const = details->the_const;
  for(i=0;i<cell->num_atoms;i++)
    if( cell->atoms[i].symb[0] != '&' ) swap_parms(&(cell->atoms[i]),&(dn_atoms[i]));
  details->the_const = const_dn;
  if( soc_diff_zetas ){
    if( !SdnR.mat ){
      SdnR.dim = num_orbs;
      /* the size of Overlap_R (memory.c): one matrix for a molecule */
      SdnR.mat = (real *)calloc((size_t)num_orbs*num_orbs*
                                (details->Execution_Mode == MOLECULAR || tot_overlaps < 1 ? 1 : tot_overlaps),
                                sizeof(real));
      if( !SdnR.mat ) fatal("Can't allocate the spin-down overlap matrices.");
    }
    R_space_overlap_matrix(cell,details,SdnR,num_orbs,tot_overlaps,orbital_lookup_table,0);
    full_R_space_Hamiltonian(cell,details,SdnR,HdnR,num_orbs,orbital_lookup_table,0);
  } else {
    full_R_space_Hamiltonian(cell,details,overlapR,HdnR,num_orbs,orbital_lookup_table,0);
  }
  details->the_const = save_const;
  for(i=0;i<cell->num_atoms;i++)
    if( cell->atoms[i].symb[0] != '&' ) swap_parms(&(cell->atoms[i]),&(dn_atoms[i]));
}

/****************************************************************************
 *
 *                   Procedure soc_setup
 *
 *   Reads the SOC file, builds the dense on-site H_SOC (2n x 2n), allocates
 *   the spinor work space and refuses the options not done with SOC yet.
 *   Called once, after read_inputfile.
 *
 ****************************************************************************/
void soc_setup(cell_type *cell, detail_type *details, int num_orbs,
               int *orbital_lookup_table)
{
  soc_entry_type entries[SOC_MAX_ENTRIES];
  int num_entries, i, e, a, b, l, n2, off, lwork;
  zc_type Lx[7][7], Ly[7][7], Lz[7][7];
  char err[MAX_STR_LEN];

  if( !soc_active() ) return;
#ifndef USE_LAPACK
  fatal("SOC needs the program built with LAPACK (-DUSE_BLAS_LAPACK=ON).");
#endif

  /******** what is not done with SOC yet (see yaehmop-dev.md) ********/
#define SOC_REFUSE(cond,what) if( cond ){ \
    sprintf(err,"%s is not implemented with SOC yet.",what); fatal(err); }
#undef SOC_REFUSE

  if( soc_parm_down ) soc_setup_two_tables(cell,details,num_orbs);

  num_entries = read_soc_file(soc_file_name,entries);
  max_occup = 1.0;
  num_spin_comp = 2;
  soc_dim = n2 = 2*num_orbs;

  socR = (real *)calloc(n2*n2,sizeof(real));
  socI = (real *)calloc(n2*n2,sizeof(real));
  Hs = (real *)calloc(n2*n2,sizeof(real));
  Ss = (real *)calloc(n2*n2,sizeof(real));
  zH = (complex *)calloc(n2*n2,sizeof(complex));
  zS = (complex *)calloc(n2*n2,sizeof(complex));
  lwork = n2*n2;
  zwork = (complex *)calloc(lwork,sizeof(complex));
  soc_vals = (real *)calloc(n2,sizeof(real));
  soc_rwork = (real *)calloc(3*n2,sizeof(real));
  vupR = (real *)calloc(n2*num_orbs,sizeof(real));
  vupI = (real *)calloc(n2*num_orbs,sizeof(real));
  vdnR = (real *)calloc(n2*num_orbs,sizeof(real));
  vdnI = (real *)calloc(n2*num_orbs,sizeof(real));
  soc_chg = (real *)calloc(n2*num_orbs,sizeof(real));
  soc_occ = (real *)calloc(n2,sizeof(real));
  SdR = (real *)calloc(num_orbs*num_orbs,sizeof(real));
  SdI = (real *)calloc(num_orbs*num_orbs,sizeof(real));
  soc_tmp = (real *)calloc(n2*n2,sizeof(real));
  soc_spin = (real *)calloc(n2*num_orbs,sizeof(real));
  if( !socR || !socI || !Hs || !Ss || !zH || !zS || !zwork || !soc_vals || !soc_rwork ||
      !vupR || !vupI || !vdnR || !vdnI || !soc_chg || !soc_occ || !SdR || !SdI || !soc_tmp || !soc_spin )
    fatal("Can't allocate memory for the SOC matrices.");

  fprintf(output_file,"\n; Spin-orbit coupling (on-site lambda L.S), file %s\n",soc_file_name);
  fprintf(output_file,"#SOC_STATES: %d\n",n2);
  fprintf(output_file,"; spinor basis: orbitals 1..%d spin up, %d..%d the same orbitals spin down\n",
          num_orbs,num_orbs+1,n2);
  fprintf(output_file,"; atom   shell  lambda (eV)\n");

  for(i=0;i<cell->num_atoms;i++){
    atom_type *atom = &(cell->atoms[i]);
    int begin = orbital_lookup_table[i], shell_start[4];
    if( begin < 0 ) continue;              /* dummy atom */
    off = begin;
    shell_start[0] = shell_start[1] = shell_start[2] = shell_start[3] = -1;
    if( atom->ns ) off += 1;
    if( atom->np ){ shell_start[1] = off; off += 3; }
    if( atom->nd ){ shell_start[2] = off; off += 5; }
    if( atom->nf ){ shell_start[3] = off; off += 7; }
    for(e=0;e<num_entries;e++){
      if( !same_symbol(entries[e].symb,atom->symb) ) continue;
      l = entries[e].l;
      if( shell_start[l] < 0 ){
        sprintf(err,"SOC file: %s has no %c shell in this basis; ignored.",
                atom->symb,"spdf"[l]);
        if( !entries[e].used ) error(err);
        entries[e].used = 1;
        continue;
      }
      entries[e].used = 1;
      if( entries[e].lambda == 0.0 ) continue;
      fprintf(output_file,"%4s%-4d  %c   %10.6lf\n",atom->symb,i+1,"spdf"[l],
              entries[e].lambda);
      angular_momentum(l,Lx,Ly,Lz);
      {
        real h = 0.5*entries[e].lambda;
        int s0 = shell_start[l], nl = 2*l+1;
        for(a=0;a<nl;a++) for(b=0;b<nl;b++){
          int up_a = s0+a, up_b = s0+b;
          int dn_a = num_orbs+s0+a, dn_b = num_orbs+s0+b;
          /* up-up: +h Lz ; dn-dn: -h Lz */
          socR[up_a*n2+up_b] += h*Lz[a][b].re;  socI[up_a*n2+up_b] += h*Lz[a][b].im;
          socR[dn_a*n2+dn_b] -= h*Lz[a][b].re;  socI[dn_a*n2+dn_b] -= h*Lz[a][b].im;
          /* up-dn: h (Lx - i Ly) ; dn-up: h (Lx + i Ly) */
          socR[up_a*n2+dn_b] += h*(Lx[a][b].re + Ly[a][b].im);
          socI[up_a*n2+dn_b] += h*(Lx[a][b].im - Ly[a][b].re);
          socR[dn_a*n2+up_b] += h*(Lx[a][b].re - Ly[a][b].im);
          socI[dn_a*n2+up_b] += h*(Lx[a][b].im + Ly[a][b].re);
        }
      }
    }
  }
  for(e=0;e<num_entries;e++){
    if( !entries[e].used ){
      sprintf(err,"SOC file: element %s is not in the structure; ignored.",entries[e].symb);
      error(err);
    }
  }

  /* hermiticity check of H_SOC */
  for(a=0;a<n2;a++) for(b=0;b<n2;b++){
    if( fabs(socR[a*n2+b]-socR[b*n2+a]) > 1e-10 || fabs(socI[a*n2+b]+socI[b*n2+a]) > 1e-10 )
      FATAL_BUG("H_SOC is not hermitian.");
  }

  /* text copy of H_SOC, to rebuild the spinor H(R) outside the program */
  {
    char name[600];
    FILE *f;
    sprintf(name,"%s.SOC",details->filename);
    f = fopen(name,"w");
    if( !f ) fatal("Can't open the .SOC file.");
    fprintf(f,"# on-site spin-orbit matrix H_SOC (eV) from %s\n",soc_file_name);
    fprintf(f,"# spinor basis: 1..%d spin up, %d..%d spin down (YAeHMOP orbital order)\n",
            num_orbs,num_orbs+1,n2);
    if( soc_two_tables )
      fprintf(f,"# H_spinor(k) = diag(H_up(k), H_down(k)) + H_SOC,  S_spinor(k) = %s\n",
              soc_diff_zetas ? "diag(S_up(k), S_down(k))" : "I2 (x) S(k)");
    else
      fprintf(f,"# H_spinor(k) = I2 (x) H(k) + H_SOC,  S_spinor(k) = I2 (x) S(k)\n");
    fprintf(f,"# nonzero elements: i j Re Im (1-based)\n");
    fprintf(f,"%d\n",n2);
    for(a=0;a<n2;a++) for(b=0;b<n2;b++){
      if( fabs(socR[a*n2+b]) > 1e-12 || fabs(socI[a*n2+b]) > 1e-12 )
        fprintf(f,"%d %d %.12f %.12f\n",a+1,b+1,socR[a*n2+b],socI[a*n2+b]);
    }
    fclose(f);
  }

  if( details->Execution_Mode == MOLECULAR && details->use_symmetry ){
    fprintf(output_file,"; SOC: the symmetry labels of the MOs are not determined (double group)\n");
  }
}

/* physical element (i,j) of a packed hermitian matrix of dimension n */
static void packed_get(const real *mat, int n, int i, int j, real *re, real *im)
{
  if( i == j ){ *re = mat[i*n+i]; *im = 0.0; }
  else if( i < j ){ *re = mat[i*n+j]; *im = -mat[j*n+i]; }
  else { *re = mat[j*n+i]; *im = mat[i*n+j]; }
}

/* packed spinor H(k), S(k) (dimension 2n) from the spinless ones */
static void build_spinor(hermetian_matrix_type hamilK, hermetian_matrix_type overlapK,
                         int num_orbs, k_point_type *kpoint)
{
  int n2 = soc_dim, I, J;
  /* two spin tables: spin-down H(k) (and S(k), with different zetas) from
     the R space matrices of the spin-down channel */
  if( soc_two_tables ){
    hermetian_matrix_type sdn = overlapK;
    if( soc_diff_zetas ){
      if( details->Execution_Mode == MOLECULAR )
        memcpy(SdnK.mat,SdnR.mat,num_orbs*num_orbs*sizeof(real));
      else
        build_k_overlap_FAT(unit_cell,kpoint,SdnR,SdnK,num_orbs);
      sdn = SdnK;
    }
    if( details->Execution_Mode == MOLECULAR )
      memcpy(HdnK.mat,HdnR.mat,num_orbs*num_orbs*sizeof(real));
    else
      build_k_hamil_FAT(unit_cell,HdnR,HdnK,sdn,num_orbs);
  }
  for(I=0;I<n2;I++){
    for(J=I;J<n2;J++){
      int a = I/num_orbs, b = J/num_orbs, p = I%num_orbs, q = J%num_orbs;
      real hr = 0.0, hi = 0.0, sr = 0.0, si = 0.0;
      if( a == b ){
        if( a == 1 && soc_two_tables ) packed_get(HdnK.mat,num_orbs,p,q,&hr,&hi);
        else packed_get(hamilK.mat,num_orbs,p,q,&hr,&hi);
        if( a == 1 && soc_two_tables && soc_diff_zetas ) packed_get(SdnK.mat,num_orbs,p,q,&sr,&si);
        else packed_get(overlapK.mat,num_orbs,p,q,&sr,&si);
      }
      hr += socR[I*n2+J];
      hi += socI[I*n2+J];
      if( I == J ){
        Hs[I*n2+I] = hr;
        Ss[I*n2+I] = sr;
      } else {
        Hs[I*n2+J] = hr;  Hs[J*n2+I] = -hi;
        Ss[I*n2+J] = sr;  Ss[J*n2+I] = -si;
      }
    }
  }
}

/* eigenvalues of the packed spinor pair (Hs, Ss) into soc_vals; with
   want_vectors the eigenvectors are left in zH (column J = state J) */
static void diag_spinor(detail_type *details, int want_vectors)
{
#ifdef USE_LAPACK
  int n2 = soc_dim, j, k, itype = 1, info = 0, lwork = n2*n2;
  char jobz = want_vectors ? 'V' : 'N', uplo = 'L';
  for(j=0;j<n2;j++){
    int jtab = j*n2;
    for(k=j+1;k<n2;k++){
      int ktab = k*n2;
      zH[jtab+k].r = Hs[jtab+k];  zH[jtab+k].i = Hs[ktab+j];
      zS[jtab+k].r = Ss[jtab+k];  zS[jtab+k].i = Ss[ktab+j];
      zH[ktab+j].r = zH[ktab+j].i = 0.0;
      zS[ktab+j].r = zS[ktab+j].i = 0.0;
    }
    zH[jtab+j].r = Hs[jtab+j];  zH[jtab+j].i = 0.0;
    zS[jtab+j].r = Ss[jtab+j];  zS[jtab+j].i = 0.0;
  }
  if( !details->diag_wo_overlap ){
    zhegv((integer *)&itype,&jobz,&uplo,(integer *)&n2,zH,(integer *)&n2,zS,
          (integer *)&n2,soc_vals,zwork,(integer *)&lwork,soc_rwork,(integer *)&info);
  } else {
    zheev(&jobz,&uplo,(integer *)&n2,zH,(integer *)&n2,soc_vals,zwork,
          (integer *)&lwork,soc_rwork,(integer *)&info);
  }
  fprintf(status_file,"Error value from SOC Diagonalization (0 is good): %d\n",info);
  if( info != 0 ) error("Problems in the SOC diagonalization, try more overlaps.");
#endif
}

/* split the eigenvectors in zH into spin-up and spin-down parts */
static void extract_vectors(int num_orbs)
{
  int n2 = soc_dim, J, mu;
  for(J=0;J<n2;J++){
    complex *col = &zH[J*n2];
    for(mu=0;mu<num_orbs;mu++){
      vupR[J*num_orbs+mu] = col[mu].r;
      vupI[J*num_orbs+mu] = col[mu].i;
      vdnR[J*num_orbs+mu] = col[num_orbs+mu].r;
      vdnI[J*num_orbs+mu] = col[num_orbs+mu].i;
    }
  }
}

/* charge matrix of 2n spinor states with coefficients uR+i uI (spin up) and
   dR+i dI (spin down), n2 x n each, over the packed overlap Sp (n x n):
     chg[J][mu]  = sum over spins of Re( conj(c_mu) (S c)_mu )
     spin[J][mu] = the same, spin up minus spin down
   (S_spinor = I2 (x) S, so the spins don't mix); each row of chg adds up to 1. */
static void spinor_weights(const real *Sp, const real *Sp_dn, int num_orbs, const real *uR,
                           const real *uI, const real *dR, const real *dI, real *chg, real *spin)
{
  int n2 = soc_dim, J, mu, nu, sp;
  for(mu=0;mu<num_orbs;mu++) for(nu=0;nu<num_orbs;nu++)
    packed_get(Sp,num_orbs,mu,nu,&SdR[mu*num_orbs+nu],&SdI[mu*num_orbs+nu]);
  /* spin-down overlap (different zetas with two spin tables) */
  if( Sp_dn != Sp )
    for(mu=0;mu<num_orbs;mu++) for(nu=0;nu<num_orbs;nu++)
      packed_get(Sp_dn,num_orbs,mu,nu,&SdR2[mu*num_orbs+nu],&SdI2[mu*num_orbs+nu]);
  for(J=0;J<n2;J++){
    for(mu=0;mu<num_orbs;mu++) chg[J*num_orbs+mu] = 0.0;
    for(sp=0;sp<2;sp++){
      const real *vr = sp ? &dR[J*num_orbs] : &uR[J*num_orbs];
      const real *vi = sp ? &dI[J*num_orbs] : &uI[J*num_orbs];
      for(mu=0;mu<num_orbs;mu++){
        real scr = 0.0, sci = 0.0;
        const real *SR = (sp && Sp_dn != Sp) ? SdR2 : SdR, *SI = (sp && Sp_dn != Sp) ? SdI2 : SdI;
        for(nu=0;nu<num_orbs;nu++){
          real sr = SR[mu*num_orbs+nu], si = SI[mu*num_orbs+nu];
          scr += sr*vr[nu] - si*vi[nu];
          sci += sr*vi[nu] + si*vr[nu];
        }
        chg[J*num_orbs+mu] += vr[mu]*scr + vi[mu]*sci;
        if( sp == 0 ) spin[J*num_orbs+mu] = vr[mu]*scr + vi[mu]*sci;
        else spin[J*num_orbs+mu] -= vr[mu]*scr + vi[mu]*sci;
      }
    }
  }
}

static void spinor_charge_matrix(hermetian_matrix_type overlapK, int num_orbs)
{
  spinor_weights(overlapK.mat,(soc_two_tables && soc_diff_zetas) ? SdnK.mat : overlapK.mat,
                 num_orbs,vupR,vupI,vdnR,vdnI,soc_chg,soc_spin);
}

/* dense real and imaginary parts of a packed matrix of dimension n */
static void print_spinor_matrix(const char *what, real *packed, int n, detail_type *details)
{
  int i, j;
  real re, im;
  for(i=0;i<n;i++) for(j=0;j<n;j++){
    packed_get(packed,n,i,j,&re,&im);
    soc_tmp[i*n+j] = re;
  }
  fprintf(output_file,";\t\t --- %s, real part (spinor basis: 1..%d up, %d..%d down) ---\n",
          what,n/2,n/2+1,n);
  printmat(soc_tmp,n,n,output_file,1e-4,0,details->line_width);
  for(i=0;i<n;i++) for(j=0;j<n;j++){
    packed_get(packed,n,i,j,&re,&im);
    soc_tmp[i*n+j] = im;
  }
  fprintf(output_file,";\t\t --- %s, imaginary part ---\n",what);
  printmat(soc_tmp,n,n,output_file,1e-4,0,details->line_width);
}

/****************************************************************************
 *
 *                   Procedure soc_postprocess
 *
 *   What postprocess_results does at one k point (or for a molecule), for
 *   the spinor states: wave functions, energies and occupations, Mulliken
 *   overlap populations and net charges, charge matrix.
 *
 ****************************************************************************/
static void soc_postprocess(cell_type *cell, detail_type *details,
                            hermetian_matrix_type overlapK, int num_orbs)
{
  int n2 = soc_dim, i, j, J, a, begin, end;
  eigenset_type tmpset;
  real tot, accumR;

  if( details->wave_fn_PRT ){
    char transp = details->wave_fn_PRT & PRT_TRANSPOSE_FLAG;
    fprintf(output_file,";\t*** Spinor wavefunctions *** (states in rows, AO's in columns, "
            "energies INCREASE down)\n");
    fprintf(output_file,";\t***> SPIN UP, REAL:\n");
    print_labelled_mat(vupR,n2,num_orbs,output_file,1e-4,cell->atoms,cell->num_atoms,
                       orbital_lookup_table,num_orbs,transp,LABEL_COLS,details->line_width);
    fprintf(output_file,";\t***> SPIN UP, IMAGINARY:\n");
    print_labelled_mat(vupI,n2,num_orbs,output_file,1e-4,cell->atoms,cell->num_atoms,
                       orbital_lookup_table,num_orbs,transp,LABEL_COLS,details->line_width);
    fprintf(output_file,";\t***> SPIN DOWN, REAL:\n");
    print_labelled_mat(vdnR,n2,num_orbs,output_file,1e-4,cell->atoms,cell->num_atoms,
                       orbital_lookup_table,num_orbs,transp,LABEL_COLS,details->line_width);
    fprintf(output_file,";\t***> SPIN DOWN, IMAGINARY:\n");
    print_labelled_mat(vdnI,n2,num_orbs,output_file,1e-4,cell->atoms,cell->num_atoms,
                       orbital_lookup_table,num_orbs,transp,LABEL_COLS,details->line_width);
  }

  /* occupations of the 2n states (one electron each) */
  tmpset.dim = n2;
  tmpset.val = soc_vals;
  tmpset.vectR = tmpset.vectI = 0;
  for(J=0;J<n2;J++) soc_occ[J] = 0.0;
  calc_occupations(details,cell->num_electrons,n2,soc_occ,tmpset);

  if( details->levels_PRT || details->Execution_Mode == MOLECULAR ){
    real total = 0.0;
    fprintf(output_file,"\n#\t******* Energies (in eV)  and Occupation Numbers *******\n");
    fprintf(output_file,"; SOC: %d spinor states, one electron each\n",n2);
    for(J=0;J<n2;J++){
      fprintf(output_file,"%d:--->  %8.6lg  [%4.3lf Electrons]\n",J+1,soc_vals[J],soc_occ[J]);
      total += soc_occ[J]*soc_vals[J];
    }
    fprintf(output_file,"Total_Energy: %8.6lg\n",total);
    properties.total_E = total;
  }

  /* Mulliken analysis (same formulas as eval_mulliken, adding both spins) */
  /* (zeta variation takes the net charges of a molecule from here, as
     postprocess_results does without SOC) */
  if( details->OP_mat_PRT || details->ROP_mat_PRT || details->net_chg_PRT || details->vary_zeta ||
      details->mod_OP_mat_PRT || details->mod_ROP_mat_PRT || details->mod_net_chg_PRT ){
    real *OP = properties.OP_mat, *acc = work1;
    for(i=0;i<num_orbs;i++) acc[i] = 0.0;
    for(i=0;i<num_orbs;i++){
      for(j=0;j<num_orbs;j++){
        real sr, si;
        accumR = 0.0;
        if( soc_two_tables && soc_diff_zetas ){
          /* each spin part with the overlaps of its own functions */
          real au = 0.0, ad = 0.0, srd, sid;
          for(J=0;J<n2;J++){
            if( soc_occ[J] == 0.0 ) continue;
            au += soc_occ[J]*(vupR[J*num_orbs+i]*vupR[J*num_orbs+j] +
                              vupI[J*num_orbs+i]*vupI[J*num_orbs+j]);
            ad += soc_occ[J]*(vdnR[J*num_orbs+i]*vdnR[J*num_orbs+j] +
                              vdnI[J*num_orbs+i]*vdnI[J*num_orbs+j]);
          }
          packed_get(overlapK.mat,num_orbs,i,j,&sr,&si);
          packed_get(SdnK.mat,num_orbs,i,j,&srd,&sid);
          OP[i*num_orbs+j] = (i != j ? 2.0 : 1.0)*(au*sr + ad*srd);
          acc[i] += au*sr + ad*srd;
          continue;
        }
        for(J=0;J<n2;J++){
          if( soc_occ[J] == 0.0 ) continue;
          accumR += soc_occ[J]*(vupR[J*num_orbs+i]*vupR[J*num_orbs+j] +
                                vupI[J*num_orbs+i]*vupI[J*num_orbs+j] +
                                vdnR[J*num_orbs+i]*vdnR[J*num_orbs+j] +
                                vdnI[J*num_orbs+i]*vdnI[J*num_orbs+j]);
        }
        packed_get(overlapK.mat,num_orbs,i,j,&sr,&si);
        OP[i*num_orbs+j] = (i != j ? 2.0 : 1.0)*accumR*sr;
        acc[i] += accumR*sr;
      }
    }
    for(a=0;a<cell->num_atoms;a++){
      find_atoms_orbs(num_orbs,cell->num_atoms,a,orbital_lookup_table,&begin,&end);
      if( begin >= 0 ){
        real net = VALENCE_FRACTION*cell->atoms[a].num_valence;
        for(i=begin;i<end;i++) net -= acc[i];
        properties.net_chgs[a] = net;
      }
    }
    if( details->ROP_mat_PRT ){
      reduced_mulliken(cell->num_atoms,num_orbs,orbital_lookup_table,
                       properties.OP_mat,properties.ROP_mat);
    }
    /* modified Mulliken: a function of the (spin-summed) OP matrix only */
    if( details->mod_OP_mat_PRT || details->mod_ROP_mat_PRT || details->mod_net_chg_PRT ){
      modified_mulliken(cell,tmpset,overlapK,num_orbs,soc_occ,orbital_lookup_table,
                        properties.OP_mat,properties.mod_OP_mat,properties.mod_net_chgs,work1);
      if( details->mod_ROP_mat_PRT )
        reduced_mulliken(cell->num_atoms,num_orbs,orbital_lookup_table,
                         properties.mod_OP_mat,properties.mod_ROP_mat);
    }
  }
  if( details->OP_mat_PRT ){
    fprintf(output_file,
            "\n\n; \t\tq-q-q-q-q-q-q-q  Mulliken Analysis q-q-q-q-q-q-q-q\n");
    fprintf(output_file,"; SOC: populations add the spin-up and spin-down parts\n");
    fprintf(output_file,
            "\n;    Mulliken Overlap Population Matrix for %6.3lf electrons:\n",
            cell->num_electrons);
    print_labelled_mat(properties.OP_mat,num_orbs,num_orbs,output_file,1e-5,
                       cell->atoms,cell->num_atoms,orbital_lookup_table,
                       num_orbs,details->OP_mat_PRT & PRT_TRANSPOSE_FLAG,LABEL_BOTH,
                       details->line_width);
  }
  if( details->ROP_mat_PRT ){
    fprintf(output_file,
            "\n;   Reduced Mulliken Overlap Population Matrix for %6.3lf electrons:\n",
            cell->num_electrons);
    print_sym_mat(properties.ROP_mat,cell->num_atoms,cell->num_atoms,output_file,(char *)0,
                  (char *)0,details->line_width);
  }
  if( details->net_chg_PRT ){
    fprintf(output_file,"\n;    Net Atomic Charges for %6.3lf electrons:\n",
            cell->num_electrons);
    tot = 0.0;
    for(j=0;j<cell->num_atoms;j++){
      tot += properties.net_chgs[j];
      fprintf(output_file,"%d %s: %8.6lf\n", j+1, cell->atoms[j].symb,
              properties.net_chgs[j]);
    }
    fprintf(output_file,";      Total Charge is: %8.6lf\n",tot);

    /* spin moments (z): Mulliken spin-up minus spin-down populations */
    fprintf(output_file,"\n;    Spin moments m_z (N_up - N_down, Bohr magnetons) with SOC:\n");
    tot = 0.0;
    for(a=0;a<cell->num_atoms;a++){
      real m = 0.0;
      find_atoms_orbs(num_orbs,cell->num_atoms,a,orbital_lookup_table,&begin,&end);
      if( begin < 0 ) continue;
      for(J=0;J<n2;J++){
        if( soc_occ[J] == 0.0 ) continue;
        for(i=begin;i<end;i++) m += soc_occ[J]*soc_spin[J*num_orbs+i];
      }
      tot += m;
      fprintf(output_file,"%d %s: %8.6lf\n",a+1,cell->atoms[a].symb,m);
    }
    fprintf(output_file,";      Total moment is: %8.6lf\n",tot);
  }

  if( details->mod_OP_mat_PRT ){
    fprintf(output_file,
            "\n\n; \t\tq-q-q-q-q-q  Modified Mulliken Analysis q-q-q-q-q-q\n");
    fprintf(output_file,"; SOC: populations add the spin-up and spin-down parts\n");
    fprintf(output_file,
            "\n;    Modified Mulliken Overlap Population Matrix for %6.3lf electrons:\n",
            cell->num_electrons);
    print_labelled_mat(properties.mod_OP_mat,num_orbs,num_orbs,output_file,1e-5,
                       cell->atoms,cell->num_atoms,orbital_lookup_table,
                       num_orbs,details->mod_OP_mat_PRT & PRT_TRANSPOSE_FLAG,
                       LABEL_BOTH,details->line_width);
  }
  if( details->mod_ROP_mat_PRT ){
    fprintf(output_file,
            "\n;   Reduced Modified Mulliken Overlap Population Matrix for %6.3lf electrons:\n",
            cell->num_electrons);
    print_sym_mat(properties.mod_ROP_mat,cell->num_atoms,cell->num_atoms,
                  output_file,(char *)0,(char *)0,details->line_width);
  }
  if( details->mod_net_chg_PRT ){
    fprintf(output_file,"\n;   Modified Mulliken Net Atomic Charges for %6.3lf electrons:\n",
            cell->num_electrons);
    tot = 0.0;
    for(j=0;j<cell->num_atoms;j++){
      tot += properties.mod_net_chgs[j];
      fprintf(output_file,"%d %s: %8.6lf\n", j+1, cell->atoms[j].symb,
              properties.mod_net_chgs[j]);
    }
    fprintf(output_file,";      Total Charge is: %8.6lf\n",tot);
  }

  /* charge matrix: 2n spinor states x n orbitals */
  if( details->chg_mat_PRT || details->Rchg_mat_PRT ){
    fprintf(output_file,"\n\n; \t\tq-q-q-q-q-q-q-q  Charge Matrix q-q-q-q-q-q-q-q\n");
    fprintf(output_file,"; SOC: %d spinor states (rows), spin-up + spin-down weight of each AO\n",n2);
  }
  if( details->chg_mat_PRT ){
    fprintf(output_file,";  Complete Charge Matrix Independant of Occupation\n");
    print_labelled_mat(soc_chg,n2,num_orbs,output_file,1e-5,
                       cell->atoms,cell->num_atoms,orbital_lookup_table,
                       num_orbs,details->chg_mat_PRT & PRT_TRANSPOSE_FLAG,
                       LABEL_COLS,details->line_width);
  }
  if( details->Rchg_mat_PRT ){
    /* atoms in rows, states in columns, as reduced_charge_mat */
    for(a=0;a<cell->num_atoms;a++){
      find_atoms_orbs(num_orbs,cell->num_atoms,a,orbital_lookup_table,&begin,&end);
      for(J=0;J<n2;J++){
        real t = 0.0;
        if( begin >= 0 ) for(i=begin;i<end;i++) t += soc_chg[J*num_orbs+i];
        soc_tmp[a*n2+J] = t;
      }
    }
    fprintf(output_file,"\n;   Reduced Charge Matrix Independant of Occupation\n");
    printmat(soc_tmp,cell->num_atoms,n2,output_file,1e-5,(char)0,details->line_width);
  }
}

/****************************************************************************
 *
 *                   Procedure soc_fmo
 *
 *   FMO (molecules) and FCO (crystals) analysis of the spinor states.  The
 *   fragment orbitals are those of the spinless EHT hamiltonian (spin-up
 *   table with two tables), each holding a spin-up and a spin-down
 *   function, so the spin-up and spin-down parts of a spinor state are
 *   transformed separately with the usual transformation.  Populations
 *   and charge matrices add both spins, as the orbital projections do.  In
 *   the .FMO/.FCO files every fragment orbital appears twice (spin up, spin
 *   down) and the 2n x 2n charge matrix has the weight of each spin.
 *
 ****************************************************************************/
static int open_fco_file(const char *name)
{
  int fd;
#ifndef _MSC_VER
  fd = open(name, O_RDWR | O_TRUNC | O_CREAT, S_IRUSR | S_IWUSR);
#else
  fd = open(name, O_RDWR | O_TRUNC | O_CREAT, _S_IREAD | _S_IWRITE);
#endif
  if( fd == -1 ) fatal("Can't open .FCO file for binary I/O");
  return fd;
}
static real *fupR = 0, *fupI = 0, *fdnR = 0, *fdnI = 0, *fchg = 0, *fspin = 0;
static void soc_fmo(cell_type *cell, detail_type *details, hermetian_matrix_type overlapK,
                    hermetian_matrix_type hamilK, int num_orbs, int ik, int num_KPOINTS)
{
  static int fco_fd = -1;
  int n = num_orbs, n2 = soc_dim, sp, b, J, f, g, i;
  eigenset_type tmp;
  const char *what = details->num_FMO_frags ? "FMO" : "FCO";
  FMO_prop_type *P = details->FMO_props;

  if( !fupR ){
    fupR = (real *)calloc(6*n2*n,sizeof(real));
    if( !fupR ) fatal("Can't allocate the SOC FMO arrays.");
    fupI = fupR + n2*n; fdnR = fupI + n2*n; fdnI = fdnR + n2*n;
    fchg = fdnI + n2*n; fspin = fchg + n2*n;
  }

  /* FCO of a crystal: fragment orbitals at this k point (molecules: driver.c) */
  if( details->num_FCO_frags && details->Execution_Mode != MOLECULAR ){
    build_FMO_overlap(details,n,cell->num_atoms,overlapK,orbital_lookup_table);
    build_FMO_hamil(details,n,cell->num_atoms,hamilK,orbital_lookup_table);
    diagonalize_FMO(details,work1,work2,work3,cmplx_hamil,cmplx_overlap,cmplx_work);
    gen_FMO_tform_matrices(details);
  }

  /* coefficients in the fragment basis: two blocks of n states per spin part */
  tmp.dim = n;
  tmp.val = soc_vals;
  for(sp=0;sp<2;sp++) for(b=0;b<2;b++){
    tmp.vectR = (sp ? vdnR : vupR) + b*n*n;
    tmp.vectI = (sp ? vdnI : vupI) + b*n*n;
    tform_wavefuncs_to_FMO_basis(details,n,cell->num_atoms,tmp,orbital_lookup_table);
    memcpy((sp ? fdnR : fupR) + b*n*n,P->eigenset.vectR,n*n*sizeof(real));
    memcpy((sp ? fdnI : fupI) + b*n*n,P->eigenset.vectI,n*n*sizeof(real));
  }
  tform_hermetian_matrix_to_FMO_basis(details,n,cell->num_atoms,overlapK,work3,P->chg_mat,
                                      P->overlap,orbital_lookup_table);
  spinor_weights(P->overlap.mat,P->overlap.mat,n,fupR,fupI,fdnR,fdnI,fchg,fspin);

  if( details->wave_fn_PRT ){
    char transp = details->wave_fn_PRT & PRT_TRANSPOSE_FLAG;
    fprintf(output_file,"\n;\t Spinor wavefunctions in %s basis (spinor states in rows)\n",what);
    fprintf(output_file,";\t***> SPIN UP, REAL:\n");
    printmat(fupR,n2,n,output_file,1e-5,transp,details->line_width);
    fprintf(output_file,";\t***> SPIN UP, IMAGINARY:\n");
    printmat(fupI,n2,n,output_file,1e-5,transp,details->line_width);
    fprintf(output_file,";\t***> SPIN DOWN, REAL:\n");
    printmat(fdnR,n2,n,output_file,1e-5,transp,details->line_width);
    fprintf(output_file,";\t***> SPIN DOWN, IMAGINARY:\n");
    printmat(fdnI,n2,n,output_file,1e-5,transp,details->line_width);
    fprintf(output_file,"\n");
  }
  if( details->overlap_mat_PRT ){
    fprintf(output_file,";\t\t --- Overlap Matrix in %s basis %s ---\n",what,
            details->Execution_Mode == MOLECULAR ? "FMO_S(R)" : "FCO_S(K)");
    printmat(P->overlap.mat,n,n,output_file,1e-4,details->overlap_mat_PRT & PRT_TRANSPOSE_FLAG,
             details->line_width);
  }
  if( details->hamil_PRT ){
    tform_hermetian_matrix_to_FMO_basis(details,n,cell->num_atoms,hamilK,work3,P->chg_mat,
                                        P->hamil,orbital_lookup_table);
    fprintf(output_file,";\t\t --- Hamiltonian Matrix (without SOC) in %s basis ---\n",what);
    printmat(P->hamil.mat,n,n,output_file,1e-4,details->hamil_PRT & PRT_TRANSPOSE_FLAG,
             details->line_width);
  }

  /* Mulliken populations in the fragment basis, both spins:
     OP[f][g] = sum_J occ_J sum_spin Re( conj(c_f) c_g S_fg ) (x2 off the diagonal) */
  if( details->OP_mat_PRT || details->ROP_mat_PRT ){
    for(f=0;f<n;f++) for(g=0;g<n;g++){
      real sr, si, acc = 0.0;
      packed_get(P->overlap.mat,n,f,g,&sr,&si);
      for(J=0;J<n2;J++){
        if( soc_occ[J] == 0.0 ) continue;
        for(sp=0;sp<2;sp++){
          real ar = (sp ? fdnR : fupR)[J*n+f], ai = (sp ? fdnI : fupI)[J*n+f];
          real br = (sp ? fdnR : fupR)[J*n+g], bi = (sp ? fdnI : fupI)[J*n+g];
          acc += soc_occ[J]*((ar*br + ai*bi)*sr - (ar*bi - ai*br)*si);
        }
      }
      P->OP_mat[f*n+g] = (f != g ? 2.0 : 1.0)*acc;
    }
    if( details->ROP_mat_PRT )
      FMO_reduced_mulliken(details,cell->num_atoms,n,P->OP_mat,P->ROP_mat);
  }
  if( details->OP_mat_PRT ){
    fprintf(output_file,"\n\n; \tq-q-q-q-q-q-q-q  Mulliken Analysis in %s basis q-q-q-q-q-q-q-q\n",what);
    fprintf(output_file,"; SOC: populations add the spin-up and spin-down parts\n");
    fprintf(output_file,"\n;    %s Mulliken Overlap Population Matrix for %6.3lf electrons:\n",
            what,cell->num_electrons);
    printmat(P->OP_mat,n,n,output_file,1e-05,details->OP_mat_PRT & PRT_TRANSPOSE_FLAG,
             details->line_width);
  }
  if( details->ROP_mat_PRT ){
    fprintf(output_file,"\n;   Reduced %s Mulliken Overlap Population Matrix for %6.3lf electrons:\n",
            what,cell->num_electrons);
    print_sym_mat(P->ROP_mat,details->num_FMO_frags ? details->num_FMO_frags : details->num_FCO_frags,
                  details->num_FMO_frags ? details->num_FMO_frags : details->num_FCO_frags,
                  output_file,(char *)0,(char *)0,details->line_width);
  }
  if( details->chg_mat_PRT || details->Rchg_mat_PRT )
    fprintf(output_file,"\n\n# \tq-F-q-F-q-F-q-F  Charge Matrix in %s basis F-q-F-q-F-q-F-q\n",what);
  if( details->chg_mat_PRT ){
    fprintf(output_file,";  Complete %s Charge Matrix Independant of Occupation\n",what);
    fprintf(output_file,"; SOC: %d spinor states (rows), spin-up + spin-down weight of each fragment orbital\n",n2);
    printmat(fchg,n2,n,output_file,1e-05,details->chg_mat_PRT & PRT_TRANSPOSE_FLAG,
             details->line_width);
  }

  /* .FMO file (molecules): spinor energies with occupations and the charge
     matrix with the weight of each spin (columns f up, f down, ...) */
  if( details->num_FMO_frags && details->Execution_Mode == MOLECULAR ){
    fprintf(FMO_file,"; Molecular orbital energies (spinor states) and occupations\n");
    for(J=0;J<n2;J++) fprintf(FMO_file,"%lg %lg\n",soc_vals[J],soc_occ[J]);
    fprintf(FMO_file,"; charge matrix (columns: each fragment orbital spin up, spin down)\n");
    for(J=0;J<n2;J++){
      for(f=0;f<n;f++){
        real c = fchg[J*n+f], m = fspin[J*n+f];
        real up = 0.5*(c+m), dn = 0.5*(c-m);
        fprintf(FMO_file,"%6.4lg %6.4lg ",fabs(up) >= 1e-6 ? up : 0.0,fabs(dn) >= 1e-6 ? dn : 0.0);
      }
      fprintf(FMO_file,"\n");
    }
  }

  /* .FCO file (crystals), the binary layout read by fit_FCO, doubled */
  if( details->num_FCO_frags ){
    int nfo2;
    if( fco_fd < 0 ){
      char name[600];
      int test_int = 4231, nf = details->num_FCO_frags;
      real tot_K_weight = 0.0;
      sprintf(name,"%s.FCO",details->filename);
      fco_fd = open_fco_file(name);
      for(i=0;i<details->num_KPOINTS;i++) tot_K_weight += details->K_POINTS[i].weight;
      write(fco_fd,(const char *)&test_int,sizeof(int));
      write(fco_fd,(const char *)&(details->num_KPOINTS),sizeof(int));
      write(fco_fd,(const char *)&tot_K_weight,sizeof(real));
      write(fco_fd,(const char *)&n2,sizeof(int));
      write(fco_fd,(const char *)&nf,sizeof(int));
      for(i=0;i<nf;i++){
        nfo2 = 2*details->FMO_frags[i].num_orbs;
        write(fco_fd,(const char *)&nfo2,sizeof(int));
      }
    }
    write(fco_fd,(const char *)soc_vals,n2*sizeof(real));
    for(i=0;i<details->num_FCO_frags;i++){
      for(f=0;f<details->FMO_frags[i].num_orbs;f++){
        write(fco_fd,(const char *)&(details->FMO_frags[i].eigenset.val[f]),sizeof(real));
        write(fco_fd,(const char *)&(details->FMO_frags[i].eigenset.val[f]),sizeof(real));
      }
    }
    for(J=0;J<n2;J++){
      for(f=0;f<n;f++){
        real c = fchg[J*n+f], m = fspin[J*n+f], w[2];
        w[0] = 0.5*(c+m); w[1] = 0.5*(c-m);
        write(fco_fd,(const char *)w,2*sizeof(real));
      }
    }
    if( ik == num_KPOINTS-1 && details->walsh_details.num_steps <= 1 ){
      close(fco_fd); fco_fd = -1;
    }
  }

  /* keep the fragment coefficients and charge matrix for the average
     properties (projected DOS, occupations and COOPs on fragment orbitals) */
  if( details->avg_props && avg_prop_info ){
    avg_prop_info_type *api = &(avg_prop_info[ik]);
    for(i=0;i<n2*n;i++){
      api->FMO_orbs[i] = (float)fupR[i];
      api->FMO_orbsI[i] = (float)fupI[i];
      api->FMO_orbs_dn[i] = (float)fdnR[i];
      api->FMO_orbsI_dn[i] = (float)fdnI[i];
      api->FMO_chg_mat[i] = (float)fchg[i];
    }
  }
}

static int open_dump(const char *name, int nk, int n)
{
  int fd;
#ifndef _MSC_VER
  fd = open(name, O_RDWR | O_TRUNC | O_CREAT, S_IRUSR | S_IWUSR);
#else
  fd = open(name, O_RDWR | O_TRUNC | O_CREAT, _S_IREAD | _S_IWRITE);
#endif
  if( fd == -1 ) fatal("Can't open a SOC matrix dump file.");
  write(fd,(const char *)&nk,sizeof(int));
  write(fd,(const char *)&n,sizeof(int));
  return fd;
}

/****************************************************************************
 *
 *                   Procedure soc_kpoint
 *
 *   Work at one k point of loop_over_k_points (or for the molecule) when
 *   SOC is on: spinor matrices, printing and binary dumps (.HAM/.OV, same
 *   format as without SOC but with dimension 2n), diagonalization, the
 *   per-k properties (soc_postprocess) and the storage of the 2n states
 *   for the average properties (Fermi level, DOS, projected DOS, COOP,
 *   occupations and charges).
 *
 ****************************************************************************/
void soc_kpoint(cell_type *cell, detail_type *details, int ik, int num_KPOINTS,
                hermetian_matrix_type hamilK, hermetian_matrix_type overlapK,
                int num_orbs)
{
  char name[600];
  int j, mu, n2 = soc_dim, want_vectors;

  build_spinor(hamilK,overlapK,num_orbs,
               details->Execution_Mode == MOLECULAR ? (k_point_type *)0 : &(details->K_POINTS[ik]));

  if( details->overlap_mat_PRT )
    print_spinor_matrix(details->Execution_Mode == MOLECULAR ? "Spinor overlap S(R)" :
                        "Spinor overlap S(K)",Ss,n2,details);
  if( details->hamil_PRT )
    print_spinor_matrix(details->Execution_Mode == MOLECULAR ? "Spinor hamiltonian H(R)" :
                        "Spinor hamiltonian H(K)",Hs,n2,details);

  if( details->dump_overlap ){
    if( ik == 0 ){
      sprintf(name,"%s.OV",details->filename);
      soc_ov_file = open_dump(name,num_KPOINTS,soc_dim);
    }
    dump_hermetian_mat(soc_ov_file,Ss,soc_dim);
    if( ik == num_KPOINTS-1 ){ close(soc_ov_file); soc_ov_file = -1; }
  }
  if( details->dump_hamil ){
    if( ik == 0 ){
      sprintf(name,"%s.HAM",details->filename);
      soc_ham_file = open_dump(name,num_KPOINTS,soc_dim);
    }
    dump_hermetian_mat(soc_ham_file,Hs,soc_dim);
    if( ik == num_KPOINTS-1 ){ close(soc_ham_file); soc_ham_file = -1; }
  }
  if( details->dump_sparse_mats ){
    /* same text format as without SOC, spinor matrices of dimension 2n */
    static FILE *sp_ov = 0, *sp_ham = 0;
    if( ik == 0 ){
      sprintf(name,"%s.SPARSE.OV",details->filename);
      sp_ov = fopen(name,"w+");
      sprintf(name,"%s.SPARSE.HAM",details->filename);
      sp_ham = fopen(name,"w+");
      if( !sp_ov || !sp_ham ) fatal("Can't open the .SPARSE.OV/.SPARSE.HAM files.");
    }
    dump_sparse_mat(sp_ov,Ss,soc_dim,1e-08);
    dump_sparse_mat(sp_ham,Hs,soc_dim,1e-08);
    if( ik == num_KPOINTS-1 ){ fclose(sp_ov); fclose(sp_ham); sp_ov = sp_ham = 0; }
  }
  if( details->just_matrices ) return;

  want_vectors = !details->just_avgE;
  diag_spinor(details,want_vectors);

  if( want_vectors ){
    extract_vectors(num_orbs);
    spinor_charge_matrix(overlapK,num_orbs);
    soc_postprocess(cell,details,overlapK,num_orbs);
    if( details->num_FMO_frags || details->num_FCO_frags )
      soc_fmo(cell,details,overlapK,hamilK,num_orbs,ik,num_KPOINTS);
    if( details->num_MOs_to_print != 0 ){
      /* .MO file for viewkel: print_MOs takes the spinor states from soc.c */
      eigenset_type noset;
      memset(&noset,0,sizeof(noset));
      print_MOs(details,num_orbs,noset,ik,unique_atoms,num_unique_atoms,
                cell->num_atoms,orbital_lookup_table);
    }
  }

  /* keep the states for the average properties */
  if( details->avg_props && avg_prop_info ){
    avg_prop_info_type *api = &(avg_prop_info[ik]);
    for(j=0;j<n2;j++) api->energies[j] = (float)soc_vals[j];
    if( want_vectors ){
      for(j=0;j<n2*num_orbs;j++){
        api->orbs[j] = (float)vupR[j];
        api->orbsI[j] = (float)vupI[j];
        api->orbs_dn[j] = (float)vdnR[j];
        api->orbsI_dn[j] = (float)vdnI[j];
        api->chg_mat[j] = (float)soc_chg[j];
        api->chg_spin[j] = (float)soc_spin[j];
      }
    }
  }
  (void)mu;
}

/* SOC eigenvalues of one band structure k point (from the spinless H(k),S(k)) */
real *soc_band_energies(detail_type *details, hermetian_matrix_type hamilK,
                        hermetian_matrix_type overlapK, int num_orbs, k_point_type *kpoint)
{
  build_spinor(hamilK,overlapK,num_orbs,kpoint);
  diag_spinor(details,0);
  return soc_vals;
}

/* energy of spinor state J and charge matrix element (state J, orbital mu)
   at the last k point diagonalized (or the molecule), for the Walsh output */
real soc_state_energy(int J)
{
  if( J < 0 || J >= soc_dim ) fatal("Walsh: SOC state number out of range.");
  return soc_vals[J];
}
/* sum of occupation x energy of the spinor states of the last diagonalization
   (the molecule): the extended Hueckel energy used by eval_electrostatics */
real soc_occupied_energy(void)
{
  int J;
  real e = 0.0;
  for(J=0;J<soc_dim;J++) e += soc_occ[J]*soc_vals[J];
  return e;
}
/* coefficients of spinor state J (spin-up part ur+i ui, spin-down part dr+i di,
   n each) with the global phase fixed: the largest coefficient (over both
   spins) is real and positive.  Used for the .MO file of viewkel, which draws
   the real part. */
void soc_state_vector(int J, real *ur, real *ui, real *dr, real *di)
{
  int n = soc_dim/2, mu, best = 0;
  real big = -1.0, pr, pi, m;
  if( J < 0 || J >= soc_dim ) fatal("SOC state number out of range.");
  for(mu=0;mu<2*n;mu++){
    int k = J*n + (mu % n);
    real a = mu < n ? vupR[k] : vdnR[k], b = mu < n ? vupI[k] : vdnI[k];
    if( a*a + b*b > big + 1e-12 ){ big = a*a + b*b; best = mu; }
  }
  {
    int k = J*n + (best % n);
    real a = best < n ? vupR[k] : vdnR[k], b = best < n ? vupI[k] : vdnI[k];
    m = sqrt(a*a + b*b);
    pr = m > 0 ? a/m : 1.0;    /* multiply by conj(phase) */
    pi = m > 0 ? b/m : 0.0;
  }
  for(mu=0;mu<n;mu++){
    int k = J*n + mu;
    ur[mu] = vupR[k]*pr + vupI[k]*pi;  ui[mu] = vupI[k]*pr - vupR[k]*pi;
    dr[mu] = vdnR[k]*pr + vdnI[k]*pi;  di[mu] = vdnI[k]*pr - vdnR[k]*pi;
  }
}

/* reduced charge matrix element: weight of atom 'atom in spinor state J */
real soc_rchg_element(int atom, int J)
{
  int n = soc_dim/2, mu, begin, end;
  real t = 0.0;
  if( J < 0 || J >= soc_dim || atom < 0 || atom >= unit_cell->num_atoms )
    fatal("Walsh: SOC reduced charge matrix element out of range.");
  find_atoms_orbs(n,unit_cell->num_atoms,atom,orbital_lookup_table,&begin,&end);
  if( begin >= 0 ) for(mu=begin;mu<end;mu++) t += soc_chg[J*n+mu];
  return t;
}
/* modulus of the coefficient of AO mu in spinor state J, both spins:
   sqrt(|c_up|^2 + |c_down|^2).  Unlike the coefficient itself it does not
   depend on the basis LAPACK picks inside a Kramers pair. */
real soc_coeff_modulus(int J, int mu)
{
  int n = soc_dim/2, k = J*(soc_dim/2)+mu;
  if( J < 0 || J >= soc_dim || mu < 0 || mu >= n )
    fatal("Walsh: SOC coefficient out of range.");
  return sqrt(vupR[k]*vupR[k] + vupI[k]*vupI[k] + vdnR[k]*vdnR[k] + vdnI[k]*vdnI[k]);
}
real soc_chg_element(int J, int mu)
{
  int n = soc_dim/2;
  if( J < 0 || J >= soc_dim || mu < 0 || mu >= n )
    fatal("Walsh: SOC charge matrix element out of range.");
  return soc_chg[J*n+mu];
}

/****************************************************************************
 *
 *                   Procedure soc_print_spin_moments
 *
 *   Average spin moments m_z per atom with SOC (Mulliken populations of
 *   the spin-up minus the spin-down parts of the occupied states), with the
 *   same k point averaging as calc_avg_occups.
 *
 ****************************************************************************/
void soc_print_spin_moments(detail_type *details, cell_type *cell, int num_orbs,
                            K_orb_ptr_type *orbital_ordering, avg_prop_info_type *api)
{
  real *mom, tot_num_K = 0.0, num_occup_bands = 0.0, tot = 0.0, m;
  int i, j, a, begin, end;
  if( !soc_active() || details->just_avgE ) return;
  mom = (real *)calloc(num_orbs,sizeof(real));
  if( !mom ) fatal("Can't allocate the spin moments.");
  i = 0;
  while( orbital_ordering[i].occup > .0001 ){
    int k = orbital_ordering[i].Kpoint, MO = orbital_ordering[i].MO;
    float *w = &(api[k].chg_spin[MO*num_orbs]);
    for(j=0;j<num_orbs;j++)
      mom[j] += (real)w[j]*orbital_ordering[i].occup*details->K_POINTS[k].weight;
    i++;
  }
  for(i=0;i<details->num_KPOINTS;i++){
    tot_num_K += details->K_POINTS[i].weight*details->K_POINTS[i].num_filled_bands;
    num_occup_bands += details->K_POINTS[i].num_filled_bands;
  }
  fprintf(output_file,"\n# AVERAGE SPIN MOMENTS (SOC, m_z = N_up - N_down, Bohr magnetons)\n");
  for(a=0;a<cell->num_atoms;a++){
    find_atoms_orbs(num_orbs,cell->num_atoms,a,orbital_lookup_table,&begin,&end);
    if( begin < 0 ) continue;
    m = 0.0;
    for(j=begin;j<end;j++) m += mom[j];
    m *= num_occup_bands/(tot_num_K*details->num_KPOINTS);
    tot += m;
    fprintf(output_file,"%d %s %lf\n",a+1,cell->atoms[a].symb,m);
  }
  fprintf(output_file,"#TOTAL_SPIN_MOMENT: %lf\n",tot);
  free(mom);
}
