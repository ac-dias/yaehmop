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
 * Input file (--soc FILE), one line per element and shell, ';' or '#' start
 * comments:
 *      Mo  d  0.080
 *      S   p  0.050
 *******************************************************/
#include "bind.h"

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
static avg_prop_info_type *soc_avg = 0;
static K_orb_ptr_type *soc_ordering = 0;
static int soc_ham_file = -1, soc_ov_file = -1;

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
  if( details->Execution_Mode == MOLECULAR )
    fatal("SOC is only implemented for extended systems (not molecules) for now.");
  if( details->Execution_Mode == THIN )
    fatal("SOC is not implemented in THIN mode.");
#define SOC_REFUSE(cond,what) if( cond ){ \
    sprintf(err,"%s is not implemented with SOC yet.\n\tWith SOC: band structure, " \
            "energies, Fermi level, total DOS and the H/S dumps.",what); fatal(err); }
  SOC_REFUSE(details->num_proj_DOS,"Projected DOS");
  SOC_REFUSE(details->the_COOPS != 0,"COOP/COHP");
  SOC_REFUSE(details->num_FMO_frags || details->num_FCO_frags,"FMO/FCO");
  SOC_REFUSE(details->do_chg_it || details->do_muller_it,"Charge iteration");
  SOC_REFUSE(details->walsh_details.num_vars,"Walsh diagram");
  SOC_REFUSE(details->num_MOs_to_print,"MO printing");
  SOC_REFUSE(details->vary_zeta,"Zeta variation");
  SOC_REFUSE(details->hamil_PRT || details->overlap_mat_PRT,"Printing H or S");
  SOC_REFUSE(details->wave_fn_PRT,"Wave function printing");
  SOC_REFUSE(details->OP_mat_PRT || details->ROP_mat_PRT || details->avg_OP_mat_PRT ||
             details->avg_ROP_mat_PRT || details->mod_OP_mat_PRT || details->mod_ROP_mat_PRT,
             "Overlap population");
  SOC_REFUSE(details->chg_mat_PRT || details->Rchg_mat_PRT,"Charge matrix");
  SOC_REFUSE(details->net_chg_PRT || details->mod_net_chg_PRT,"Net charges");
  SOC_REFUSE(details->num_occup_AVG,"Alternate occupations");
  SOC_REFUSE(details->dump_sparse_mats,"Sparse matrix dump");
  SOC_REFUSE(details->do_moments,"Moments analysis");
#ifdef INCLUDE_NETCDF_SUPPORT
  SOC_REFUSE(details->do_netCDF,"netCDF output");
#endif
#undef SOC_REFUSE

  num_entries = read_soc_file(soc_file_name,entries);
  max_occup = 1.0;
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
  if( !socR || !socI || !Hs || !Ss || !zH || !zS || !zwork || !soc_vals || !soc_rwork )
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
    fprintf(f,"# H_spinor(k) = I2 (x) H(k) + H_SOC,  S_spinor(k) = I2 (x) S(k)\n");
    fprintf(f,"# nonzero elements: i j Re Im (1-based)\n");
    fprintf(f,"%d\n",n2);
    for(a=0;a<n2;a++) for(b=0;b<n2;b++){
      if( fabs(socR[a*n2+b]) > 1e-12 || fabs(socI[a*n2+b]) > 1e-12 )
        fprintf(f,"%d %d %.12f %.12f\n",a+1,b+1,socR[a*n2+b],socI[a*n2+b]);
    }
    fclose(f);
  }

  /* energies at every k point for the Fermi level and the total DOS */
  if( details->avg_props && !details->just_matrices ){
    soc_avg = (avg_prop_info_type *)calloc(details->num_KPOINTS,sizeof(avg_prop_info_type));
    soc_ordering = (K_orb_ptr_type *)calloc(n2*details->num_KPOINTS,sizeof(K_orb_ptr_type));
    if( !soc_avg || !soc_ordering ) fatal("Can't allocate the SOC energies.");
    for(i=0;i<details->num_KPOINTS;i++){
      soc_avg[i].energies = (float *)calloc(n2,sizeof(float));
      if( !soc_avg[i].energies ) fatal("Can't allocate the SOC energies.");
    }
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
                         int num_orbs)
{
  int n2 = soc_dim, I, J;
  for(I=0;I<n2;I++){
    for(J=I;J<n2;J++){
      int a = I/num_orbs, b = J/num_orbs, p = I%num_orbs, q = J%num_orbs;
      real hr = 0.0, hi = 0.0, sr = 0.0, si = 0.0;
      if( a == b ){
        packed_get(hamilK.mat,num_orbs,p,q,&hr,&hi);
        packed_get(overlapK.mat,num_orbs,p,q,&sr,&si);
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

/* eigenvalues of the packed spinor pair (Hs, Ss) into soc_vals */
static void diag_spinor(detail_type *details)
{
#ifdef USE_LAPACK
  int n2 = soc_dim, j, k, itype = 1, info = 0, lwork = n2*n2;
  char jobz = 'N', uplo = 'L';
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
 *   Work at one k point of loop_over_k_points when SOC is on: spinor
 *   matrices, binary dumps (.HAM/.OV, same format as without SOC but with
 *   dimension 2n), eigenvalues (to the .out if requested and stored for the
 *   Fermi level / DOS).
 *
 ****************************************************************************/
void soc_kpoint(cell_type *cell, detail_type *details, int ik, int num_KPOINTS,
                hermetian_matrix_type hamilK, hermetian_matrix_type overlapK,
                int num_orbs)
{
  char name[600];
  int j;

  build_spinor(hamilK,overlapK,num_orbs);

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
  if( details->just_matrices ) return;

  diag_spinor(details);

  if( details->energies_PRT || details->levels_PRT ){
    k_point_type *kp = &(details->K_POINTS[ik]);
    fprintf(output_file,"\n; SOC energies at k point %d (%lf %lf %lf)\n",ik+1,
            kp->loc.x,kp->loc.y,kp->loc.z);
    for(j=0;j<soc_dim;j++) fprintf(output_file,"%d:---> %10.6lf\n",j+1,soc_vals[j]);
  }
  if( soc_avg ){
    for(j=0;j<soc_dim;j++) soc_avg[ik].energies[j] = (float)soc_vals[j];
  }
}

/* SOC eigenvalues of one band structure k point (from the spinless H(k),S(k)) */
real *soc_band_energies(detail_type *details, hermetian_matrix_type hamilK,
                        hermetian_matrix_type overlapK, int num_orbs)
{
  build_spinor(hamilK,overlapK,num_orbs);
  diag_spinor(details);
  return soc_vals;
}

/****************************************************************************
 *
 *                   Procedure soc_avg_props
 *
 *   Fermi level (one electron per spinor state) and total DOS with SOC.
 *
 ****************************************************************************/
void soc_avg_props(cell_type *cell, detail_type *details, real *Fermi_E)
{
  if( !soc_avg ) return;
  sort_avg_prop_info(details,soc_dim,soc_avg,soc_ordering);
  find_crystal_occupations(details,cell->num_electrons,soc_dim,soc_ordering,Fermi_E);
  if( !details->no_total_DOS_PRT || !details->just_avgE ){
    gen_total_DOS(details,cell,soc_dim,soc_avg,soc_ordering);
    fprintf(output_file,"# END OF DOS\n\n");
  }
  fprintf(output_file,"\n;  The Fermi Level was determined (with SOC) for %d K points based on\n",
          details->num_KPOINTS);
  fprintf(output_file,";     an ordering of %d spinor crystal orbitals occupied by %lf electrons\n",
          soc_dim*details->num_KPOINTS,cell->num_electrons);
  fprintf(output_file,"#Fermi_Energy:  %lf\n",*Fermi_E);
}
