#!/bin/bash
# YAeHMOP -> WanTiBEXOS recipe for MoS2 with SOC (see README.md). Run it inside this folder.
set -e
source /etc/profile.d/modules.sh 2>/dev/null || true
module load yaehmop-dev wtb-teste          # adapt to your machine: yaehmop and wtb.x in the PATH
ulimit -s unlimited
export OMP_NUM_THREADS=1 KMP_STACKSIZE=500mb

# 1) YAeHMOP: bands + H(k), S(k) on the 11x11x1 mesh, with spin-orbit coupling
cd yaehmop
yaehmop --soc MoS2_soc.dat MoS2.bind MoS2_params.dat > yaehmop.log

# 2) conversion to the WanTiBEXOS format, with the position matrix (~35 s), into out-interface/
python3 ../../yaehmop2wtb.py MoS2.bind -o ../out-interface --rmatrix --params MoS2_params.dat

# 3) WanTiBEXOS band structure
cp ../out-interface/huckel-SOC.dat ../wtb/
cd ../wtb
mkdir -p out
wtb.x < input.dat > wtb.log
cat out/gap-kpath.dat

# 4) figure
python3 plot_bands.py
