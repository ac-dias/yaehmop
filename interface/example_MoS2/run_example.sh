#!/bin/bash
# Full example: YAeHMOP (yaehmop-dev) -> yaehmop2wtb.py -> WanTiBEXOS (wtb.x) -> comparison.
# Run it inside interface/example_MoS2.
set -e
source /etc/profile.d/modules.sh 2>/dev/null || true
module load yaehmop-dev wtb-teste
ulimit -s unlimited
export OMP_NUM_THREADS=1 KMP_STACKSIZE=500mb
here=$(pwd)

# 1) YAeHMOP runs (with and without SOC): band path + H(k), S(k) on the 11x11x1 mesh
for case in soc nosoc; do
  mkdir -p yaehmop_$case && cp MoS2.bind yaehmop_$case/
  cd yaehmop_$case
  if [ $case = soc ]; then
    yaehmop --soc ../MoS2_soc.dat MoS2.bind ../MoS2_params.dat > /dev/null
  else
    yaehmop MoS2.bind ../MoS2_params.dat > /dev/null
  fi
  cd "$here"
done

# 2) conversion to the WanTiBEXOS format and 3) wtb.x band structure
for case in soc nosoc; do
  python3 ../yaehmop2wtb.py yaehmop_$case/MoS2.bind -o wtb_$case
  cd wtb_$case
  mkdir -p out
  kind=$(ls input-*.dat | sed 's/input-//; s/.dat//')
  wtb.x < input-$kind.dat > wtb.log
  cd "$here"
done

# 4) comparison of the bands and figure
python3 compare_bands.py
