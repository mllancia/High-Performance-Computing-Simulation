#!/bin/bash
# verify.sh <serial_binary> <mpi_binary> <workdir>
# For each lattice: run serial, then MPI at several (ranks, threads), compare
#  - bulk.dat line by line (energy, <phi>, fluctuation)
#  - every field snapshot written with -S, byte for byte
CPU=$1; MPI=$2; W=$3
mkdir -p $W; cd $W
cases=("128 1 1" "128 2 1" "128 4 1" "128 8 1" "128 4 2" \
       "128x64 4 1" "64x256 8 1" "96x64 3 1" "256x32 4 1")
printf "%-8s %5s %4s | %-14s | %-18s | %s\n" lattice ranks thr "bulk.dat" "field snapshots" "max |dE| rel"
for c in "${cases[@]}"; do
  set -- $c; N=$1; P=$2; T=$3
  ref=ref_$N
  if [ ! -d $ref ]; then
    mkdir -p $ref && (cd $ref && OMP_NUM_THREADS=1 $CPU -N $N -S seq >/dev/null 2>&1)
  fi
  d=mpi_${N}_${P}_${T}; rm -rf $d; mkdir -p $d
  (cd $d && mpirun --allow-run-as-root --oversubscribe -np $P -x OMP_NUM_THREADS=$T $MPI -N $N -S seq >/dev/null 2>&1)
  if cmp -s $ref/bulk.dat $d/bulk.dat; then b="identical"; else b="DIFFERS"; fi
  nref=$(ls $ref/seq | wc -l); nbad=0
  for f in $ref/seq/seq.*; do cmp -s $f $d/seq/$(basename $f) || nbad=$((nbad+1)); done
  [ "$(ls $d/seq 2>/dev/null | wc -l)" = "$nref" ] || nbad=$((nbad+1000))
  if [ $nbad = 0 ]; then s="$nref/$nref identical"; else s="$nbad of $nref differ"; fi
  rel=$(paste $ref/bulk.dat $d/bulk.dat | awk '{r=($4-$8)/$4; if(r<0)r=-r; if(r>m)m=r} END{printf "%.1e", m}')
  printf "%-8s %5s %4s | %-14s | %-18s | %s\n" $N $P $T "$b" "$s" "$rel"
done
