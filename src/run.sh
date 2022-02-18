#!/usr/bin/env bash

mkdir -p ../trajs ../rst ../models ../imgs ../results

# catch args
pdb=$1
max_iteration=${2:-10}
patience=${3:-1}

# preliminary MD
python prepare.py $pdb

# iterative sampling
for iteration in $(seq 1 $max_iteration)
do
  # check converge
  python converge.py $pdb $patience
  exit_code=$?
  if [[ $exit_code != 0 ]]; then
    break
  fi
  python simulation.py $pdb $iteration
  python training.py $pdb $iteration
  python pick.py $pdb $iteration
done
