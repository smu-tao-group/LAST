#!/usr/bin/env bash

mkdir -p trajs rst models imgs results

pdb=$1

python prepare.py $pdb

for i in {0..30}
do
  python simulation.py $pdb $i
  python training.py $pdb $i
  python pick.py $pdb $i
done
