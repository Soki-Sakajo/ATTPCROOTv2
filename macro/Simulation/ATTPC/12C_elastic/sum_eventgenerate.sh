#!/bin/bash

for i in {30..99}; do
    root -l -q "C12_12C12C_sim_130matm.C(100000,$i)"
done
for i in {30..99}; do
    root -l -q "C12_pp_sim_130matm.C(100000,$i)"
done

