#!/bin/bash

for i in {0..29}; do
    root -l -q "C12_12C12C_sim_80matm.C(100000,$i)"
done
for i in {0..29}; do
    root -l -q "C12_pp_sim_80matm.C(100000,$i)"
done

