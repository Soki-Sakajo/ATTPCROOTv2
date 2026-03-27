#!/bin/bash

for i in {0..29}; do
    root -l -q "sim_macro_12c12c_0502_120matm.C($i)"
done
for i in {0..29}; do
    root -l -q "sim_macro_pp_0502_120matm.C($i)"
done

