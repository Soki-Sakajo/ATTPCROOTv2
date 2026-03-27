#!/bin/bash

#root -l -q "unpack_rcnp.C(70)"
#root -l -q "unpack_rcnp.C(72)"
for i in {7..112}; do
#for i in {7..176}; do
    root -l -q "unpack_rcnp.C($i)"
done

