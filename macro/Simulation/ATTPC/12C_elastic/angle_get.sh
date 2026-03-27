#!/bin/bash

for i in {0..99}; do
    root -l -q "angle_get.C($i,\"12c\")"
done
for i in {0..99}; do
    root -l -q "angle_get.C($i,\"p\")"
done

