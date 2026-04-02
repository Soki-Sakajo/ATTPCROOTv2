#!/bin/bash

for i in {28..87}; do
	root -l -q "test_one_track.C($i)"
done
