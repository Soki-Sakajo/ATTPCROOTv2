#!/bin/bash

for i in {119..138}; do
	root -l -q "transfer_root.C($i)"
done
