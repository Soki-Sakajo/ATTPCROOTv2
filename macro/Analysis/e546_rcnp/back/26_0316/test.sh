#!/bin/bash

for pid in /proc/[0-9]*; do
  awk -v pid=$(basename "$pid") '
    /VmSwap/ && $2 > 0 { print $2 " KB\tPID=" pid }
  ' "$pid/status" 2>/dev/null
done | sort -nr | head
