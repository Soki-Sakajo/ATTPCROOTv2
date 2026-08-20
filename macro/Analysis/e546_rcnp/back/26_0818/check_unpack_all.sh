#!/bin/bash

check_log () {
  local name="$1"
  local logfile="$2"
  local total="$3"

  if [ ! -f "$logfile" ]; then
    echo "=== $name ==="
    echo "log file not found: $logfile"
    echo
    return
  fi

#  done_jobs=$(wc -l < "$logfile")
  done_jobs=$(awk 'NR>1' "$logfile" | wc -l)
#  success=$(awk '$7==0' "$logfile" | wc -l)
  success=$(awk 'NR>1 && $7==0' "$logfile" | wc -l)
#  failed=$(awk '$7!=0' "$logfile" | wc -l)
  failed=$(awk 'NR>1 && $7!=0' "$logfile" | wc -l)

  read min avg max <<< $(awk '
    {
      h=$4/3600
      if (NR==1 || h<min) min=h
      if (NR==1 || h>max) max=h
      sum+=h
    }
    END {
      if (NR>0)
        printf "%.2f %.2f %.2f", min, sum/NR, max
    }' "$logfile")

  echo "=== $name ==="
  echo "total jobs   : $total"
  echo "done jobs    : $done_jobs"
  echo "success      : $success"
  echo "failed       : $failed"
  if [ "$done_jobs" -gt 0 ]; then
    echo "elapsed (h)  : min=$min  avg=$avg  max=$max"
  fi
  echo
}

# setting
check_log "76matm"  "joblog_76matm.log"  106
check_log "129matm" "joblog_129matm.log" 64
