#!/usr/bin/env bash
# Sample system memory every 2 s: epoch used_MB avail_MB swap_used_MB top_rss_MB top_cmd
# usage: memsample.sh <out.log> &   (kill when done)
out=$1
while true; do
  m=$(free -m | awk '/^Mem:/{u=$3;a=$7} /^Swap:/{s=$3} END{print u, a, s}')
  top=$(ps -eo rss=,args= --sort=-rss | head -1 | awk '{printf "%d %s", $1/1024, $2" "$3}')
  echo "$(date +%s) $m $top" >> "$out"
  sleep 2
done
