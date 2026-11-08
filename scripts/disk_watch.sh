#!/bin/bash
# Every 10 min: free space on /data2; if below 150 GB, append an alert line to logs/jobs.log (the 2026-10-20 disk-full took down 9 jobs).
while true; do f=$(df -BG /data2 | awk "NR==2{gsub(\"G\",\"\",\$4); print \$4}"); echo "$(date -u +%FT%TZ) /data2 free ${f}G" >> /data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/disk_watch.log; [ "$f" -lt 150 ] && echo "ALERT $(date): /data2 has only ${f}G free" >> /data2/projects/Chitin_Metabolism_Analysis/01.Genome_Assembly/new-version/logs/jobs.log; sleep 600; done
