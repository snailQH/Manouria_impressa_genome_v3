#!/bin/bash
# Every 10 min: free space on <shared-storage>; if below 150 GB, append an alert line to logs/jobs.log (the 2026-10-20 disk-full took down 9 jobs).
while true; do f=$(df -BG <shared-storage> | awk "NR==2{gsub(\"G\",\"\",\$4); print \$4}"); echo "$(date -u +%FT%TZ) <shared-storage> free ${f}G" >> $PROJECT/01.Genome_Assembly/new-version/logs/disk_watch.log; [ "$f" -lt 150 ] && echo "ALERT $(date): <shared-storage> has only ${f}G free" >> $PROJECT/01.Genome_Assembly/new-version/logs/jobs.log; sleep 600; done
