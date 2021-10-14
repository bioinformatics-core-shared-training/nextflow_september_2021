#!/bin/bash

export TMPDIR=temp
mkdir -p "$TMPDIR"

function clean_up
{
    rm -rf "$TMPDIR"
    exit $1
}

trap clean_up SIGHUP SIGINT SIGTERM

!{params.java} -Djava.io.tmpdir="$TMPDIR" \
-Xms!{task.memory.mega}m -Xmx!{task.memory.mega}m \
-jar !{params.picard} MarkDuplicates \
!{'INPUT=' + inBams.join(' INPUT=')} \
OUTPUT="!{outBam}" \
METRICS_FILE="!{metrics}" \
ASSUME_SORTED=true \
MAX_RECORDS_IN_RAM=1000000 \
CREATE_INDEX=true \
COMPRESSION_LEVEL=5 \
VALIDATION_STRINGENCY=SILENT \
TMP_DIR="$TMPDIR"

clean_up $?
