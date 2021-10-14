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
-jar !{params.picard} SortSam \
INPUT="!{inBam}" \
OUTPUT="!{outBam}" \
SORT_ORDER=coordinate \
MAX_RECORDS_IN_RAM=1000000 \
CREATE_INDEX=false \
COMPRESSION_LEVEL=1 \
VALIDATION_STRINGENCY=SILENT \
TMP_DIR="$TMPDIR"

clean_up $?
