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
-jar !{params.picard} CollectAlignmentSummaryMetrics \
INPUT=!{inBam} \
OUTPUT="!{metrics}" \
REFERENCE_SEQUENCE="!{params.referenceFasta}" \
VALIDATION_STRINGENCY=SILENT \
TMP_DIR="$TMPDIR"

clean_up $?
