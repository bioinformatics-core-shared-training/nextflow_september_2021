#!/bin/bash

set -o pipefail

!{params.bwa} \
    mem \
    -t !{task.cpus} \
    -o "!{outSam}" \
    "!{params.bwaIndex}" \
    "!{fastqFile}"
