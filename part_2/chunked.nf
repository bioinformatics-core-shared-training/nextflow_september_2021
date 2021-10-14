#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/**
 * This is a little Groovy function to get the name of the file without the .fq.gz suffix.
 * Don't worry too much about what is says. It takes in a path from a channel and returns
 * the name of the file without the suffix and, where present, a chunk number.
 * So if you say:
 *
 * def nameAndChunk = extractBaseName(theFile)
 * 
 * you can access the name with "nameAndChunk.name" and the chunk number with "nameAndChunk.chunk".
 * There might not be a chunk (in which case, it will be null).
 */
def extractBaseName(fastqFile)
{
    def m = fastqFile.name =~ /(?i)^(.+?)(\.(\d+))?\.fq(\.gz)?$/
    assert m.matches() : "Cannot extract base name from ${fastqFile}"
    def info = new Expando()
    info.name = m[0][1]
    info.chunk = m[0][3] as Integer
    return info
}


process bwa_mem
{
    cpus 4
    memory 8.GB
    time = 1.hour
    
    input:
        tuple val(basename), val(chunk), path(fastqFile)

    output:
        tuple val(basename), val(chunk), path(outSam)

    shell:
        outSam = "${basename}.${chunk}.bwamem.sam"
        template "bwa/bwamem.sh"
}

process picard_sortsam
{
    cpus = 1
    memory = 4.GB
    time = 4.hour
   
    input:
        tuple val(basename), val(chunk), path(inBam)

    output:
        tuple val(basename), val(chunk), path(outBam)

    shell:
        outBam = "${basename}.${chunk}.sorted.bam"
        template "picard/SortSam.sh"
}

process picard_markduplicates
{
    cpus = 1
    memory = 8.GB
    time = 4.hour
   
    publishDir params.bamDir, mode: "link"

    input:
        tuple val(basename), path(inBams)

    output:
        tuple val(basename), path(outBam), emit: merged_bam
        path metrics

    shell:
        outBam = "${basename}.bam"
        metrics = "${basename}.duplication.txt"

        template "picard/MarkDuplicates.sh"
}

process picard_alignmentmetrics
{
    cpus = 1
    memory = 4.GB
    time = 4.hour
   
    publishDir params.bamDir, mode: "link"

    input:
        tuple val(basename), path(inBam)

    output:
        path metrics

    shell:
        metrics = "${basename}.alignment.txt"
        template "picard/CollectAlignmentSummaryMetrics.sh"
}


workflow
{
    fastq_channel = channel.fromPath("${params.fastqDir}/*.fq.gz")
        .splitFastq(file: true, by: 10000)
        .map
        {
            fastq ->
            nameChunkPair = extractBaseName(fastq)
            tuple nameChunkPair.name, nameChunkPair.chunk, fastq 
        }

    bwa_mem(fastq_channel) | picard_sortsam 
    
    grouped_chunks = picard_sortsam.out.map { name, chunk, file -> tuple name, file }.groupTuple()
    
    picard_markduplicates(grouped_chunks)
    
    picard_alignmentmetrics(picard_markduplicates.out.merged_bam)
}
