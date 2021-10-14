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

// Put processes here.

workflow
{
}
