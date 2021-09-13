# Nextflow training session

# Exercise Two - Alignment

The second exercise for learning Nextflow is a simplified version of a real
life pipeline: alignment. We are going to create a pipeline that aligns
single reads (not paired end), sorts them, marks duplicates, and calculates
alignment metrics on them. We'll use BWAmem for alignment and the Picard
tools for the other three steps.

You will need to refer to the Nextflow documentation quite a bit.
It is available at https://www.nextflow.io/docs/latest/

## Getting the exercise

We have a shell pipeline project to get you started, and four FASTQs for
you to work with:

* http://internal-bioinformatics.cruk.cam.ac.uk/training/nextflow/nextflow-alignment.tar.gz (the shell pipeline)
* http://internal-bioinformatics.cruk.cam.ac.uk/training/nextflow/nftraining-alignment-data.zip (the data)

You'll need to activate VPN outside the building to fetch
these, but you won't need the VPN after that.

Expand these archives into two different directories. You'll develop the
pipeline in the first directory and run it from the data directory.

## Before you start

The pipeline directory already contains the file `nextflow.config`. We'll
talk about this in the training, but for now you need to make a couple of tweaks.

There is a parameter that sets the location of java. Its default is the path
on Linux, so if you're on a Mac you'll need to change this to its location on
your computer (`which java` will give the path).

BWA and Picard are supplied in the pipeline zip file. Picard is Java so is
platform agnostic, but BWA is compiled code. There are two directories `bin_linux`
and `bin_mac`. Create a symbolic link to the appropriate one for your system
as `bin`; for example:

```BASH
ln -s bin_mac bin
```

We also need a reference genome. The exercise is set to use GRCh38, and the
defaults as supplied are the location of the reference on our cluster. If you
want to work offline on your laptop, you will need to copy the FASTA file and
the "`bwa-0.7.17`" directory from our references onto your local disk. If you
do this, change the "`referenceFasta`" and "`bwaIndex`" parameters in
`nextflow.config` too.

## Part One – A Simple Alignment Pipeline

For the first part of the exercise, we’ll run FASTQ as they are through the
align, sort, mark duplicates and calculate metrics steps.

Before you start, make a copy of the `start.nf` file as `simple.nf`. You'll
be editing `simple.nf`. To run Nextflow with your pipeline at each step, the
command is:

```
nextflow run <path to pipeline>/simple.nf
```

Run this from the data directory (wherever you have expanded the data zip).

One thing to note when using the templates. Declare the executable block in
the process using the "_shell_" executor rather than "_script_". So this block
will look like:

```
shell:
	template "<path>"
```

Have a little read through the section
https://www.nextflow.io/docs/latest/process.html#script

### Step 1: Alignment

1. Create a channel that will pick up all the FASTQ files in
"`params.fastqDir`" in the work flow.
    1. Use the function "`extractBaseName`" in the file to extract
    a "base name" from the file name. This is the file name without
    `.fq.gz` on the end.
    2. The output from this channel should be the base name and the
    FASTQ file.
2. Create a process to run _BWAmem_.
    1. Input needs to be a tuple of the base name and the FASTQ file.
    Output needs to be a tuple of the base name and the SAM file.
    2. Use the template `bwa/bwamem.sh` to run BWAmem. You'll need
    to call your FASTQ input "`fastqFile`" and out SAM output
    "`outSam`" in your process, or you can tweak the template to
    use your names.
3. In the work flow section, call your process passing it the FASTQ channel.

_Where does the work get done? How can you see what is actually being executed?_
_Where does the generated SAM file live?_

### Step 2: Sort

1. Create a process to run Picard's _SortSam_.
    1. Input needs to be a tuple of the base name and SAM/BAM file. Output
    needs to be a tuple of the base name and BAM file.
    2. Use the template `picard/SortSam.sh`. Your input needs to be named
    "`inBam`" and your output "`outBam`", or again tweak the template to
    use your names.
2. Connect the output of _BWAmem_ as the input to _SortSam_.

### Step 3: Mark Duplicates

1. Create a process to run Picard's _MarkDuplicates_.
    1. Input needs to be a tuple of the base name and BAM file. Output also
    needs to be a tuple of the base name and BAM file.
    2. Use the template `picard/MarkDuplicates.sh`. Your input needs to be
    named "`inBams`" (note the plural) and your output "`outBam`", or again
    tweak the template to use your names.
    3. _MarkDuplicates_ produces a metrics file as well as the BAM file.
    Create a second output for the metrics file called "`metrics`", or change
    the template to use your name.
2. Connect the output of _SortSam_ as the input to _MarkDuplicates_.

### Step 4: Alignment Metrics

1. Create a process to run Picard's _CollectAlignmentSummaryMetrics_.
    1. Input needs to be a tuple of the base name and BAM file. Output also
    needs to be a tuple of the base name and BAM file.
    2. Use the template `picard/ CollectAlignmentSummaryMetrics.sh`. Your
    input needs to be named "`inBam`", or tweak the template.
    3. The output from this process is a metrics file without a BAM file.
    Create an output for the metrics called "`metrics`".
2. Connect the output from _MarkDuplicates_ to _CollectAlignmentSummaryMetrics_.
    1. _MarkDuplicates_ has two outputs though: the basename with output BAM
    and the metrics file. How can we tell Nextflow which of these outputs to
    use as the input to alignment metrics?
    (_Hint: look in the "DSL2 - Process" section of the documentation._)

At this point, you have a pipeline that will do everything we need. However,
it's not clear where your final outputs are: they’ll be buried in the
work directory.

### Step 5: Publish outputs

1. Change _MarkDuplicates_ to publish its outputs into the directory specified by
"`params.bamDir`".
2. Change _CollectAlignmentSummaryMetrics_ to publish its outputs into the
same directory.

It seems that publishing outputs creates symbolic links that will break if the
work directory is removed.
_How can we change this to put proper files in the BAM directory?_

In real life the BAM files will be large, and creating a copy will be slow and
use extra space.
_Is there an alternative to copying them that will be quick and still use the_
_same space, while not being breakable symbolic links?_
_(Clue: this is a Unix mechanism that Nextflow allows us to use.)_

### Step 6: Specify process requirements

Nextflow uses the process' memory and CPU resources to know how many jobs to
run at once (in local mode). We’ve not set these.

1. Set all four tasks to use 4GB.
2. Set the three Picard tasks to use one CPU.
3. Set _BWAmem_ to use four CPUs.

We might want to also consider running this pipeline on the cluster, where
there is wall time.

1. Set all processes to have a four hour wall time.

## Part Two – Optimising for Parallel Tasks

The real alignment pipeline makes use of the cluster for BWA & BWAmem by
splitting the FASTQ files into chunks of reads, aligning each chunk, then
merging back into a single BAM file for each FASTQ file. This exercise adds
this feature to the pipeline you've written in the previous part as an example
of splitting and merging.

Before making changes, make a copy of your `single.nf` file as `chunked.nf`.
Make the changes to your pipeline in `chunked.nf`. The command to launch the
pipeline thus becomes:

```BASH
nextflow run <path to pipeline>/chunked.nf
```

### Step 1: Split the FASTQ file

1. Change the initial FASTQ channel that finds the files to split the FASTQ
file into chunks of 10,000 reads.
    1. _Hint: look at the Nextflow documentation under "Splitting Operators"._
    2. Your channel should now be tuples of three elements: the file's base
    name, the chunk number, and the FASTQ file (one chunk).

### Step 2: Adapt processes to work with chunks

The _BWAmem_ and _SortSam_ processes will now work with the chunks of FASTQ
rather than the whole file.

1. Change the _BWAmem_ process to accept the chunk number as a value in its
input tuple and pass it on in its output tuple.
2. Make the same change to _SortSam_.

### Step 3: Change the channel to remove chunk numbers

1. Change the output channel from _SortSam_ to remove the chunk number.
2. Group this channel by base name.
    1. _Hint: see "Combining Operators" in the documentation._

### Step 4: Adapt MarkDuplicates to work accept all chunks for a file

_MarkDuplicates_ merges multiple input BAM files into one as part of its
processing if required. We will use this step to combine the chunks for each
source file into one BAM file. The process as written will accept more than
one input BAM file.

1. Change the work flow so the adapted channel is what is passed into
_MarkDuplicates_.
