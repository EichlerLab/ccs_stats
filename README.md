# PacBio CCS Stats Pipeline #

Computes summary statistics for unaligned PacBio CCS reads including.


## Run ##

Each run of this tool generates statistics for one sample, which may include one or more sequence files. To run stats
on several samples, invoke the tool once per sample.

### Input file ###

The only required input file is a list of sequence data files. Input files must be ".subreads.bam" files. RS II data
can be handled if the BAX are converted to BAMs first (see PacBio's bax2bam).

This pipeline is up to date as of Sequel 2.2 (1M chip). It is not sensitive to the sequencing platform or chemistry
as long as the structure of the input files remains consistent.

The file name should be "SAMPLE.fofn" where "SAMPLE" is the sample name. If it follows this convention, then the
pipeline can get the sample name from the FOFN file name instead of requiring that it be specified. Using sample names
allows several samples to be run in the same run directory. Alternatively, if the FOFN file does not contain the sample
name, then the sample name may be specified after `--config` as `sample=SAMPLE`.

This file name is provided to the pipeline after `--config` on the command line (see below).

### Execute ###

Define a variable that gives the full path to the pipeline code, which is the directory that contains `Snakefile`
and this `README.md` file. The pipeline itself does not use the variable, but commands in this README will.

Example:
`PIPELINE_DIR=/net/eichler/vol27/projects/structural_variation/nobackups/pipelines/subread_stats/202001`

Load required modules (may work with later versions of these modules):
```
module load miniconda/4.5.11
```

Run distributed (SGE):
`mkdir -p log; snakemake -s ${PIPELINE_DIR}/Snakefile --ri -j 30 -k --jobname "{rulename}.{jobid}" --drmaa " -V -cwd -e ./log -o ./log -pe serial {cluster.cpu} -l mfree={cluster.mem} -l h_rt={cluster.rt} -l gpfsstate=0 -j y -w n -S /bin/bash" -w 60 -u ${PIPELINE_DIR}/config/cluster.json --config fofn=/path/to/SAMPLE.fofn`

Run on a single node:
`snakemake -s ${PIPELINE_DIR}/Snakefile --ri -k --config fofn=/path/to/SAMPLE.fofn`

...where "SAMPLE.fofn" is the input file name.

#### Results ####

The results for each sample will be placed in a subdirectory named after the sample. Tables will appear as
tab-delimited files (".tab") and Excel files (".xlsx").

* cell_summary: Summary stats table per cell
* sample_summary: Summary stats table for all reads in the sample
* plot/cdf: Cumulative distribution plots of all subreads or the longest subread per read (ROI).
* cells/CELL-NAME/zmw_summary.tab.gz: Data each ZMW including number of subreads and subread lengths. These stats are
  used to compute the cell and sample summaries. If the raw input are not needed, these may be deleted after the stats
  are generated.


#### Specifying the sample name ####

The sample name may be explicitly provided on the command line by specifying `sample=SAMPLE` on the command line
anywhere after `--config`. If it is not specified, the FOFN file name is used to determine the sample name.

