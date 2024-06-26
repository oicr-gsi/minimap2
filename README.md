# minimap2

Workflow to align the fastq file to a reference genome

## Overview

## Dependencies

* [minimap2](https://github.com/lh3/minimap2)
* [samtools](https://github.com/samtools/samtools)


## Usage

### Cromwell
```
java -jar cromwell.jar run minimap2.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`ref`|String|the reference file name used for alignment
`fastqFile`|File|a fastq file to be aligned


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`outputFileNamePrefix`|String|basename(basename(fastqFile,".gz"),".fastq")|Variable used to set the name of the outputfile
`additionalParameters`|String?|None|Additional parameters to be added to the minimap2 command


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`align.minimap2`|String|"minimap2"|minimap2 module name to use.
`align.modules`|String|"minimap2/2.17"|Environment module names and version to load (space separated) before command execution.
`align.memory`|Int|31|Memory (in GB) allocated for job.
`align.timeout`|Int|24|Runtime for the job in hours.
`sam2Bam.samtools`|String|"samtools"|samtools module name to use.
`sam2Bam.modules`|String|"samtools/1.9"|Modules to load for the task
`sam2Bam.memory`|Int|31|Memory (in GB) allocated for job.
`sam2Bam.timeout`|Int|24|Runtime for the job in hours.


### Outputs

Output | Type | Description | Labels
---|---|---|---
`bam`|File|Alignments, bam file.|vidarr_label: bam
`bamIndex`|File|Alignments, bai index|vidarr_label: bamIndex


## Commands
 This section lists command(s) run by minimap2 workflow
 
 * Running minimap2
 
 ### Run minimap2
 
 ```
         ~{minimap2} \
         -ax map-ont ~{ref} \
         --MD \
         ~{additionalParameters} \
         ~{fastqFile} > alignment.sam
 ```
 ### Converting to bam, sorting and indexing
 
 ```
         ~{samtools} view -S -b ~{samfile} > alignment.bam
         ~{samtools} sort alignment.bam -o ~{outputFileNamePrefix}.bam
         ~{samtools} index ~{outputFileNamePrefix}.bam
 ```
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
