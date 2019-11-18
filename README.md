# minimap2

Workflow to run align the fastq file to a reference genome

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
`fastqFile1`|File|a fastq file to be sequenced
`outputFileNamePrefix`|String|Variable used to set the name of the outputfile


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`fastqFile2`|File?|None|an optional second fastq file to be sequenced


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`convert2Sam.minimap2`|String?|"minimap2"|minimap2 module name to use.
`convert2Sam.modules`|String?|"minimap2/2.17"|Environment module names and version to load (space separated) before command execution.
`convert2Sam.memory`|Int?|31|Memory (in GB) allocated for job.
`sam2Bam.samtools`|String?|"samtools"|samtools module name to use.
`sam2Bam.modules`|String?|"samtools/1.9"|
`sam2Bam.memory`|Int?|31|Memory (in GB) allocated for job.


### Outputs

Output | Type | Description
---|---|---
`bam`|File|output of samtools in .bam format
`bamIndex`|File|index of the bam file


## Niassa + Cromwell

This WDL workflow is wrapped in a Niassa workflow (https://github.com/oicr-gsi/pipedev/tree/master/pipedev-niassa-cromwell-workflow) so that it can used with the Niassa metadata tracking system (https://github.com/oicr-gsi/niassa).

* Building
```
mvn clean install
```

* Testing
```
mvn clean verify -Djava_opts="-Xmx1g -XX:+UseG1GC -XX:+UseStringDeduplication" -DrunTestThreads=2 -DskipITs=false -DskipRunITs=false -DworkingDirectory=/path/to/tmp/ -DschedulingHost=niassa_oozie_host -DwebserviceUrl=http://niassa-url:8080 -DwebserviceUser=niassa_user -DwebservicePassword=niassa_user_password -Dcromwell-host=http://cromwell-url:8000
```

## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with wdl_doc_gen (https://github.com/oicr-gsi/wdl_doc_gen/)_
