version 1.0

workflow minimap2 {
    input {
        String ref
        File fastqFile
        String outputFileNamePrefix = basename( basename(fastqFile, ".gz"), ".fastq")
        String? additionalParameters
    }
    parameter_meta {
        ref: "the reference file name used for alignment"
        fastqFile: "a fastq file to be aligned"
        outputFileNamePrefix: "Variable used to set the name of the outputfile"
        additionalParameters: "Additional parameters to be added to the minimap2 command"
    }

    meta {
        author: "Matthew Wong"
        description: "Workflow to align the fastq file to a reference genome"
        dependencies: [{
            name: "minimap2",
            url: "https://github.com/lh3/minimap2"
        }, {
            name: "samtools",
            url: "https://github.com/samtools/samtools"            
        }]
        output_meta: {
          bam: {
            description: "Alignments, bam file.",
            vidarr_label: "bam"
          },
          bamIndex: {
            description: "Alignments, bai index",
            vidarr_label: "bamIndex"
          }
        }
    }
    call align {
        input:
            ref = ref,
            fastqFile = fastqFile,
            additionalParameters = additionalParameters
    }
    call sam2Bam {
        input:
            samfile = align.alignment,
            outputFileNamePrefix = outputFileNamePrefix,
    }
    output {
        File bam = sam2Bam.bam
        File bamIndex = sam2Bam.bamIndex
    }
}

task align {
    input {
        String minimap2 = "minimap2"
        String ref
        String? additionalParameters
        File fastqFile
        String modules = "minimap2/2.17"
        Int memory = 31
        Int timeout = 24
    }
    parameter_meta {
        minimap2: "minimap2 module name to use."
        ref: "the reference file name used for alignment"
        fastqFile: "a fastq file to be aligned"
        modules: "Environment module names and version to load (space separated) before command execution."
        memory: "Memory (in GB) allocated for job."
        timeout: "Runtime for the job in hours."
        additionalParameters: "Additional parameters to be added to the minimap2 command"
    }
    meta {
        output_meta : {
            alignment: "output of minimap2 in .sam format"
        }
    }

    command <<<
        ~{minimap2} \
        -ax map-ont ~{ref} \
        --MD \
        ~{additionalParameters} \
        ~{fastqFile} > alignment.sam
    >>>

    output {
        File alignment = "alignment.sam"
    }
    runtime {
        modules: "~{modules}"
        memory: "~{memory} G"
        timeout: "~{timeout}"
    }
}

task sam2Bam {
    input {
        File samfile
        String samtools = "samtools"
        String modules = "samtools/1.9"
        String outputFileNamePrefix
        Int memory = 31
        Int timeout = 24
    }
    parameter_meta {
        modules: "Modules to load for the task"
        samtools: "samtools module name to use."
        samfile: "path to samfile"
        memory: "Memory (in GB) allocated for job."
        timeout: "Runtime for the job in hours."
        outputFileNamePrefix: "Variable used to set the name of the outputfile"
    }
    meta {
        output_meta : {
            bam: "output of samtools in .bam format",
            bamIndex: "index of the bam file"
        }
    }

    command <<<
        ~{samtools} view -S -b ~{samfile} > alignment.bam
        ~{samtools} sort alignment.bam -o ~{outputFileNamePrefix}.bam
        ~{samtools} index ~{outputFileNamePrefix}.bam
    >>>

    output {
        File bam = "~{outputFileNamePrefix}.bam"
        File bamIndex = "~{outputFileNamePrefix}.bam.bai"
    }

    runtime {
        modules: "~{modules}"
        memory: "~{memory} G"
        timeout: "~{timeout}"
    }
}
