version 1.0

workflow minimap2 {
    input {
        String ref
        File fastqFile1
        File? fastqFile2
    }
    parameter_meta {
        ref: "path to the reference file used for alignment"
        fastqFile1: "a fastq file to be sequenced"
        fastqFile2: "an optional second fastq file to be sequenced"
    }

    meta {
        author: "Matthew Wong"
        email: "m2wong@oicr.on.ca"
        description: "Workflow to run align the fastq file to a reference genome"
        dependencies: [{
            name: "minimap2",
            url: "https://github.com/lh3/minimap2"
        }, {
            name: "samtools",
            url: "https://github.com/samtools/samtools"            
        }]
    }
    call convert2Sam {
        input:
            ref = ref,
            fastqFile1 = fastqFile1,
            fastqFile2 = fastqFile2
    }
    call sam2Bam {
        input:
            samfile = convert2Sam.alignment
    }
    output {
        File bam = sam2Bam.bam
        File bamIndex = sam2Bam.bamIndex
    }
}

task convert2Sam {
    input {
        String? minimap2 = "minimap2"
        String ref
        File fastqFile1
        File? fastqFile2
        String? modules = "minimap2/2.17 hg19/p13"
        Int? memory = 31
    }
    parameter_meta {
        minimap2: "minimap2 module name to use."
        ref: "the reference file name used for alignment"
        fastqFile1: "A fastq file to be sequenced"
        fastqFile2: "Optional second fastq file to be sequenced"
        modules: "Environment module names and version to load (space separated) before command execution."
        memory: "Memory (in GB) allocated for job."
    }
    meta {
        output_meta : {
            alignment: "output of minimap2 in .sam format"
        }
    }

    command <<<
        ~{minimap2} \
        -ax map-ont $HG19_ROOT/~{ref} \
        --MD \
        ~{fastqFile1} ~{fastqFile2} > alignment.sam
    >>>

    output {
        File alignment = "alignment.sam"
    }
    runtime {
        modules: "~{modules}"
        memory: "~{memory} G"
    }
}

task sam2Bam {
    input {
        File samfile
        String? samtools = "samtools"
        String? modules = "samtools/1.9"
        Int? memory = 31
    }
    parameter_meta {
        samtools: "samtools module name to use."
        samfile: "path to samfile"
        memory: "Memory (in GB) allocated for job."
    }
    meta {
        output_meta : {
            bam: "output of samtools in .bam format",
            bamIndex: "index of the bam file"
        }
    }

    command <<<
        ~{samtools} view -S -b ~{samfile} > alignment.bam
        ~{samtools} sort alignment.bam -o alignment.bam
        ~{samtools} index alignment.bam
    >>>

    output {
        File bam = "alignment.bam"
        File bamIndex = "alignment.bam.bai"
    }

    runtime {
        modules: "~{modules}"
        memory: "~{memory} G"
    }
}