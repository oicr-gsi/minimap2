version 1.0

workflow minimap2 {
    input {
        String ref
        String fastqPath
    }
    call convert2Sam {
        input:
            ref = ref,
            fastqPath = fastqPath
    }
    call sam2Bam {
        input:
            samfile = convert2Sam.alignment
    }
    output {
        File bam = sam2Bam.bam
    }
}

task convert2Sam {
    input {
        String? minimap2 = "minimap2"
        String ref
        String fastqPath
        String? modules = "minimap2/2.17"
        Int? memory = 31
    }

    command <<<
        ~{minimap2} \
        -a ~{ref} \
        ~{fastqPath} > alignment.sam
    >>>

    output {
        File alignment = "./alignment.sam"
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

    command <<<
        ~{samtools} view -S -b ~{samfile} > alignment.bam
    >>>

    output {
        File bam = "./alignment.bam"
    }

    runtime {
        modules: "~{modules}"
        memory: "~{memory} G"
    }
}