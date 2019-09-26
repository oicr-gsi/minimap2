version 1.0

workflow minimap2 {
    input {
        File ref
        File fastqFile1
        File? fastqFile2
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
    }
}

task convert2Sam {
    input {
        String? minimap2 = "minimap2"
        File ref
        File fastqFile1
        File? fastqFile2
        String? modules = "minimap2/2.17"
        Int? memory = 31
    }

    command <<<
        ~{minimap2} \
        -a ~{ref} \
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

    command <<<
        ~{samtools} view -S -b ~{samfile} > alignment.bam
    >>>

    output {
        File bam = "alignment.bam"
    }

    runtime {
        modules: "~{modules}"
        memory: "~{memory} G"
    }
}