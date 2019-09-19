version 1.0

workflow minimap2 {
    input {
        String ref
        String fastqPath
    }
    call convert2BAM {
        input:
            ref = ref,
            fastqPath = fastqPath
    }
    output {
        File alignment = convert2BAM.alignment
    }
}

task convert2BAM {
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