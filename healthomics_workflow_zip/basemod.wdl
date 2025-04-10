version 1.0

task DSS_DMR {
    input {
        File tumor_bed
        File normal_bed
        String sname
        Int threads
    }

    Float file_size = ceil(size(tumor_bed, "GB") + size(normal_bed, "GB") + 10)

    command <<<
        set -euxo pipefail

        # Extract methylation info from bed files for DSS
        cut -f1,2,6,7 ~{tumor_bed} > ~{basename(tumor_bed) + ".tmp"}
        cut -f1,2,6,7 ~{normal_bed} > ~{basename(normal_bed) + ".tmp"}

        # Run DSS
        Rscript --vanilla /app/DSS_tumor_normal.R \
            ~{basename(tumor_bed) + ".tmp"} \
            ~{basename(normal_bed) + ".tmp"} \
            ~{sname + ".DMR.tsv"} \
            ~{threads}

        rm -f ~{basename(tumor_bed) + ".tmp"} ~{basename(normal_bed) + ".tmp"}
    >>>

    output {
        File output_DMR = sname + ".DMR.tsv"
    }

    runtime {
        docker: "860660336427.dkr.ecr.us-east-1.amazonaws.com/hifisomatic:somatic_r_tools"
        cpu: threads
        memory: "~{threads * 4} GB"
        disk: file_size + " GB"
        maxRetries: 2
        preemptible: 1
    }
}

task annotate_DMR {
    input {
        File DMR
        String sname
        Int threads
    }

    Float file_size = ceil(size(DMR, "GB"))

    command <<<
        set -euxo pipefail

        # Annotate DMRs
        Rscript --vanilla /app/annotatr_dmr.R \
            ~{DMR} \
            ~{sname} \
            ~{threads}
    >>>

    output {
        Array[File]+ output_DMR_annotated = glob(sname + "*.tsv.gz")
    }

    runtime {
        docker: "860660336427.dkr.ecr.us-east-1.amazonaws.com/hifisomatic:somatic_r_tools"
        cpu: threads
        memory: "~{threads * 4} GB"
        disk: file_size + " GB"
        maxRetries: 2
        preemptible: 1
    }
}
