version 1.0

# Use Severus to call SV
task Severus_sv {
  input {
    File tumor_bam
    File tumor_bam_index
    File normal_bam
    File normal_bam_index
    File trf_bed
    File? phased_vcf
    String sname
    Int threads
    Int min_supp_reads
  }

  Float file_size = ceil(size(tumor_bam, "GB") + size(normal_bam, "GB") + size(phased_vcf, "GB") + 10)

  command <<<
    set -euxo pipefail
    
    echo "Running Severus for ~{sname}"

    severus --version

    severus \
      --target-bam ~{tumor_bam} \
      --control-bam ~{normal_bam} \
      ~{"--phasing-vcf " + phased_vcf} \
      --out-dir ~{sname + "_severus"} \
      -t ~{threads} \
      --vntr-bed ~{trf_bed} \
      --min-support ~{min_supp_reads}
  >>>

  output {
    File output_vcf = sname + "_severus/somatic_SVs/severus_somatic" + ".vcf"
  }

  runtime {
    docker: "860660336427.dkr.ecr.us-east-1.amazonaws.com/hifisomatic:severus"
    cpu: threads
    memory: "~{threads * 4} GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  }
}

# SAVANA SV calling
task SAVANA_sv {
  input {
    File tumor_bam
    File tumor_bam_index
    File normal_bam
    File normal_bam_index
    File ref_fasta
    File ref_fasta_index
    File? phased_vcf
    String sname
    Int threads
    Int? min_supp_reads
    Float? min_af
    Int svlength = 50
  }

  command <<<
    set -euxo pipefail

    echo -e "chr1\nchr2\nchr3\nchr4\nchr5\nchr6\nchr7\nchr8\nchr9\nchr10\nchr11\nchr12\nchr13\nchr14\nchr15\nchr16\nchr17\nchr18\nchr19\nchr20\nchr21\nchr22\nchrX\nchrY" > contig.txt
    
    echo "Running SAVANA for ~{sname}"
    savana \
      --tumour ~{tumor_bam} \
      --normal ~{normal_bam} \
      --ref ~{ref_fasta} \
      --threads ~{threads} \
      --outdir ~{sname + "_savana"} \
      --contigs contig.txt \
      --pb \
      --single_bnd \
      --no_blacklist \
      --sample ~{sname} \
      --length ~{svlength} \
      ~{"--phased_vcf " + phased_vcf} \
      ~{"--min_support " + min_supp_reads} \
      ~{"--min_af " + min_af}

  >>>

  output {
    File savana_vcf = sname + "_savana/" + sname + ".classified.somatic.vcf"
    File cnv_segments = sname + "_savana/" + sname + "_read_counts_mnorm_log2r_segmented.tsv"
    File cnv_absolute_cn = sname + "_savana/" + sname + "_segmented_absolute_copy_number.tsv"
    File purity_ploidy = sname + "_savana/" + sname + "_fitted_purity_ploidy.tsv"
    File purity_ploidy_solutions = sname + "_savana/" + sname + "_ranked_solutions.tsv"
    Array[File] savana_output = [savana_vcf, cnv_segments, cnv_absolute_cn, purity_ploidy, purity_ploidy_solutions]

  }

  runtime {
    container: "860660336427.dkr.ecr.us-east-1.amazonaws.com/hifisomatic:savana"
    cpu: threads
    memory: "~{threads * 4} GB"
  }
}