version 1.0

# Optionally, call SNVs and indels with ClairS
task ClairS {
  input {
    String sname
    File tumor_bam
    File tumor_bam_index
    File normal_bam
    File normal_bam_index
    File ref_fasta
    File ref_fasta_index
    File contig
    String platform
    Int threads
  }

  Float file_size = ceil(size(tumor_bam, "GB") * 2 + size(normal_bam, "GB") * 2 + size(ref_fasta, "GB") + size(contig, "GB") + 20)
  
  command <<<
    set -euxo pipefail

    echo "Running tumor normal variant calling for sample ~{sname} using ~{tumor_bam} and ~{normal_bam}"

    /opt/bin/run_clairs --version
    
    # If using old version of Singularity (e.g. 2.6), the path for python doesn't work and ClairS
    # will fail with python not found error. In that case, uncomment this:
    # export PATH=/opt/conda/envs/clairs/bin:$PATH
    /opt/bin/run_clairs \
      --normal_bam_fn ~{normal_bam} \
      --tumor_bam_fn ~{tumor_bam} \
      --ref_fn ~{ref_fasta} \
      --threads ~{threads} \
      --platform ~{platform} \
      --output_dir ~{sname + "_clairs"} \
      --conda_prefix /opt/conda/envs/clairs \
      --enable_indel_calling \
      --bed_fn ~{contig} \
      --enable_clair3_germline_output 
    
    # Instead of using ClairS remove intermediate flag,
    # we remove BAMs in tmp directory, but keep everything else in case
    # we need to troubleshoot
    if [ -d ~{sname + "_clairs/tmp/clair3_output/phased_output"} ]
      then rm -f ~{sname + "_clairs/tmp/clair3_output/phased_output/*.bam"} 
    fi
  >>>

  output {
    File output_snv_vcf = sname + "_clairs/snv.vcf.gz"
    File output_indel_vcf = sname + "_clairs/indel.vcf.gz"
    File output_snv_vcf_index = sname + "_clairs/snv.vcf.gz.tbi"
    File output_indel_vcf_index = sname + "_clairs/indel.vcf.gz.tbi"
    File output_snv_vcf_germline_tumor = sname + "_clairs/clair3_tumor_germline_output.vcf.gz"
    File output_snv_vcf_germline_normal = sname + "_clairs/clair3_normal_germline_output.vcf.gz"
  }

  runtime {
    docker: "860660336427.dkr.ecr.us-east-1.amazonaws.com/hifisomatic:clair3"
    cpu: threads
    memory: "~{threads * 4} GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  }
}

# Collect SNPs and indel into a single file
task gather_ClairS {
  input {
    Array[File] snv_vcf
    Array[File] snv_vcf_index
    Array[File] indel_vcf
    Array[File] indel_vcf_index
    File ref_fasta_index
    String sname
    Int threads
    Int snv_qual
    Int indel_qual
  }

  Float file_size = ceil(size(snv_vcf, "GB") * 2 + size(indel_vcf, "GB") * 2 + size(ref_fasta_index, "GB") + 10)

  command <<<
  set -euxo pipefail

  echo "Merging somatic small variant for ~{sname}"

  bcftools --version

  # Put the path of files of snv and indel into a single file
  # All the empty files stuff below are there to handle
  # the case of test dataset where there's no variants but I
  # want the workflow to complete
  touch ~{sname + ".snv_indel_vcf.txt"}
  for file in $(echo ~{sep=' ' snv_vcf} ~{sep=' ' indel_vcf});\
    do \
    # If file is not empty, then add to filelist
    if [[ -s ${file} ]]; then echo -e "${file}" >> ~{sname + ".snv_indel_vcf.txt"}; fi; done

  # if snv_indel_vcf.txt is not empty, then index and concat
  if [[ -s ~{sname + ".snv_indel_vcf.txt"} ]];
    then
      # Sometimes the indices failed to be localized properly. 
      # Reindex here to make it more robust...
      count=0
      for file in $(cat ~{sname + ".snv_indel_vcf.txt"})
      do
        ln -s ${file} ${count}.vcf.gz
        tabix ${count}.vcf.gz
        count=$((count+1))
      done

      # Create dummy header from ref fasta index to sort contig
      echo -e "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" > header.vcf
      bcftools reheader -f ~{ref_fasta_index} header.vcf | bgzip -c > header.reheader.vcf.gz
      tabix header.reheader.vcf.gz

      ls header.reheader.vcf.gz > file_list.txt

      ls *.vcf.gz | grep -v "reheader" >> file_list.txt

      bcftools concat -a -f file_list.txt |\
        bcftools sort -Oz -o ~{sname + ".somatic_small_variants.vcf.gz"}

      # Filter QUAL for snv and indel separately
      bcftools filter -i '(type=="indel" && QUAL>~{indel_qual}) || (type=="snp" && QUAL>~{snv_qual})' ~{sname + ".somatic_small_variants.vcf.gz"} |\
        bcftools sort -Oz -o ~{sname + ".somatic_small_variants.PASS.vcf.gz"}

      tabix ~{sname + ".somatic_small_variants.vcf.gz"}
      tabix ~{sname + ".somatic_small_variants.PASS.vcf.gz"}
    else
      echo "No variant to concat. Is this a test data with no variant?"
      touch ~{sname + ".somatic_small_variants.PASS.vcf.gz"}
      touch ~{sname + ".somatic_small_variants.PASS.vcf.gz.tbi"}
  fi
  >>>

  output {
    File output_vcf = sname + ".somatic_small_variants.PASS.vcf.gz"
    File output_vcf_index = sname + ".somatic_small_variants.PASS.vcf.gz.tbi"
  }

  runtime {
    docker: "860660336427.dkr.ecr.us-east-1.amazonaws.com/hifisomatic:bcftools"
    cpu: threads
    memory: "~{threads * 4} GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  }
}

# Collect germline variants
# Collect SNPs and indel into a single file
task gather_ClairS_germline {
  input {
    Array[File] tumor_vcf
    Array[File] normal_vcf
    String sname
    Int threads
  }

  Float file_size = ceil(size(tumor_vcf, "GB") * 2 + size(normal_vcf, "GB") * 2 + 10)

  command <<<
    set -euxo pipefail

    echo "Merging germline small variant for ~{sname}"

    bcftools --version

    # Put the path of files of snv and indel into a single file
    # All the empty files stuff below are there to handle
    # the case of test dataset where there's no variants but I
    # want the workflow to complete
    touch ~{sname + ".tumor_germline_vcf.txt"}
    touch ~{sname + ".normal_germline_vcf.txt"}
    for file in $(echo ~{sep=' ' tumor_vcf});\
      do \
      # If file is not empty, then add to filelist
      if [[ -s ${file} ]]; then echo -e "${file}" >> ~{sname + ".tumor_germline_vcf.txt"}; fi; done

    for file in $(echo ~{sep=' ' normal_vcf});\
      do \
      # If file is not empty, then add to filelist
      if [[ -s ${file} ]]; then echo -e "${file}" >> ~{sname + ".normal_germline_vcf.txt"}; fi; done

    # if snv_indel_vcf.txt is not empty, then index and concat
    if [[ -s ~{sname + ".tumor_germline_vcf.txt"} ]];
      then
        # Sometimes the indices failed to be localized properly. 
        # Reindex here to make it more robust...
        count=0
        for file in $(cat ~{sname + ".tumor_germline_vcf.txt"})
        do
          ln -s ${file} ${count}.tumor.vcf.gz
          tabix ${count}.tumor.vcf.gz
          count=$((count+1))
        done

        ls *.tumor.vcf.gz > file_list_tumor.txt

        bcftools concat -a -f file_list_tumor.txt |\
          bcftools sort |\
          bcftools norm -d exact |\
          bcftools view -f PASS -Oz -o ~{sname + ".tumor.germline_small_variants.vcf.gz"}

        tabix ~{sname + ".tumor.germline_small_variants.vcf.gz"}
      else
        echo "No tumor germline variant to concat. Is this a test data with no variant?"
        touch ~{sname + ".tumor.germline_small_variants.vcf.gz"}
        touch ~{sname + ".tumor.germline_small_variants.vcf.gz.tbi"}
    fi

    if [[ -s ~{sname + ".normal_germline_vcf.txt"} ]];
      then
        # Sometimes the indices failed to be localized properly. 
        # Reindex here to make it more robust...
        count=0
        for file in $(cat ~{sname + ".normal_germline_vcf.txt"})
        do
          ln -s ${file} ${count}.normal.vcf.gz
          tabix ${count}.normal.vcf.gz
          count=$((count+1))
        done

        ls *.normal.vcf.gz > file_list_normal.txt

        bcftools concat -a -f file_list_normal.txt |\
          bcftools sort |\
          bcftools norm -d exact |\
          bcftools view -f PASS -Oz -o ~{sname + ".normal.germline_small_variants.vcf.gz"}

        tabix ~{sname + ".normal.germline_small_variants.vcf.gz"}
      else
        echo "No normal germline variant to concat. Is this a test data with no variant?"
        touch ~{sname + ".normal.germline_small_variants.vcf.gz"}
        touch ~{sname + ".normal.germline_small_variants.vcf.gz.tbi"}
    fi
  >>>

  output {
    File output_tumor_germline_vcf = sname + ".tumor.germline_small_variants.vcf.gz"
    File output_tumor_germline_vcf_index = sname + ".tumor.germline_small_variants.vcf.gz.tbi"
    File output_normal_germline_vcf = sname + ".normal.germline_small_variants.vcf.gz"
    File output_normal_germline_vcf_index = sname + ".normal.germline_small_variants.vcf.gz.tbi"
  }

  runtime {
    docker: "860660336427.dkr.ecr.us-east-1.amazonaws.com/hifisomatic:bcftools"
    cpu: threads
    memory: "~{threads * 4} GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  }
}
