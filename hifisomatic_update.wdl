version 1.0

import "alignment.wdl" as alignment
import "clairs.wdl" as clairs
import "cnvkit.wdl" as cnvkit
import "common.wdl" as common
import "structural_variants.wdl" as structural_variants
import "phasing.wdl" as phasing
import "basemod.wdl" as basemod
import "annotation.wdl" as annotation
import "prioritization.wdl" as prioritization
import "clonality.wdl" as clonality
import "deepsomatic.wdl" as deepsomatic

struct Patient {
    String patient_names
    Array[File] normal_bams
    Array[File] tumor_bams

    String? sex
}

struct TumorOnlyPatient {
    String patient_names
    Array[File] tumor_bams
}

struct TumorOnlyCohort {
    Array[TumorOnlyPatient] patients
}

struct Cohort {
    Array[Patient] patients
}

struct IndexData {
  File data
  File data_index
}

workflow hifisomatic {
  input {
    Cohort cohort
    # Files are not aligned by default. If aligned, set skip_align to true
    Boolean skip_align = false    
    # Strip kinetics or not
    Boolean strip_kinetics = false
    File ref_fasta
    File ref_fasta_dict
    File ref_gff
    # FASTA index refers to the standard faidx. For mmi index use ref_fasta
    File ref_fasta_index
    # Can also define FASTA MMI for pbmm2
    File? ref_fasta_mmi
    File ref_bed
    File cnvkit_refflat
    Int cnvkit_threads = 8
    Int pbmm2_threads = 24
    Int samtools_threads = 8
    Int merge_bam_threads = 4
    Int tumor_pileup_mincov = 5
    Int normal_pileup_mincov = 5
    Int cpg_pileup_threads = 8
    Int dss_threads = 16
    # Call small variants with deepsomatic
    Boolean use_deepsomatic = true
    Boolean call_small_variants = true
    String clairs_platform = "hifi_revio_ssrs"
    Int clairs_threads = 16
    Int clairs_snv_qual = 2
    Int clairs_indel_qual = 11
    ## Deepsomatic threads per task
    Int deepsomatic_threads = 16
    # Mutational signature max delta for MutationalPattern fitting
    Float mutsig_max_delta = 0.004
    # SV-related 
    File trf_bed
    Int sv_threads = 8
    Int wakhan_threads = 16
    File control_vcf
    File control_vcf_index
    File? severus_pon_tsv
    # Minimum reads support for Severus
    Int severus_min_reads = 3
    # AnnotSV cache can be downloaded using install script from https://lbgi.fr/AnnotSV/.
    # After the database is downloaded, zip the folder and provide the path to the zip file.
    # E.g. by default this is $ANNOTSV/share/AnnotSV
    # If this is not specified in input, AnnotSV will not be run
    File? annotsv_cache
    Int annotsv_threads = 8
    # Default number of threads for misc tasks (4GB per thread assigned for almost all steps)
    Int def_threads = 2
    # Scatter small variants calling into equal chunk per chromosome to make use of multiple nodes. Default of 75 Mbp per chromosome (total of 42 chunks for hg38)
    Int chunk_size = 75000000
    # Calculate DMR?
    Boolean calculate_DMR = true
    # Annotate VCF. If vep_cache is not specified in JSON, VEP will not be run
    # VEP cache can be downloaded from https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_refseq_vep_110_GRCh38.tar.gz
    File? vep_cache
    Int vep_threads = 8
    # Annotate germline variants?
    Boolean annotate_germline = false
    # Run SAVANA sv call too?
    Boolean run_savana = false
    # Minimum number of CG to prioritize
    Int ncg_to_filter = 50
    # Amber, Cobalt and Purple
    # Higher value = smoother CNV, but lower resolution
    Int pcf_gamma_value = 1000
    File? ensembl_data_dir_tarball
    Int hmftools_threads = 8
    # Visualization script
    File? report_vis_script

    File? prephased_tumor_vcf
    File? prephased_tumor_vcf_index
    File? prephased_normal_vcf
    File? prephased_normal_vcf_index
  }

  scatter (individual in cohort.patients) {
    String patient = individual.patient_names
    Array[File] patient_tumor_bam_files = individual.tumor_bams
    Array[File] patient_normal_bam_files = individual.normal_bams
    File ref_to_use = select_first([ref_fasta_mmi, ref_fasta])

    call alignment.align_all_bams as align_tumor {
      input:
        patient = patient,
        bam_files = patient_tumor_bam_files,
        ref_fasta = ref_to_use,
        ref_fasta_index = ref_fasta_index,
        pbmm2_threads = pbmm2_threads,
        merge_bam_threads = merge_bam_threads,
        samtools_threads = samtools_threads,
        skip_align = skip_align,
        strip_kinetics = strip_kinetics,
        suffix = "tumor"
    }

    call alignment.align_all_bams as align_normal {
      input:
        patient = patient,
        bam_files = patient_normal_bam_files,
        ref_fasta = ref_to_use,
        ref_fasta_index = ref_fasta_index,
        pbmm2_threads = pbmm2_threads,
        merge_bam_threads = merge_bam_threads,
        samtools_threads = samtools_threads,
        skip_align = skip_align,
        strip_kinetics = strip_kinetics,
        suffix = "normal"
    }

    ... (existing calls continue here)

    if(size(select_first([prephased_tumor_vcf, run_deepsomatic.clair3_tumor_vcf, gather_ClairS_germline.output_tumor_germline_vcf])) > 0){
      call phasing.hiphase_with_somatic as phaseTumorBam {
        input:
          bam = align_tumor.bam_final,
          bam_index = align_tumor.bam_final_index,
          vcf = select_first([prephased_tumor_vcf, run_deepsomatic.clair3_tumor_vcf, gather_ClairS_germline.output_tumor_germline_vcf]),
          vcf_index = select_first([prephased_tumor_vcf_index, run_deepsomatic.clair3_tumor_vcf_tbi, gather_ClairS_germline.output_tumor_germline_vcf_index]),
          somatic_SNP_indel_vcf = select_first([run_deepsomatic.deepsomatic_vcf, gather_ClairS.output_vcf]),
          somatic_SNP_indel_vcf_index = select_first([run_deepsomatic.deepsomatic_vcf_tbi, gather_ClairS.output_vcf_index]),
          pname = patient,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          threads = samtools_threads
      }
    }

    if(size(select_first([prephased_normal_vcf, run_deepsomatic.clair3_normal_vcf, gather_ClairS_germline.output_normal_germline_vcf])) > 0) {
      call phasing.hiphase as phaseNormalBam {
        input:
          bam = align_normal.bam_final,
          bam_index = align_normal.bam_final_index,
          vcf = select_first([prephased_normal_vcf, run_deepsomatic.clair3_normal_vcf, gather_ClairS_germline.output_normal_germline_vcf]),
          vcf_index = select_first([prephased_normal_vcf_index, run_deepsomatic.clair3_normal_vcf_tbi, gather_ClairS_germline.output_normal_germline_vcf_index]),
          pname = patient,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          threads = samtools_threads
      }
    }

    ... (rest of the workflow continues unchanged)
  }
}
