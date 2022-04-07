#!/usr/bin/bash
sed -i 's/docker run --user root:root/singularity exec/g' arriba_calling_rule
sed -i 's/docker run --user root:root/singularity exec/g' arriba_drawing_rule
sed -i 's/docker run --user root:root/singularity exec/g' arriba_normalization_rule
sed -i 's/docker run --user root:root/singularity exec/g' bcftools_concat_cancer_rule
sed -i 's/docker run --user root:root/singularity exec/g' bcftools_concat_wes_rule
sed -i 's/docker run --user root:root/singularity exec/g' bcftools_concat_wgs_rule
sed -i 's/docker run --user root:root/singularity exec/g' bcftools_merge_wts_wes_rule
sed -i 's/docker run --user root:root/singularity exec/g' bcftools_reheader_cancer_rule
sed -i 's/docker run --user root:root/singularity exec/g' bedtools_wes_rule
sed -i 's/docker run --user root:root/singularity exec/g' bwa_mem_cancer_rule
sed -i 's/docker run --user root:root/singularity exec/g' bwa_mem_cancer_rule.bak
sed -i 's/docker run --user root:root/singularity exec/g' fastq_concat_rule
sed -i 's/docker run --user root:root/singularity exec/g' fastq_concat_rule
sed -i 's/docker run --user root:root/singularity exec/g' fusion_catcher_rule
sed -i 's/docker run --user root:root/singularity exec/g' fusion_normalisation_script_rule
sed -i 's/docker run --user root:root/singularity exec/g' gatk_applybqsr_wgs_cancer_rule
sed -i 's/docker run --user root:root/singularity exec/g' gatk_applybqsr_wgs_rule
sed -i 's/docker run --user root:root/singularity exec/g' gatk_applybqsr_wts_rule
sed -i 's/docker run --user root:root/singularity exec/g' gatk_base_recalibrator_wgs_cancer_rule
sed -i 's/docker run --user root:root/singularity exec/g' gatk_base_recalibrator_wgs_rule
sed -i 's/docker run --user root:root/singularity exec/g' gatk_base_recalibrator_wts_rule
sed -i 's/docker run --user root:root/singularity exec/g' gatk_filter_mutect_calls_wes_rule
sed -i 's/docker run --user root:root/singularity exec/g' gatk_filter_mutect_calls_wes_tumor_only_rule
sed -i 's/docker run --user root:root/singularity exec/g' gatk_haplotype_caller_wgs_cancer_rule
sed -i 's/docker run --user root:root/singularity exec/g' gatk_haplotype_caller_wts_rule
sed -i 's/docker run --user root:root/singularity exec/g' gatk_mergevcf_mutect2_rule
sed -i 's/docker run --user root:root/singularity exec/g' gatk_mutect2_wes_rule
sed -i 's/docker run --user root:root/singularity exec/g' gatk_mutect2_wes_tumor_only_rule
sed -i 's/docker run --user root:root/singularity exec/g' gatk_select_variants_mutect_wes_rule
sed -i 's/docker run --user root:root/singularity exec/g' msisensor_rule
sed -i 's/docker run --user root:root/singularity exec/g' picard_collect_insert_size_metrics_wes_rule
sed -i 's/docker run --user root:root/singularity exec/g' picard_collect_insert_size_metrics_wgs_cancer_rule
sed -i 's/docker run --user root:root/singularity exec/g' picard_collect_insert_size_metrics_wts_rule
sed -i 's/docker run --user root:root/singularity exec/g' picard_markdup_wgs_rule
sed -i 's/docker run --user root:root/singularity exec/g' picard_markdup_wts_rule
sed -i 's/docker run --user root:root/singularity exec/g' quality_wes_changement_rule
sed -i 's/docker run --user root:root/singularity exec/g' rsync_vcf_tmb_rule
sed -i 's/docker run --user root:root/singularity exec/g' sambamba_flagstat_wgs_rule
sed -i 's/docker run --user root:root/singularity exec/g' sambamba_flagstat_wts_rule
sed -i 's/docker run --user root:root/singularity exec/g' sambamba_index_bam_rule
sed -i 's/docker run --user root:root/singularity exec/g' sambamba_index_genome_rule
sed -i 's/docker run --user root:root/singularity exec/g' sambamba_index_wgs_markdup_rule
sed -i 's/docker run --user root:root/singularity exec/g' sambamba_index_wgs_rule
sed -i 's/docker run --user root:root/singularity exec/g' sambamba_index_wts_markdup_rule
sed -i 's/docker run --user root:root/singularity exec/g' sambamba_index_wts_rule
sed -i 's/docker run --user root:root/singularity exec/g' sambamba_merge_facet_rule
sed -i 's/docker run --user root:root/singularity exec/g' sambamba_merge_wgs_rule
sed -i 's/docker run --user root:root/singularity exec/g' sambamba_merge_wts_rule
sed -i 's/docker run --user root:root/singularity exec/g' sambamba_slice_wgs_rule
sed -i 's/docker run --user root:root/singularity exec/g' sambamba_slice_wts_rule
sed -i 's/docker run --user root:root/singularity exec/g' signature_mut_matrix_decomp_rule
sed -i 's/docker run --user root:root/singularity exec/g' signature_mut_matrix_ext_rule
sed -i 's/docker run --user root:root/singularity exec/g' signature_mut_matrix_gen_rule
sed -i 's/docker run --user root:root/singularity exec/g' snpeff_tmb_rule
sed -i 's/docker run --user root:root/singularity exec/g' snpeff_wts_wes_rule
sed -i 's/docker run --user root:root/singularity exec/g' snpsift_annotate_leaves_trio_wgs_rule
sed -i 's/docker run --user root:root/singularity exec/g' snpsift_annotate_TMB_rule
sed -i 's/docker run --user root:root/singularity exec/g' snpsift_annotate_wes_rule
sed -i 's/docker run --user root:root/singularity exec/g' snpsift_dbnsfp_trio_GRCh38_wgs_rule
sed -i 's/docker run --user root:root/singularity exec/g' snpsift_dbnsfp_trio_wgs_rule
sed -i 's/docker run --user root:root/singularity exec/g' snpsift_phastcons_trio_wgs_rule
sed -i 's/docker run --user root:root/singularity exec/g' star_arriba_rule
sed -i 's/docker run --user root:root/singularity exec/g' star_cancer_rule
sed -i 's/docker run --user root:root/singularity exec/g' star_fusion_inspector_rule
sed -i 's/docker run --user root:root/singularity exec/g' star_fusion_rule
sed -i 's/docker run --user root:root/singularity exec/g' tabix_cancer_vcf_wes_rule
sed -i 's/docker run --user root:root/singularity exec/g' tabix_cancer_vcf_wgs_rule
sed -i 's/docker run --user root:root/singularity exec/g' tabix_cancer_vcf_wts_rule
sed -i 's/docker run --user root:root/singularity exec/g' tabix_wes_rule
sed -i 's/docker run --user root:root/singularity exec/g' tmb_curie_rule
sed -i 's/docker run --user root:root/singularity exec/g' vt_decompose_cancer_database_rule
sed -i 's/docker run --user root:root/singularity exec/g' vt_decompose_cancer_tmb_rule
sed -i 's/docker run --user root:root/singularity exec/g' vt_decompose_cancer_wes_rule
sed -i 's/docker run --user root:root/singularity exec/g' vt_decompose_cancer_wgs_rule
sed -i 's/docker run --user root:root/singularity exec/g' vt_decompose_cancer_wts_rule
sed -i 's/docker run --user root:root/singularity exec/g' vt_normalize_uniq_cancer_database_rule
sed -i 's/docker run --user root:root/singularity exec/g' vt_normalize_uniq_cancer_tmb_rule
sed -i 's/docker run --user root:root/singularity exec/g' vt_normalize_uniq_cancer_wes_rule
sed -i 's/docker run --user root:root/singularity exec/g' vt_normalize_uniq_cancer_wgs_rule
sed -i 's/docker run --user root:root/singularity exec/g' vt_normalize_uniq_cancer_wts_rule
sed -i 's/-v/-B/g' arriba_calling_rule
sed -i 's/-v/-B/g' arriba_drawing_rule
sed -i 's/-v/-B/g' arriba_normalization_rule
sed -i 's/-v/-B/g' bcftools_concat_cancer_rule
sed -i 's/-v/-B/g' bcftools_concat_wes_rule
sed -i 's/-v/-B/g' bcftools_concat_wgs_rule
sed -i 's/-v/-B/g' bcftools_merge_wts_wes_rule
sed -i 's/-v/-B/g' bcftools_reheader_cancer_rule
sed -i 's/-v/-B/g' bedtools_wes_rule
sed -i 's/-v/-B/g' bwa_mem_cancer_rule
sed -i 's/-v/-B/g' bwa_mem_cancer_rule.bak
sed -i 's/-v/-B/g' fastq_concat_rule
sed -i 's/-v/-B/g' fastq_concat_rule
sed -i 's/-v/-B/g' fusion_catcher_rule
sed -i 's/-v/-B/g' fusion_normalisation_script_rule
sed -i 's/-v/-B/g' gatk_applybqsr_wgs_cancer_rule
sed -i 's/-v/-B/g' gatk_applybqsr_wgs_rule
sed -i 's/-v/-B/g' gatk_applybqsr_wts_rule
sed -i 's/-v/-B/g' gatk_base_recalibrator_wgs_cancer_rule
sed -i 's/-v/-B/g' gatk_base_recalibrator_wgs_rule
sed -i 's/-v/-B/g' gatk_base_recalibrator_wts_rule
sed -i 's/-v/-B/g' gatk_filter_mutect_calls_wes_rule
sed -i 's/-v/-B/g' gatk_filter_mutect_calls_wes_tumor_only_rule
sed -i 's/-v/-B/g' gatk_haplotype_caller_wgs_cancer_rule
sed -i 's/-v/-B/g' gatk_haplotype_caller_wts_rule
sed -i 's/-v/-B/g' gatk_mergevcf_mutect2_rule
sed -i 's/-v/-B/g' gatk_mutect2_wes_rule
sed -i 's/-v/-B/g' gatk_mutect2_wes_tumor_only_rule
sed -i 's/-v/-B/g' gatk_select_variants_mutect_wes_rule
sed -i 's/-v/-B/g' msisensor_rule
sed -i 's/-v/-B/g' picard_collect_insert_size_metrics_wes_rule
sed -i 's/-v/-B/g' picard_collect_insert_size_metrics_wgs_cancer_rule
sed -i 's/-v/-B/g' picard_collect_insert_size_metrics_wts_rule
sed -i 's/-v/-B/g' picard_markdup_wgs_rule
sed -i 's/-v/-B/g' picard_markdup_wts_rule
sed -i 's/-v/-B/g' quality_wes_changement_rule
sed -i 's/-v/-B/g' rsync_vcf_tmb_rule
sed -i 's/-v/-B/g' sambamba_flagstat_wgs_rule
sed -i 's/-v/-B/g' sambamba_flagstat_wts_rule
sed -i 's/-v/-B/g' sambamba_index_bam_rule
sed -i 's/-v/-B/g' sambamba_index_genome_rule
sed -i 's/-v/-B/g' sambamba_index_wgs_markdup_rule
sed -i 's/-v/-B/g' sambamba_index_wgs_rule
sed -i 's/-v/-B/g' sambamba_index_wts_markdup_rule
sed -i 's/-v/-B/g' sambamba_index_wts_rule
sed -i 's/-v/-B/g' sambamba_merge_facet_rule
sed -i 's/-v/-B/g' sambamba_merge_wgs_rule
sed -i 's/-v/-B/g' sambamba_merge_wts_rule
sed -i 's/-v/-B/g' sambamba_slice_wgs_rule
sed -i 's/-v/-B/g' sambamba_slice_wts_rule
sed -i 's/-v/-B/g' signature_mut_matrix_decomp_rule
sed -i 's/-v/-B/g' signature_mut_matrix_ext_rule
sed -i 's/-v/-B/g' signature_mut_matrix_gen_rule
sed -i 's/-v/-B/g' snpeff_tmb_rule
sed -i 's/-v/-B/g' snpeff_wts_wes_rule
sed -i 's/-v/-B/g' snpsift_annotate_leaves_trio_wgs_rule
sed -i 's/-v/-B/g' snpsift_annotate_TMB_rule
sed -i 's/-v/-B/g' snpsift_annotate_wes_rule
sed -i 's/-v/-B/g' snpsift_dbnsfp_trio_GRCh38_wgs_rule
sed -i 's/-v/-B/g' snpsift_dbnsfp_trio_wgs_rule
sed -i 's/-v/-B/g' snpsift_phastcons_trio_wgs_rule
sed -i 's/-v/-B/g' star_arriba_rule
sed -i 's/-v/-B/g' star_cancer_rule
sed -i 's/-v/-B/g' star_fusion_inspector_rule
sed -i 's/-v/-B/g' star_fusion_rule
sed -i 's/-v/-B/g' tabix_cancer_vcf_wes_rule
sed -i 's/-v/-B/g' tabix_cancer_vcf_wgs_rule
sed -i 's/-v/-B/g' tabix_cancer_vcf_wts_rule
sed -i 's/-v/-B/g' tabix_wes_rule
sed -i 's/-v/-B/g' tmb_curie_rule
sed -i 's/-v/-B/g' vt_decompose_cancer_database_rule
sed -i 's/-v/-B/g' vt_decompose_cancer_tmb_rule
sed -i 's/-v/-B/g' vt_decompose_cancer_wes_rule
sed -i 's/-v/-B/g' vt_decompose_cancer_wgs_rule
sed -i 's/-v/-B/g' vt_decompose_cancer_wts_rule
sed -i 's/-v/-B/g' vt_normalize_uniq_cancer_database_rule
sed -i 's/-v/-B/g' vt_normalize_uniq_cancer_tmb_rule
sed -i 's/-v/-B/g' vt_normalize_uniq_cancer_wes_rule
sed -i 's/-v/-B/g' vt_normalize_uniq_cancer_wgs_rule
sed -i 's/-v/-B/g' vt_normalize_uniq_cancer_wts_rule
sed -i 's/:z / /g' arriba_calling_rule
sed -i 's/:z / /g' arriba_drawing_rule
sed -i 's/:z / /g' arriba_normalization_rule
sed -i 's/:z / /g' bcftools_concat_cancer_rule
sed -i 's/:z / /g' bcftools_concat_wes_rule
sed -i 's/:z / /g' bcftools_concat_wgs_rule
sed -i 's/:z / /g' bcftools_merge_wts_wes_rule
sed -i 's/:z / /g' bcftools_reheader_cancer_rule
sed -i 's/:z / /g' bedtools_wes_rule
sed -i 's/:z / /g' bwa_mem_cancer_rule
sed -i 's/:z / /g' bwa_mem_cancer_rule.bak
sed -i 's/:z / /g' fastq_concat_rule
sed -i 's/:z / /g' fastq_concat_rule
sed -i 's/:z / /g' fusion_catcher_rule
sed -i 's/:z / /g' fusion_normalisation_script_rule
sed -i 's/:z / /g' gatk_applybqsr_wgs_cancer_rule
sed -i 's/:z / /g' gatk_applybqsr_wgs_rule
sed -i 's/:z / /g' gatk_applybqsr_wts_rule
sed -i 's/:z / /g' gatk_base_recalibrator_wgs_cancer_rule
sed -i 's/:z / /g' gatk_base_recalibrator_wgs_rule
sed -i 's/:z / /g' gatk_base_recalibrator_wts_rule
sed -i 's/:z / /g' gatk_filter_mutect_calls_wes_rule
sed -i 's/:z / /g' gatk_filter_mutect_calls_wes_tumor_only_rule
sed -i 's/:z / /g' gatk_haplotype_caller_wgs_cancer_rule
sed -i 's/:z / /g' gatk_haplotype_caller_wts_rule
sed -i 's/:z / /g' gatk_mergevcf_mutect2_rule
sed -i 's/:z / /g' gatk_mutect2_wes_rule
sed -i 's/:z / /g' gatk_mutect2_wes_tumor_only_rule
sed -i 's/:z / /g' gatk_select_variants_mutect_wes_rule
sed -i 's/:z / /g' msisensor_rule
sed -i 's/:z / /g' picard_collect_insert_size_metrics_wes_rule
sed -i 's/:z / /g' picard_collect_insert_size_metrics_wgs_cancer_rule
sed -i 's/:z / /g' picard_collect_insert_size_metrics_wts_rule
sed -i 's/:z / /g' picard_markdup_wgs_rule
sed -i 's/:z / /g' picard_markdup_wts_rule
sed -i 's/:z / /g' quality_wes_changement_rule
sed -i 's/:z / /g' rsync_vcf_tmb_rule
sed -i 's/:z / /g' sambamba_flagstat_wgs_rule
sed -i 's/:z / /g' sambamba_flagstat_wts_rule
sed -i 's/:z / /g' sambamba_index_bam_rule
sed -i 's/:z / /g' sambamba_index_genome_rule
sed -i 's/:z / /g' sambamba_index_wgs_markdup_rule
sed -i 's/:z / /g' sambamba_index_wgs_rule
sed -i 's/:z / /g' sambamba_index_wts_markdup_rule
sed -i 's/:z / /g' sambamba_index_wts_rule
sed -i 's/:z / /g' sambamba_merge_facet_rule
sed -i 's/:z / /g' sambamba_merge_wgs_rule
sed -i 's/:z / /g' sambamba_merge_wts_rule
sed -i 's/:z / /g' sambamba_slice_wgs_rule
sed -i 's/:z / /g' sambamba_slice_wts_rule
sed -i 's/:z / /g' signature_mut_matrix_decomp_rule
sed -i 's/:z / /g' signature_mut_matrix_ext_rule
sed -i 's/:z / /g' signature_mut_matrix_gen_rule
sed -i 's/:z / /g' snpeff_tmb_rule
sed -i 's/:z / /g' snpeff_wts_wes_rule
sed -i 's/:z / /g' snpsift_annotate_leaves_trio_wgs_rule
sed -i 's/:z / /g' snpsift_annotate_TMB_rule
sed -i 's/:z / /g' snpsift_annotate_wes_rule
sed -i 's/:z / /g' snpsift_dbnsfp_trio_GRCh38_wgs_rule
sed -i 's/:z / /g' snpsift_dbnsfp_trio_wgs_rule
sed -i 's/:z / /g' snpsift_phastcons_trio_wgs_rule
sed -i 's/:z / /g' star_arriba_rule
sed -i 's/:z / /g' star_cancer_rule
sed -i 's/:z / /g' star_fusion_inspector_rule
sed -i 's/:z / /g' star_fusion_rule
sed -i 's/:z / /g' tabix_cancer_vcf_wes_rule
sed -i 's/:z / /g' tabix_cancer_vcf_wgs_rule
sed -i 's/:z / /g' tabix_cancer_vcf_wts_rule
sed -i 's/:z / /g' tabix_wes_rule
sed -i 's/:z / /g' tmb_curie_rule
sed -i 's/:z / /g' vt_decompose_cancer_database_rule
sed -i 's/:z / /g' vt_decompose_cancer_tmb_rule
sed -i 's/:z / /g' vt_decompose_cancer_wes_rule
sed -i 's/:z / /g' vt_decompose_cancer_wgs_rule
sed -i 's/:z / /g' vt_decompose_cancer_wts_rule
sed -i 's/:z / /g' vt_normalize_uniq_cancer_database_rule
sed -i 's/:z / /g' vt_normalize_uniq_cancer_tmb_rule
sed -i 's/:z / /g' vt_normalize_uniq_cancer_wes_rule
sed -i 's/:z / /g' vt_normalize_uniq_cancer_wgs_rule
sed -i 's/:z / /g' vt_normalize_uniq_cancer_wts_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' arriba_calling_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' arriba_drawing_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' arriba_normalization_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' bcftools_concat_cancer_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' bcftools_concat_wes_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' bcftools_concat_wgs_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' bcftools_merge_wts_wes_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' bcftools_reheader_cancer_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' bedtools_wes_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' bwa_mem_cancer_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' bwa_mem_cancer_rule.bak
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' fastq_concat_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' fastq_concat_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' fusion_catcher_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' fusion_normalisation_script_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' gatk_applybqsr_wgs_cancer_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' gatk_applybqsr_wgs_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' gatk_applybqsr_wts_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' gatk_base_recalibrator_wgs_cancer_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' gatk_base_recalibrator_wgs_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' gatk_base_recalibrator_wts_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' gatk_filter_mutect_calls_wes_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' gatk_filter_mutect_calls_wes_tumor_only_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' gatk_haplotype_caller_wgs_cancer_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' gatk_haplotype_caller_wts_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' gatk_mergevcf_mutect2_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' gatk_mutect2_wes_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' gatk_mutect2_wes_tumor_only_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' gatk_select_variants_mutect_wes_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' msisensor_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' picard_collect_insert_size_metrics_wes_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' picard_collect_insert_size_metrics_wgs_cancer_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' picard_collect_insert_size_metrics_wts_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' picard_markdup_wgs_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' picard_markdup_wts_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' quality_wes_changement_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' rsync_vcf_tmb_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' sambamba_flagstat_wgs_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' sambamba_flagstat_wts_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' sambamba_index_bam_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' sambamba_index_genome_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' sambamba_index_wgs_markdup_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' sambamba_index_wgs_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' sambamba_index_wts_markdup_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' sambamba_index_wts_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' sambamba_merge_facet_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' sambamba_merge_wgs_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' sambamba_merge_wts_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' sambamba_slice_wgs_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' sambamba_slice_wts_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' signature_mut_matrix_decomp_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' signature_mut_matrix_ext_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' signature_mut_matrix_gen_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' snpeff_tmb_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' snpeff_wts_wes_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' snpsift_annotate_leaves_trio_wgs_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' snpsift_annotate_TMB_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' snpsift_annotate_wes_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' snpsift_dbnsfp_trio_GRCh38_wgs_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' snpsift_dbnsfp_trio_wgs_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' snpsift_phastcons_trio_wgs_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' star_arriba_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' star_cancer_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' star_fusion_inspector_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' star_fusion_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' tabix_cancer_vcf_wes_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' tabix_cancer_vcf_wgs_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' tabix_cancer_vcf_wts_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' tabix_wes_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' tmb_curie_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' vt_decompose_cancer_database_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' vt_decompose_cancer_tmb_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' vt_decompose_cancer_wes_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' vt_decompose_cancer_wgs_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' vt_decompose_cancer_wts_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' vt_normalize_uniq_cancer_database_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' vt_normalize_uniq_cancer_tmb_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' vt_normalize_uniq_cancer_wes_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' vt_normalize_uniq_cancer_wgs_rule
sed -i 's/gitlab-bioinfo.aphp.fr:5000\/sequoia-docker-tools/\/mnt\/beegfs\/software\/seqoia\/singularity/g' vt_normalize_uniq_cancer_wts_rule
