require(tidyr)
require(dplyr)

tumor_sample_vcf_VAFs <- read.delim(
    base::as.character(x = snakemake.input[["readcount"]]),
    header = FALSE,
    fill = TRUE,
    col.names = c(
        "CHROM", "POS", "REF", "TUMOR_DP", "c5",
        "count_A", "count_C", "count_G", "count_T", "count_N"
    ),
    sep = "\t",
    stringsAsFactors = FALSE
)

tumor_sample_vcf_VAFs <- dplyr::select(
    tumor_sample_vcf_VAFs, CHROM, POS, REF,
    TUMOR_DP, count_A, count_C, count_G, count_T, count_N
)

tumor_sample_vcf_VAFs <- tumor_sample_vcf_VAFs %>% separate(
    count_A,
    c(
        "A", "A_AD", "A_avg_mapping_quality", "A_avg_basequality", 
        "A_avg_se_mapping_quality", "A_num_plus_strand", "A_num_minus_strand"
    ),
    extra = "drop",
    sep = ":"
)

tumor_sample_vcf_VAFs <- tumor_sample_vcf_VAFs %>% separate(
    count_C, 
    c(
        "C", "C_AD", "C_avg_mapping_quality", "C_avg_basequality", 
        "C_avg_se_mapping_quality", "C_num_plus_strand", "C_num_minus_strand"
    ),
    extra = "drop",
    sep = ":"
)

tumor_sample_vcf_VAFs <- tumor_sample_vcf_VAFs %>% separate(
    count_G, 
    c(
        "G", "G_AD", "G_avg_mapping_quality", "G_avg_basequality", 
        "G_avg_se_mapping_quality", "G_num_plus_strand", "G_num_minus_strand"
    ),
    extra = "drop",
    sep = ":"
)

tumor_sample_vcf_VAFs <- tumor_sample_vcf_VAFs %>% separate(
    count_T, 
    c(
        "T", "T_AD", "T_avg_mapping_quality", "T_avg_basequality", 
        "T_avg_se_mapping_quality", "T_num_plus_strand", "T_num_minus_strand"
    ),
    extra = "drop",
    sep = ":"
)

tumor_sample_vcf_VAFs <- dplyr::select(
    tumor_sample_vcf_VAFs, CHROM, POS, REF, TUMOR_DP, A_AD, C_AD, G_AD, T_AD
)
tumor_sample_vcf_VAFs$TUMOR_DP <- as.numeric(tumor_sample_vcf_VAFs$TUMOR_DP)
tumor_sample_vcf_VAFs$A_AD <- as.numeric(tumor_sample_vcf_VAFs$A_AD)
tumor_sample_vcf_VAFs$C_AD <- as.numeric(tumor_sample_vcf_VAFs$C_AD)
tumor_sample_vcf_VAFs$G_AD <- as.numeric(tumor_sample_vcf_VAFs$G_AD)
tumor_sample_vcf_VAFs$T_AD <- as.numeric(tumor_sample_vcf_VAFs$T_AD)

tumor_sample_vcf_VAFs$VAF_A <- tumor_sample_vcf_VAFs$A_AD/tumor_sample_vcf_VAFs$TUMOR_DP
tumor_sample_vcf_VAFs$VAF_C <- tumor_sample_vcf_VAFs$C_AD/tumor_sample_vcf_VAFs$TUMOR_DP
tumor_sample_vcf_VAFs$VAF_G <- tumor_sample_vcf_VAFs$G_AD/tumor_sample_vcf_VAFs$TUMOR_DP
tumor_sample_vcf_VAFs$VAF_T <- tumor_sample_vcf_VAFs$T_AD/tumor_sample_vcf_VAFs$TUMOR_DP

tumor_sample_vcf_VAFs_TUMOR_DP20 <- tumor_sample_vcf_VAFs[
    which(tumor_sample_vcf_VAFs$TUMOR_DP >= 20),
]
write.csv2(
    tumor_sample_vcf_VAFs_TUMOR_DP20,
    file = snakemake@output[["tumor_dp20"]],
    row.names = FALSE
)
## Looking for VAF columns
vaf.cols <- grep(
    pattern = '^VAF_',
    x = colnames(tumor_sample_vcf_VAFs_TUMOR_DP20)
)
names(vaf.cols) <- sub(
    pattern = "^VAF_",
    replacement = "",
    colnames(tumor_sample_vcf_VAFs_TUMOR_DP20)[vaf.cols]
)

## Looking for AD columns
ad.cols <- grep(
    pattern = "_AD$",
    x = colnames(tumor_sample_vcf_VAFs_TUMOR_DP20)
)
names(ad.cols) <- sub(
    pattern = "_AD$",
    replacement = "",
    colnames(tumor_sample_vcf_VAFs_TUMOR_DP20)[ad.cols]
)

## Preparing new columns
newdata <- as.data.frame(
    matrix(
        NA,
        ncol = 3,
        nrow = nrow(tumor_sample_vcf_VAFs_TUMOR_DP20),
        dimnames = list(NULL, c("ALT", "ALT_AD", "ALT_VAF"))
    )
)

for (k in seq_len(nrow(tumor_sample_vcf_VAFs_TUMOR_DP20))) {
    tmp.vaf.cols <- vaf.cols[
        !toupper(tumor_sample_vcf_VAFs_TUMOR_DP20$REF[k]) == names(vaf.cols)
    ]
    alt.sum <- sum(tumor_sample_vcf_VAFs_TUMOR_DP20[k, tmp.vaf.cols])
    ## If ALT exists
    if(alt.sum > 0) {
        alt.vaf.max <- max(tumor_sample_vcf_VAFs_TUMOR_DP20[k, tmp.vaf.cols])
        alt.bases <- names(tmp.vaf.cols)[
            tumor_sample_vcf_VAFs_TUMOR_DP20[k, tmp.vaf.cols] == alt.vaf.max
        ]
        newdata$ALT[k] <- paste(alt.bases, collapse = "/")
        newdata$ALT_VAF[k] <- tumor_sample_vcf_VAFs_TUMOR_DP20[
            k, vaf.cols[alt.bases[1]]
        ]
        newdata$ALT_AD[k] <- tumor_sample_vcf_VAFs_TUMOR_DP20[
            k, ad.cols[alt.bases[1]]
        ]
    }
}

tumor_sample_vcf_VAFs_TUMOR_DP20 <- cbind(
    tumor_sample_vcf_VAFs_TUMOR_DP20, newdata
)

tumor_sample_vcf_VAFs_TUMOR_DP20_called <- tumor_sample_vcf_VAFs_TUMOR_DP20[
    which(!is.na(tumor_sample_vcf_VAFs_TUMOR_DP20$ALT)),
]
tumor_sample_vcf_VAFs_TUMOR_DP20_called$REF <- toupper(
    tumor_sample_vcf_VAFs_TUMOR_DP20_called$REF
)

tumor_sample_vcf_VAFs_TUMOR_DP20_called_ALT10 <- tumor_sample_vcf_VAFs_TUMOR_DP20_called[
    which(tumor_sample_vcf_VAFs_TUMOR_DP20_called$ALT_AD>=10),
]
tumor_sample_vcf_VAFs_TUMOR_DP20_called_ALT10_VAF3 <- tumor_sample_vcf_VAFs_TUMOR_DP20_called_ALT10[
    which(tumor_sample_vcf_VAFs_TUMOR_DP20_called_ALT10$ALT_VAF>=0.03),
]
tumor_sample_vcf_VAFs_TUMOR_DP20_called_ALT10_VAF3$allele <- apply(
    tumor_sample_vcf_VAFs_TUMOR_DP20_called_ALT10_VAF3[, c("REF", "ALT")],
    1,
    paste,
    collapse = "/")
tumor_sample_vcf_VAFs_TUMOR_DP20_called_ALT10_VAF3$strand <- rep(
    "." , length(tumor_sample_vcf_VAFs_TUMOR_DP20_called_ALT10_VAF3$CHROM)
)
tumor_sample_vcf_VAFs_TUMOR_DP20_called_ALT10_VAF3$refcount <- tumor_sample_vcf_VAFs_TUMOR_DP20_called_ALT10_VAF3$TUMOR_DP - tumor_sample_vcf_VAFs_TUMOR_DP20_called_ALT10_VAF3$ALT_AD
tumor_sample_vcf_VAFs_TUMOR_DP20_called_ALT10_VAF3$identifier <- apply(
    tumor_sample_vcf_VAFs_TUMOR_DP20_called_ALT10_VAF3[, c('refcount','ALT_AD','TUMOR_DP','ALT_VAF')],
    1,
    paste,
    collapse = "|"
)

vep_input <- tumor_sample_vcf_VAFs_TUMOR_DP20_called_ALT10_VAF3[c(1,2,2,16,17,19)]
write.table(
    vep_input,
    file = base::as.character(x = snakemake@output[["filtered"]]),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    sep = "\t"
)

