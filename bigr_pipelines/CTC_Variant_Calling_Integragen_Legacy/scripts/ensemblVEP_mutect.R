require(ensemblVEP)
require(dplyr)

hotspot_list <- read.table(
    base::as.character(x = snakemake@input[["cancer_genes"]]),
    header = TRUE,
    stringsAsFactors = FALSE
)

files_to_process <- sapply(
    snakemake@input[["vcfs"]],
    base::as.character
)

for (invcf in files_to_process) {

    input_mutect <- readVcf(invcf, "hg38")
    name_normal <- meta(header(input_mutect))$SAMPLE["NORMAL", "SampleName"]
    name_tumor <- meta(header(input_mutect))$SAMPLE["TUMOR", "SampleName"]
    CHROM <- as.character(seqnames(input_mutect))
    POS <- as.character(start(input_mutect))
    REF <- as.character(granges(input_mutect)$REF)
    ALT <- as.character(unlist(granges(input_mutect)$ALT))
    csq <- parseCSQToGRanges(input_mutect, VCFRowID = rownames(input_mutect))
    csq_df <- as.data.frame(csq, row.names = NULL)
    csq_df_aggregated <- aggregate(
        csq_df, list(csq_df$VCFRowID), paste, collapse = ","
    )
    Gene_Name <- csq_df_aggregated$SYMBOL
    VARIANT_CLASS <- csq_df_aggregated$VARIANT_CLASS
    Consequence <- csq_df_aggregated$Consequence
    IMPACT <- csq_df_aggregated$IMPACT
    Protein_position <- csq_df_aggregated$Protein_position
    Amino_acids <- csq_df_aggregated$Amino_acids
    tmp<-list()

    for (i in strsplit(csq_df_aggregated$SYMBOL, split = ",")){
        tmp <- c(tmp, list(i %in% hotspot_list$Hotspot_SYMBOL))
    }
    Hotspot <- sapply(tmp, paste, collapse=",")
    Patient <- strsplit(name_tumor, split = "_")[[1]][1]
    Patient <- rep(Patient, length(CHROM))
    SampleName <- name_tumor
    SampleName <- rep(SampleName, length(CHROM))
    Method <- rep("mutect", length(CHROM))
    normal_ref_count <- matrix(
        unlist(geno(input_mutect)$AD[, "NORMAL"]), ncol = 2, byrow = TRUE
    )[, 1]
    tumor_ref_count <- matrix(
        unlist(geno(input_mutect)$AD[, "TUMOR"]), ncol = 2, byrow = TRUE
    )[, 1]
    normal_alt_count <- matrix(
        unlist(geno(input_mutect)$AD[, "NORMAL"]), ncol = 2, byrow = TRUE
    )[, 2]
    tumor_alt_count <- matrix(
        unlist(geno(input_mutect)$AD[, "TUMOR"]), ncol = 2, byrow = TRUE
    )[, 2]
    normal_depth <- rowSums(
        matrix(
            unlist(geno(input_mutect)$AD[, "NORMAL"]),
            ncol = 2,
            byrow = TRUE
        )
    )
    tumor_depth <- rowSums(
        matrix(
            unlist(geno(input_mutect)$AD[, "TUMOR"]),
            ncol = 2,
            byrow = TRUE
        )
    )
    normal_alt_freq <- normal_alt_count / normal_depth
    tumor_alt_freq <- tumor_alt_count / tumor_depth
    varID <- unique(names(csq))
    STRAND <- csq_df_aggregated$strand
    FILTER <- as.character(granges(input_mutect)$FILTER)
    Feature_RefSeq <- csq_df_aggregated$Feature
    CANONICAL <- csq_df_aggregated$CANONICAL
    pos_in_canonical_vectors <- sapply(
        strsplit(csq_df_aggregated$CANONICAL, split = ","),
        function(y) which(y == "YES")
    )
    Feature_RefSeq_splitted <- strsplit(Feature_RefSeq, split = ",")
    Canonical_NM <- NULL

    for (i in 1:length(pos_in_canonical_vectors)) {
        if (length(pos_in_canonical_vectors[[i]]) > 0) {
            for (j in 1:length(pos_in_canonical_vectors[[i]])) {
                if (is.null(Canonical_NM[i]) || is.na(Canonical_NM[i])) {
                    Canonical_NM[i] <- Feature_RefSeq_splitted[[j]][pos_in_canonical_vectors[[i]]]
                } else {
                    Canonical_NM[i] <- paste(
                        Canonical_NM[i],
                        Feature_RefSeq_splitted[[j]][pos_in_canonical_vectors[[i]]],
                        sep = ","
                    )
                }
            }
        } else {
            Canonical_NM[i] <- NA
        }
    }
    Protein_position_splitted <- strsplit(Protein_position, split = ",")
    Canonical_Protein_position <- NULL

    for (i in 1:length(pos_in_canonical_vectors)) {
        if (length(pos_in_canonical_vectors[[i]]) > 0) {
            for (j in 1:length(pos_in_canonical_vectors[[i]])) {
                if (is.null(Canonical_Protein_position[i]) || is.na(Canonical_Protein_position[i])) {
                    Canonical_Protein_position[i] <- Protein_position_splitted[[j]][pos_in_canonical_vectors[[i]]]
                }else{
                    Canonical_Protein_position[i] <- paste(
                        Canonical_Protein_position[i],
                        Protein_position_splitted[[j]][pos_in_canonical_vectors[[i]]],
                        sep = ","
                    )
                }
            }
        }else{
            Canonical_Protein_position[i] <- NA
        }
    }

    Amino_acids_splitted <- strsplit(Amino_acids, split = ",")
    Canonical_Amino_acids <- NULL

    for (i in 1:length(pos_in_canonical_vectors)) {
        if (length(pos_in_canonical_vectors[[i]]) > 0) {
            for (j in 1:length(pos_in_canonical_vectors[[i]])) {
                if (is.null(Canonical_Amino_acids[i]) || is.na(Canonical_Amino_acids[i])) {
                    Canonical_Amino_acids[i] <- Amino_acids_splitted[[j]][pos_in_canonical_vectors[[i]]]
                } else {
                    Canonical_Amino_acids[i] <- paste(
                        Canonical_Amino_acids[i],
                        Amino_acids_splitted[[j]][pos_in_canonical_vectors[[i]]],
                        sep = ","
                    )
                }
            }
        } else {
            Canonical_Amino_acids[i] <- NA
        }
    }

    IMPACT_splitted <- strsplit(IMPACT, split = ",")
    Canonical_IMPACT <- NULL
    
    for (i in 1:length(pos_in_canonical_vectors)) {
        if (length(pos_in_canonical_vectors[[i]]) > 0) {
            for (j in 1:length(pos_in_canonical_vectors[[i]])) {
                if (is.null(Canonical_IMPACT[i]) || is.na(Canonical_IMPACT[i])) {
                    Canonical_IMPACT[i] <- IMPACT_splitted[[j]][pos_in_canonical_vectors[[i]]]
                }else {
                    Canonical_IMPACT[i] <- paste(
                        Canonical_IMPACT[i],
                        IMPACT_splitted[[j]][pos_in_canonical_vectors[[i]]],
                        sep = ","
                    )
                }
            }
        } else {
            Canonical_IMPACT[i] <- NA
        }
    }

    Consequence_splitted <- strsplit(Consequence, split = ",")
    Canonical_Consequence <- NULL

    for (i in 1:length(pos_in_canonical_vectors)) {
        if (length(pos_in_canonical_vectors[[i]]) > 0) {
            for (j in 1:length(pos_in_canonical_vectors[[i]])) {
                if(is.null(Canonical_Consequence[i]) || is.na(Canonical_Consequence[i])){
                    Canonical_Consequence[i] <- Consequence_splitted[[j]][pos_in_canonical_vectors[[i]]]
                } else {
                    Canonical_Consequence[i] <- paste(
                        Canonical_Consequence[i],
                        Consequence_splitted[[j]][pos_in_canonical_vectors[[i]]],
                        sep = ","
                    )
                }
            }
        } else {
            Canonical_Consequence[i] <- NA
        }
    }

    CODING <- ifelse(IMPACT == "HIGH" | IMPACT == "MODERATE", "TRUE", "FALSE")
    EXON_splitted <- strsplit(csq_df_aggregated$EXON, split = ",")
    test <- sapply(EXON_splitted, function(y) which(y != "NA"))
    inExon <- lengths(test)
    BIOTYPE <- csq_df_aggregated$BIOTYPE
    Position_Type <- csq_df_aggregated$EXON
    SIFT <- csq_df_aggregated$SIFT
    PolyPhen <- csq_df_aggregated$PolyPhen
    dbSNP_COSMIC <- csq_df_aggregated$Existing_variation
    uniqueID <- paste(SampleName, varID, sep = "_")
    table_to_write <- data.frame(
        CHROM, POS, REF, ALT, Method, Gene_Name, VARIANT_CLASS,
        Consequence, IMPACT, Protein_position, Amino_acids, Hotspot,
        Patient, SampleName, normal_ref_count, tumor_ref_count,
        normal_alt_count, tumor_alt_count, normal_depth, tumor_depth,
        normal_alt_freq, tumor_alt_freq, varID, STRAND, FILTER,
        Feature_RefSeq, CANONICAL, Canonical_NM, Canonical_Protein_position,
        Canonical_Amino_acids, Canonical_IMPACT, Canonical_Consequence,
        CODING, inExon, BIOTYPE, Position_Type, SIFT, PolyPhen,
        dbSNP_COSMIC, uniqueID
    )
    write.table(
        table_to_write,
        file = snakemake@output[["tsv"]],
        quote = FALSE,
        sep = "\t",
        row.names = FALSE
    )
    rm(
        ALT,Amino_acids, Amino_acids_splitted, BIOTYPE, CANONICAL,
        Canonical_Amino_acids, Canonical_Consequence, Canonical_IMPACT,
        Canonical_NM, Canonical_Protein_position, CHROM, CODING, Consequence,
        Consequence_splitted, csq, csq_df, csq_df_aggregated, dbSNP_COSMIC,
        EXON_splitted, Feature_RefSeq, Feature_RefSeq_splitted, FILTER,
        Hotspot,i,IMPACT,IMPACT_splitted,inExon,j,Method, name_normal,
        name_tumor,normal_alt_count,normal_alt_freq,normal_depth,
        normal_ref_count, input_mutect,Patient,PolyPhen,POS, Gene_Name,
        pos_in_canonical_vectors,Position_Type, Protein_position,
        Protein_position_splitted,REF,SampleName,SIFT,STRAND,table_to_write,
        test,tmp,tumor_alt_count,tumor_alt_freq,tumor_depth,tumor_ref_count,
        uniqueID, VARIANT_CLASS,varID
    )
}