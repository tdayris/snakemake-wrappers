library(ensemblVEP)
library(dplyr)
library(VariantAnnotation)

hotspot_list <- read.table(
    file = base::as.character(x = snakemake@input[["cancer_genes"]]),
    header = TRUE,
    stringsAsFactors = FALSE
)

files_to_process <- sapply(
    snakemake@input[["vcfs"]],
    function(path) as.character(x = path)
)

for (invcf in files_to_process) {
    input_hc <- readVcf(invcf, "hg38")
    if(length(as.character(seqnames(input_hc)))>0){
        name_normal <- NA
        name_tumor <- rownames(colData(input_hc))
        CHROM <- as.character(x = seqnames(input_hc))
        POS <- as.character(x = start(input_hc))
        REF <- as.character(x = granges(input_hc)$REF)
        ALT <- NULL    

        for (i in 1:length(rowRanges(input_hc)$ALT)) {
            for (j in 1:length(rowRanges(input_hc)$ALT[[i]])) {
                if (is.null(ALT[i]) || is.na(ALT[i])) {
                    ALT[i] <- as.character(rowRanges(input_hc)$ALT[[i]][[j]])
                } else {
                    if (as.character(rowRanges(input_hc)$ALT[[i]][[j]]) == "") {
                        ALT[i] <- paste(ALT[i], "*", sep = ",")
                    } else {
                        ALT[i] <- paste(
                            ALT[i],
                            as.character(rowRanges(input_hc)$ALT[[i]][[j]]),
                            sep = ","
                        )
                    }
                }
            }
        }    

        csq <- parseCSQToGRanges(input_hc, VCFRowID = rownames(x = input_hc))
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
        tmp <- list()    

        for (i in strsplit(x = csq_df_aggregated$SYMBOL, split = ",")) {
            tmp <- c(tmp, list(i %in% hotspot_list$Hotspot_SYMBOL))
        }    

        Hotspot <- sapply(tmp, paste, collapse = ",")
        Patient <- strsplit(x = name_tumor, split = "_")[[1]][1]
        Patient <- rep(Patient, length(x = CHROM))
        SampleName <- name_tumor
        SampleName <- rep(SampleName, length(x = CHROM))
        Method <- rep("haplotypecaller", length(x = CHROM))
        normal_ref_count <- rep(NA, length(x = CHROM))
        normal_alt_count <- rep(NA, length(x = CHROM))
        normal_depth <- rep(NA, length(x = CHROM))
        normal_alt_freq <- rep(NA, length(x = CHROM))
        tumor_ref_count <- NULL    

        for (i in 1:length(geno(input_hc)$AD[, name_tumor])) {
            tumor_is_null <- is.null(tumor_ref_count[i])
            tumor_is_na <- is.na(tumor_ref_count[i])
            if (tumor_is_null || tumor_is_na) {
                tumor_ref_count[i] <- as.character(
                    x = geno(input_hc)$AD[, name_tumor][[i]][[1]]
                )
            } else {
                predicate <- as.character(
                    x = geno(input_hc)$AD[, name_tumor][[i]][[1]]
                )
                if (predicate == "") {
                    tumor_ref_count[i] <- paste(
                        tumor_ref_count[i], "*", sep = ","
                    )
                } else {
                    tumor_ref_count[i] <- paste(
                        tumor_ref_count[i],
                        as.character(
                            x = geno(input_hc)$AD[, name_tumor][[i]][[1]]
                        ),
                        sep = ","
                    )
                }
            }
        }    

        tumor_alt_count <- NULL
        for (i in 1:length(geno(input_hc)$AD[, name_tumor])) {
            if (is.null(tumor_alt_count[i]) || is.na(tumor_alt_count[i])) {
                tumor_alt_count[i] <- as.character(
                    x = geno(input_hc)$AD[, name_tumor][[i]][[2]]
                )
            } else {
                predicate <- as.character(
                    x = geno(input_hc)$AD[, name_tumor][[i]][[2]]
                )
                if (predicate == "") {
                    tumor_alt_count[i] <- paste(
                        tumor_alt_count[i], "*", sep = ","
                    )
                } else {
                    tumor_alt_count[i] <- paste(
                        tumor_alt_count[i],
                        as.character(
                            x = geno(input_hc)$AD[, name_tumor][[i]][[2]]
                        ),
                        sep = ","
                    )
                }
            }
        }    

        tumor_depth <- rowSums(
            cbind(
                as.numeric(x = tumor_ref_count), 
                as.numeric(x = tumor_alt_count)
            )
        )    

        tumor_alt_freq <- as.numeric(x = tumor_alt_count) / tumor_depth
        varID <- unique(names(csq))
        STRAND <- csq_df_aggregated$strand
        FILTER <- as.character(granges(input_hc)$FILTER)
        Feature_RefSeq <- csq_df_aggregated$Feature
        CANONICAL <- csq_df_aggregated$CANONICAL
        pos_in_canonical_vectors <- sapply(
            strsplit(x = csq_df_aggregated$CANONICAL, split = ","),
            function(y) which(y == "YES")
        )
        Feature_RefSeq_splitted <- strsplit(
            x = Feature_RefSeq, split = ","
        )
        Canonical_NM <- NA
        Canonical_Protein_position <- NA
        Canonical_Amino_acids <- NA
        Canonical_IMPACT <- NA 
        Canonical_Consequence <- NA
        Amino_acids_splitted <- NA
        Consequence_splitted <- NA
        IMPACT_splitted <- NA
        Protein_position_splitted <- NA   
        if(!identical(pos_in_canonical_vectors[[1]], integer(0)) && !is.integer(pos_in_canonical_vectors)){
          for (i in 1:length(pos_in_canonical_vectors)){
            if (length(pos_in_canonical_vectors[[i]]) > 0) {
              for (j in 1:length(pos_in_canonical_vectors[[i]])) {
                if (is.null(Canonical_NM[i]) || is.na(Canonical_NM[i])) {
                  Canonical_NM[i] <- Feature_RefSeq_splitted[[i]][pos_in_canonical_vectors[[i]][j]]
                } else {
                  Canonical_NM[i] <- paste(
                    Canonical_NM[i],
                    Feature_RefSeq_splitted[[i]][pos_in_canonical_vectors[[i]][j]],
                    sep = ","
                  )
                }
              }
            } else {
              Canonical_NM[i] <- NA
            }
          }    
          
          Protein_position_splitted <- strsplit(
            x = Protein_position, split = ","
          )
          Canonical_Protein_position <- NULL    
          
          for (i in 1:length(pos_in_canonical_vectors)) {
            if (length(pos_in_canonical_vectors[[i]]) > 0) {
              for (j in 1:length(pos_in_canonical_vectors[[i]])){
                if (is.null(Canonical_Protein_position[i]) || is.na(Canonical_Protein_position[i])) {
                  Canonical_Protein_position[i] <- Protein_position_splitted[[i]][pos_in_canonical_vectors[[i]][j]]
                } else {
                  Canonical_Protein_position[i] <- paste(
                    Canonical_Protein_position[i],
                    Protein_position_splitted[[i]][pos_in_canonical_vectors[[i]][j]],
                    sep = ","
                  )
                }
              }
            } else {
              Canonical_Protein_position[i] <- NA
            }
          }    
          
          Amino_acids_splitted <- strsplit(Amino_acids, split = ",")
          Canonical_Amino_acids <- NULL    
          
          for (i in 1:length(pos_in_canonical_vectors)){
            if (length(pos_in_canonical_vectors[[i]]) > 0) {
              for (j in 1:length(pos_in_canonical_vectors[[i]])){
                if (is.null(Canonical_Amino_acids[i]) || is.na(Canonical_Amino_acids[i])) {
                  Canonical_Amino_acids[i] <- Amino_acids_splitted[[i]][pos_in_canonical_vectors[[i]][j]]
                } else {
                  Canonical_Amino_acids[i] <- paste(
                    Canonical_Amino_acids[i],
                    Amino_acids_splitted[[i]][pos_in_canonical_vectors[[i]][j]],
                    sep = ","
                  )
                }
              }
            } else {
              Canonical_Amino_acids[i] <- NA
            }
          }
          
          IMPACT_splitted <- strsplit(IMPACT,",")
          Canonical_IMPACT <- NULL
          
          for (i in 1:length(pos_in_canonical_vectors)) {
            if (length(pos_in_canonical_vectors[[i]]) > 0) {
              for (j in 1:length(pos_in_canonical_vectors[[i]])) {
                if (is.null(Canonical_IMPACT[i]) || is.na(Canonical_IMPACT[i])) {
                  Canonical_IMPACT[i] <- IMPACT_splitted[[i]][pos_in_canonical_vectors[[i]][j]]
                } else {
                  Canonical_IMPACT[i] <- paste(
                    Canonical_IMPACT[i],
                    IMPACT_splitted[[i]][pos_in_canonical_vectors[[i]][j]],
                    sep = ","
                  )
                }
                
              }
            } else {
              Canonical_IMPACT[i] <- NA
            }
          }    
          
          Consequence_splitted <- strsplit(x = Consequence, split = ",")
          Canonical_Consequence <- NULL    
          
          for (i in 1:length(pos_in_canonical_vectors)) {
            if (length(pos_in_canonical_vectors[[i]]) > 0) {
              for (j in 1:length(pos_in_canonical_vectors[[i]])){
                if (is.null(Canonical_Consequence[i]) || is.na(Canonical_Consequence[i])) {
                  Canonical_Consequence[i] <- Consequence_splitted[[i]][pos_in_canonical_vectors[[i]][j]]
                } else {
                  Canonical_Consequence[i] <- paste(
                    Canonical_Consequence[i],
                    Consequence_splitted[[i]][pos_in_canonical_vectors[[i]][j]],
                    sep = ","
                  )
                }
              }
            } else {
              Canonical_Consequence[i] <- NA
            }
          }
        }

        CODING <- grepl("HIGH|MODERATE", IMPACT)
        EXON_splitted <- strsplit(x = csq_df_aggregated$EXON, split = ",")
        test <- sapply(EXON_splitted, function(y) which(y != "NA"))
        inExon <- NA
        if(!identical(test[[1]], integer(0)) && !is.integer(test)){
          inExon <- lengths(test[1,])
        }
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
            Canonical_Amino_acids, Canonical_IMPACT, Canonical_Consequence, CODING,
            inExon, BIOTYPE, Position_Type, SIFT, PolyPhen, dbSNP_COSMIC, uniqueID
        )    

        write.table(
            table_to_write,
            file = snakemake@output[["tsv"]],
            quote = FALSE,
            sep = "\t",
            row.names = FALSE
        )
        rm(
            ALT, Amino_acids, Amino_acids_splitted, BIOTYPE, CANONICAL,
            Canonical_Amino_acids, Canonical_Consequence, Canonical_IMPACT,
            Canonical_NM, Canonical_Protein_position, CHROM, CODING, Consequence,
            Consequence_splitted, csq, csq_df, csq_df_aggregated, dbSNP_COSMIC,
            EXON_splitted, Feature_RefSeq, Feature_RefSeq_splitted, FILTER,
            Gene_Name, Hotspot, i, IMPACT, IMPACT_splitted, inExon, j, Method,
            name_normal, name_tumor, normal_alt_count, normal_alt_freq,
            normal_depth, normal_ref_count, input_hc, Patient, PolyPhen, POS,
            pos_in_canonical_vectors, Position_Type, Protein_position,
            Protein_position_splitted, REF, SampleName, SIFT, STRAND,
            table_to_write, test, tmp, tumor_alt_count, tumor_alt_freq,
            tumor_depth, tumor_ref_count, uniqueID, VARIANT_CLASS, varID
        )
    }else{
        table_to_write <- data.frame(
            "CHROM", "POS", "REF", "ALT", "Method", "Gene_Name", "VARIANT_CLASS", 
            "Consequence", "IMPACT", "Protein_position", "Amino_acids", "Hotspot", 
            "Patient", "SampleName", "normal_ref_count", "tumor_ref_count", 
            "normal_alt_count", "tumor_alt_count", "normal_depth", "tumor_depth", 
            "normal_alt_freq", "tumor_alt_freq", "varID", "STRAND", "FILTER", 
            "Feature_RefSeq", "CANONICAL", "Canonical_NM", "Canonical_Protein_position",
            "Canonical_Amino_acids", "Canonical_IMPACT", "Canonical_Consequence", "CODING",
            "inExon", "BIOTYPE", "Position_Type", "SIFT", "PolyPhen", "dbSNP_COSMIC", "uniqueID"
        )    
        write.table(
            table_to_write,
            file = snakemake@output[["tsv"]],
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE
        )
    }
}
