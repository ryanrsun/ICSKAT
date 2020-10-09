# Analyze SJLIFE data
args <- commandArgs(trailingOnly=TRUE)
aID <- as.numeric(args[1])
chr <- as.numeric(args[2])
group <- as.numeric(args[3])

library(dplyr)
library(magrittr)
library(GWASTools)
library(SNPRelate)
library(data.table)
source("/rsrch3/home/biostatistics/rsun3/github/ICSKAT/R/ICSKAT.R")
source("/rsrch3/home/biostatistics/rsun3/github/ICSKAT/R/clean_geno.R")
source("/rsrch3/home/biostatistics/rsun3/github/ICSKAT/R/ICSKAT_fit_null.R")
source("/rsrch3/home/biostatistics/rsun3/github/ICSKAT/R/make_IC_dmat.R")

output_dir <- "/rsrch3/home/biostatistics/rsun3/sjlife/output/looseTeeth"
outRoot <- "looseTeethRes"
max_gene_size <- 5000

# load the gene information
setwd("/rsrch3/home/biostatistics/rsun3/github/LungCancerAssoc/Data")
load(file="ensembl_refgene_hg19_20180109.rda")
all_genes <- ensembl_refgene_hg19_20180109 %>%
	filter(Chr == chr) %>%
	filter(Notes == 0) %>%
	mutate(Diff = txEnd - txStart) %>%
	arrange(Diff)

# load outcome data
setwd("/rsrch3/home/biostatistics/rsun3/sjlife/data/outcomes")
outcomeTimes <- fread("lost6th.csv", data.table=FALSE) %>%
	filter(!is.na(left) & !is.na(right) & !is.na(sjlid) & !is.na(group_ctrl0_surv1)) %>%
	filter(group_ctrl0_surv1 == group)

# only use full covariates
# start with White only since we don't have PCs?
covarTab <- fread("cardiomyopathyCovar.txt", data.table=FALSE) %>%
	filter(as.character(sjlid) %in% as.character(outcomeTimes$sjlid)) %>%
	filter(as.character(race) %in% c("White")) %>%
	mutate(gender = ifelse(gender == "Male", 0, 1)) %>%	
	select(sjlid, agediag, gender, race)

# merge covariates, times
outcomeDat <- merge(outcomeTimes, covarTab, by="sjlid")

# Pick a section of genes to analyze based on aID argument
genes_per_chunk <- 10
start_row <- genes_per_chunk * (aID - 1) + 1
end_row <- min(genes_per_chunk * aID, nrow(all_genes))
buffer <- 5000

#------------------------------------------------------------------#
# End possibly varying parameters
#------------------------------------------------------------------#

# Cut to that chromosome
gene_info <- all_genes %>% filter(Chr == chr) %>% 
	arrange(Diff) %>%
	slice(start_row:min(nrow(.), end_row))

# if the .gds does not exist, will create a new one
# warning, this will be the same size as the WGS data, could be very big
setwd("/rsrch3/home/biostatistics/rsun3/sjlife/data/geno")
fname <- paste0("SJLIFE.GERMLINE.2364.GATKv3.4.vqsr.release.0714_chr", chr, ".PASS.decomposed.sjlid_eur")
gdsfile <- paste0(fname, ".gds")
if (!file.exists(gdsfile)) {
    bed_file <- paste0(fname, ".bed")
    fam_file <- paste0(fname, ".fam")
    bim_file <- paste0(fname, ".bim")
    snpgdsBED2GDS(bed_file, fam_file, bim_file, gdsfile, family = TRUE,
                  cvt.chr = "int", cvt.snpid = "int")
}
# read in the gds file
gds <- GdsGenotypeReader(gdsfile, YchromCode=24L, XYchromCode=25L)
# matches indices with positions
pos <- getPosition(gds)
# the scan IDs for those in the outcome dataset
all_IDs <- getScanID(gds)
scan_idx <- which(all_IDs %in% as.character(outcomeDat$sjlid))
scan_ids <- all_IDs[scan_idx]
# other parameters
n <- length(scan_ids)
impute <- FALSE
rareMAF <- 0.05

# results data frame
resultsDF <- data.frame(gene = gene_info$HGNC_name, start=gene_info$txStart,
                        end=gene_info$txEnd, q=NA, cleanq=NA, rareq=NA, skatp=NA,
                        burdenp=NA, complex=NA, Wskatp=NA, Wburdenp=NA, Wcomplex=NA,
                        Rskatp=NA, Rburdenp=NA, Rcomplex=NA, RWskatp=NA, RWburdenp=NA,
                        RWcomplex=NA)

# loop through genes
for (gene_it in 1:nrow(gene_info)) {

    # pick out the indices
    start_idx <- min(which(pos >= gene_info$txStart[gene_it] - buffer))
    end_idx <- max(which(pos <= gene_info$txEnd[gene_it] + buffer))

    # stop if only one SNP
    q <- end_idx - start_idx + 1
    resultsDF$q[gene_it] <- q
    if (q <= 1) (next)

    # extract genotypes
    geno_data <- getGenotypeSelection(gds, snp=(start_idx:end_idx), scan=scan_idx)
    alleleA <- getAlleleA(gds, index=(start_idx:(start_idx+q-1)))
    alleleB <- getAlleleB(gds, index=(start_idx:(start_idx+q-1)))

    # make a map file
    map_file <- data.frame(BP=pos[start_idx:end_idx], A=alleleA, B=alleleB)

    # remove non-SNPs
    simpleSNPs <- which(!(nchar(alleleA) > 1 | alleleA == "-" | nchar(alleleB) > 1 | alleleB == "-"))
    geno_data <- geno_data[simpleSNPs, ]
    map_file <- map_file[simpleSNPs, ]

    # clean the genotypes
    cleaned_output <- clean_geno(n=n, geno_data=geno_data, map_file=map_file)
    cleanG = cleaned_output$cleanG
    cleanMap <- cleaned_output$cleanMap
    if (nrow(cleanMap) < 2){
        resultsDF$cleanq[gene_it] <- nrow(cleanMap)
        next
    }

    # some will have MAF 0
    MAFs <- apply(cleanG, 2, mean) / 2
    MAFpos <- which(MAFs > 0)
    cleanG <- cleanG[, MAFpos]
    cleanMap <- cleanMap[MAFpos, ]
    MAFs <- MAFs[MAFpos]
    resultsDF$cleanq[gene_it] <- length(MAFs)
    if (length(MAFs) < 2 | length(MAFs) > max_gene_size) {next}
    # weights by beta
    #weights <- dbeta(MAFs, 1, 25)
    #GW <- cleanG %*% diag(weights)

    # order the genotype and outcome matrices
		cleanG <- cleanG[which(rownames(cleanG) %in% outcomeDat$sjlid), ]		
		#GW <- GW[which(rownames(GW) %in% teethData$sjlid), ]
		outcomeDat <- outcomeDat[match(rownames(cleanG), as.character(outcomeDat$sjlid)), ]

    # make design matrices
    dmats <- make_IC_dmat(X=as.matrix(outcomeDat %>% select(agediag, gender)), 
				lt=outcomeDat$left, rt=outcomeDat$right)

    # fit null model
    obs_ind <- as.numeric(outcomeDat$right < 999)
    tpos_ind <- as.numeric(outcomeDat$left > 0)
    null_fit <- ICSKAT_fit_null(init_beta=c(rep(0, 2), 1, 1, 0.5), lt=outcomeDat$left, rt=outcomeDat$right,
                                left_dmat=dmats$left_dmat,
                                right_dmat=dmats$right_dmat,
                                obs_ind=obs_ind, tpos_ind=tpos_ind)
    # get pvalue
    skat_output <- ICskat(left_dmat=dmats$left_dmat, tpos_ind=tpos_ind, obs_ind=obs_ind,
                          right_dmat=dmats$right_dmat, G=cleanG, lt=outcomeDat$left, rt=outcomeDat$right,
                          null_beta=as.numeric(null_fit$beta_fit), Itt=null_fit$Itt)
    # get weighted p-value
    #weighted_skat_output <- ICskat(left_dmat=dmats$left_dmat, tpos_ind=tpos_ind, obs_ind=obs_ind,
    #                               right_dmat=dmats$right_dmat, G=GW, lt=cmData$U, rt=cmData$V,
    #                               null_beta=as.numeric(null_fit$beta_fit), Itt=null_fit$Itt)
    # regular skat results
    resultsDF$skatp[gene_it] <- skat_output$p_SKAT
    resultsDF$burdenp[gene_it] <- skat_output$p_burden
    resultsDF$complex[gene_it] <- skat_output$complex
    # weighted results
    #resultsDF$Wskatp[gene_it] <- weighted_skat_output$p_SKAT
    #resultsDF$Wburdenp[gene_it] <- weighted_skat_output$p_burden
    #resultsDF$Wcomplex[gene_it] <- weighted_skat_output$complex

    cat(gene_it)
}

# close the gds connection
close(gds)

# write results
setwd(output_dir)
write.table(resultsDF, paste0(outRoot, aID, ".txt"), append=F, quote=F, row.names=F, col.names=T, sep="\t")


