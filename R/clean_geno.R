#' Clean gds genotype function, called by geno_sjlife_gds.R
#'
#' @param n only keep SNPs with at least n non-missing entries
#' @param geno_data genotype matrix from the getGenotype() function
#'
#' @return An n*m matrix where m depends on the number of SNPs that passed QC.

clean_geno <- function(n, geno_data, map_file) {

    # loop through SNPs, remove missing and add two halves
    Gmat <- c()
    removed <- c()
    flipped <- c()
    for (snp_it in 1:nrow(geno_data)) {
        # remove missing
        temp_snp <- geno_data[snp_it, ]
        temp_snp <- temp_snp[which(!is.na(temp_snp))]
        if (length(temp_snp) < n) {
            removed <- c(removed, snp_it)
            next
        }

        # flip to minor allele
        if (mean(temp_snp) > 1) {
            temp_snp <- 2 - temp_snp
            flipped <- c(flipped, snp_it)
        }

        # add SNP to data matrix
        Gmat <- cbind(Gmat, temp_snp[1:n])
    }

    # flip before remove in map file
    temp_new_A <- map_file$A
    temp_new_A[flipped] <- map_file$B[flipped]
    map_file$B[flipped] <- map_file$A[flipped]
    map_file$A <- temp_new_A
    # remove
    map_file <- map_file[-removed, ]

    return(list(cleanG=Gmat, cleanMap = map_file))
}

