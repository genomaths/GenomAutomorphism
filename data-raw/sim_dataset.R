# ============================================================================ #
#
# ======== Script used to generate the datasets used in the examples ========= #
#
# ============================================================================ #
library(GenomAutomorphism)
library(Biostrings)


aln <- c("ACCTATGTTGGTATT---GCGCTCCAACTCCTTGGCTCTAGCTCACTACAT", 
         "ATCTATGTTGGTATTACGACGCTCCAATTCCTTGGGTCC------CTCCTT")

aln <- DNAStringSet(aln)

usethis::use_data(aln)




