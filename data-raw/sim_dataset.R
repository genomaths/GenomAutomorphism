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

url <- paste0("https://github.com/genomaths/seqalignments/raw/master/CYCS/",
              "mammals_cytochrome_c_(CYCS)_18_sequences.fasta")

cyc_aln <- readDNAMultipleAlignment(filepath = url)

usethis::use_data(cyc_aln, overwrite = TRUE)


nams <- c("human_1", "human_2", "gorilla", "human_3", "human_4",
          "human_5", "human_6", "silvery_gibbon", "white_cheeked_gibbon",
          "françois_langur", "olive_baboon_1", "olive_baboon_2",
          "golden_monkey", "rhesus_monkeys_1", "rhesus_monkeys_2",
          "gelada_baboon_1", "gelada_baboon_2", "orangutan_1", "orangutan_2")

autm <- automorphism(filepath = url, 
                     group = "Z64", 
                     cube = c("ACGT", "TGCA"),
                     cube_alt = c("CATG", "GTAC"),
                     nms = nams)


usethis::use_data(autm, overwrite = TRUE)





