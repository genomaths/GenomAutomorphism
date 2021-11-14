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


URL <- paste0("https://github.com/genomaths/seqalignments/raw/master/BRCA1/",
              "brca1_primates_dna_repair_20_sequences.fasta")


nams <- c("human_1","human_2","gorilla_1","gorilla_2","gorilla_3",
          "chimpanzee_1","chimpanzee_2","chimpanzee_3","chimpanzee_4",
          "bonobos_1","bonobos_2","bonobos_3","bonobos_4","silvery_gibbon_1",
          "silvery_gibbon_1","silvery_gibbon_3","golden_monkey_1",
          "golden_monkey_2","gelada_baboon","bolivian_monkey")

brca1_autm <- automorphism(filepath = URL, 
                           group = "Z64", 
                           cube = c("ACGT", "TGCA"),
                           cube_alt = c("CATG", "GTAC"),
                           nms = nams)
usethis::use_data(brca1_autm, overwrite = TRUE)






URL <- paste0("https://github.com/genomaths/seqalignments/raw/master/", 
              "COVID-19/AY390556.1_265-13398_13398-21485_RNA-POL_SARS_COVI_GZ02.fas")


autms <- automorphism(filepath = URL, 
                      group = "Z64", 
                      cube = c("ACGT", "TGCA"),
                      cube_alt = c("CATG", "GTAC"))
autms

usethis::use_data(autms, overwrite = TRUE)
