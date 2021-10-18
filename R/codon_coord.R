## Copyright (C) 2021 Robersy Sanchez <https://genomaths.com/>
## Author: Robersy Sanchez This file is part of the R package
## 'GenomAutomorphism'.  'GenomAutomorphism' is free
## software: you can redistribute it and/or modify it under the
## terms of the GNU General Public License as published by the Free
## Software Foundation, either version 3 of the License, or (at
## your option) any later version.  This program is distributed in
## the hope that it will be useful, but WITHOUT ANY WARRANTY;
## without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for
## more details.  You should have received a copy of the GNU
## General Public License along with this program; if not, see
## <http://www.gnu.org/licenses/>.

#' @rdname codon_coord
#' @aliases base_coord
#' @title Represent codon (DNA base triplet) with its coordinates in a given 
#' Abelian group
#' @description Given a string denoting a codon or base from the DNA (or RNA) 
#' alphabet and a genetic-code Abelian group as given in reference (1).
#' @param codon A character string of DNA/RNA base-triplets (i.e., with
#' letter from the DNA/RNA alphabet: "A","C","G","T", and "U").
#' @param cube A character string denoting one of the 24 Genetic-code cubes,
#' as given in references (2 2 3).
#' @param group A character string denoting the group representation for the 
#' given base or codon as shown in reference (1).
#' @details Symbols "-" and "N" usually found in DNA sequence alignments to 
#' denote gaps and missing/unknown bases are represented by the number: '-1' on
#' Z4 and '0' in Z5. In Z64 the symbol 'NA' will be returned for codons 
#' including symbols "-" and "N".
#' @seealso [Symmetric Group of the Genetic-Code Cubes.](https://is.gd/pTl8Js)
#' @import Biostrings
#' @importFrom GenomicRanges makeGRangesFromDataFrame GRanges GRangesList
#' @export
#' @author Robersy Sanchez <https://genomaths.com>
#' @references 
#' \enumerate{
#'  \item Robersy Sanchez, Jesús Barreto (2021) Genomic Abelian Finite 
#'   Groups.
#'  [doi: 10.1101/2021.06.01.446543](https://doi.org/10.1101/2021.06.01.446543).
#'  \item M. V José, E.R. Morgado, R. Sánchez, T. Govezensky, The 24 possible 
#'  algebraic representations of the standard genetic code in six or in three 
#'  dimensions, Adv. Stud. Biol. 4 (2012) 119–152.[PDF](https://is.gd/na9eap). 
#'  \item R. Sanchez. Symmetric Group of the Genetic–Code Cubes. Effect of the 
#'  Genetic–Code Architecture on the Evolutionary Process MATCH Commun. Math. 
#'  Comput. Chem. 79 (2018) 527-560.
#' }

#' @aliases codon_coord
setGeneric("codon_coord",
    function(
            codon,
            ...)
        standardGeneric("codon_coord"))


setClassUnion("DNAStringSet_OR_missing",
              c("DNAStringSet", "DNAMultipleAlignment", "missing"))

#' @aliases codon_coord
#' @rdname codon_coord
setMethod("codon_coord", signature(codon = "DNAStringSet_OR_missing"),
    function(
            codon = NULL,
            filepath = NULL,
            cube = c(
                    "ACGT","AGCT","TCGA","TGCA","CATG",
                    "CTAG","GATC","GTAC","ACTG","ATCG",
                    "GTCA","GCTA","CAGT","TAGC","TGAC",
                    "CGAT","AGTC","ATGC","CGTA","CTGA",
                    "GACT","GCAT","TACG","TCAG"), 
            group = c("Z4","Z64", "Z125", "R^3"),
            start = NA,
            end = NA,
            chr = 1L,
            strand = "+") {
    
        if (nchar(codon) %% 3 != 0) {
            stop("*** 'codon' argument is not base-triplet sequence.",
                 " A base-triplet sequence is multiple of 3.")
        }
        
        cube <- match.arg(cube)
        group <- match.arg(group)
        
        codon <- base_coord(
                            base = codon, 
                            filepath = filepath,
                            cube = cube,
                            group = group,
                            start = start,
                            end = end,
                            chr = chr,
                            strand = strand)
        
        codon <- mcols(codon)
        idx_seq <- grep("seq", colnames(codon))
        idx_coord <- grep("coord", colnames(codon))
        
        idx <- seq(1, nrow(codon) - 2, 3)
        codon <- lapply(idx, function(k) {
            s <- t(as.matrix(codon[k:(k + 2), idx_seq]))
            s <- apply(s, 1, paste, collapse = "")
            
            co <- t(as.matrix(codon[k:(k + 2), idx_coord]))
            co <- apply(co, 1, paste, collapse = "")
            return(c(s, co))
        })
        codon <- do.call(rbind, codon)
        
        pos <- seq(1, nrow(codon), 1L)
        codon <- data.frame(
            seqnames = chr, 
            start = pos, 
            end = pos,
            strand = strand, 
            codon)
        
        codon <- makeGRangesFromDataFrame(codon, keep.extra.columns = TRUE)
        return(codon)
    }
)

#' @aliases base_coord
#' @rdname codon_coord
#' @export
setGeneric("base_coord",
    function(
            base = NULL,
            filepath = NULL,
            ...)
        standardGeneric("base_coord"))


#' @aliases base_coord
#' @rdname codon_coord
#' @export
setMethod("base_coord", signature(base = "DNAStringSet_OR_missing"),
    function(
        base = NULL, 
        filepath = NULL,
        cube = c(
                "ACGT","AGCT","TCGA","TGCA","CATG",
                "CTAG","GATC","GTAC","ACTG","ATCG",
                "GTCA","GCTA","CAGT","TAGC","TGAC",
                "CGAT","AGTC","ATGC","CGTA","CTGA",
                "GACT","GCAT","TACG","TCAG"),
        group = c("Z4","Z5","R^3"),
        start = NA,
        end = NA,
        chr = 1L,
        strand = "+") {
        
        cube <- match.arg(cube)
        group <- match.arg(group)
        cube <- toupper(cube)

        if (is.null(base) && is.null(filepath)) 
            stop("*** Arguments 'base' & 'filepath' cannot be",
                 " simultaneously NULL.")
        
        if (!is.null(filepath) && is.character(filepath)) 
            base <- readDNAMultipleAlignment(filepath = filepath)
        
        if (inherits(base, "DNAMultipleAlignment"))
            base <- base@unmasked

        len <- min(width(base))
        if (!is.na(start) || !is.na(end)) {
            if (!is.na(start) && start > len)
                stop("*** The 'start' argument is greater than",
                     " the 'base' length")
            if (!is.na(end) && end > len)
                stop("*** The 'end' argument is greater than",
                     "   the 'base' length")
            base <- DNAStringSet(
                                base,
                                start = start,
                                end = end)
        }
        
        if (length(base) > 1) {
            base <- t(as.matrix(base))
            colnames(base) <- NULL
            base <- data.frame(base)
        }
        else {
            base <- as.character(base)
            base <- strsplit(base, "")[[1]]
            base <- data.frame(base)
        }
        seq <- base
        base <- base_repl(base = base, cube = cube, group = group)
        if (is.na(start)) 
            start <- 1L
        if (is.na(end)) 
            end <- len
        
        pos <- seq(start, end, 1L)
        if (!is.null(dim(base))) {
            colnames(base) <- paste0("coord", seq_len(ncol(base)))
            colnames(seq) <- paste0("seq", seq_len(ncol(base)))
        }
        base <- data.frame(
            seqnames = chr, 
            start = pos, 
            end = pos,
            strand = strand, 
            seq,
            base)

        base <- makeGRangesFromDataFrame(base, keep.extra.columns = TRUE)
        return(base)
    }
)



## --------------------------- Auxiliary functions --------------------------

base_repl <- function(base, cube, group) {
    alf <- strsplit(cube, "")[[1]]
    
    if (group == "Z4") {
        base[ base == "U" ] <- "T"
        base[ base == alf[1] ] <- 0
        base[ base == alf[2] ] <- 1
        base[ base == alf[3] ] <- 2
        base[ base == alf[4] ] <- 3
        base[ base == "-" ] <- -1
        base[ base == "N" ] <- -1
    }
    
    if (group == "Z5") {
        base[ base == "U" ] <- "T"
        base[ base == alf[1] ] <- 1
        base[ base == alf[2] ] <- 2
        base[ base == alf[3] ] <- 3
        base[ base == alf[4] ] <- 4
        base[ base == "-" ] <- 0
        base[ base == "N" ] <- 0
    }
    
    if (group == "R^3") {
        base[ base == "U" ] <- "T"
        base[ base == alf[1] ] <- -2
        base[ base == alf[2] ] <- -1
        base[ base == alf[3] ] <- 1
        base[ base == alf[4] ] <- 2
        base[ base == "-" ] <- 0
        base[ base == "N" ] <- 0
    }
    
    if (is.data.frame(base)) 
        base <- apply(base, 2, as.numeric)
    else 
        base <- as.numeric(base)
    return(base)
}

CodonCoordZ4toZ64 <-  function(x) 4 * x[1] + 16 * x[2] + x[3]
CodonCoordZ5toZ125 <-  function(x) 5 * x[1] + 25 * x[2] + x[3]
Codon2CoordR <- function(codon, cube) {
    base <- strsplit(cube, "")[[1]]
    mc <- t(sapply(codon, function(x) match(strsplit(x, "")[[1]], base)))
    mc[is.na(mc)] <- 5
    #        1   2  3  4  5
    R <- c(-2, -1, 1, 2, 0)
    mc <- t(apply(mc, 1, function(k) c(R[k[1]], R[k[2]], R[k[3]])))
    mc <- data.frame(mc)
    return(mc)
}

