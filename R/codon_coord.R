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
#' @title Codon coordinates on a given a given Abelian group representation.
#' @description Given a string denoting a codon or base from the DNA (or RNA) 
#' alphabet and a genetic-code Abelian group as given in reference (1).
#' @param codon,base An object from a \code{\link[Biostrings]{DNAStringSet}} or 
#' \code{\link[Biostrings]{DNAMultipleAlignment}} class carrying the DNA
#' pairwise alignment of two sequences. 
#' @param filepath A character vector containing the path to a file in 
#' \emph{\strong{fasta}} format to be read. This argument must be given if 
#' \emph{codon & base} arguments are not provided.
#' @param cube A character string denoting one of the 24 Genetic-code cubes,
#' as given in references (2 2 3).
#' @param start,end,chr,strand Optional parameters required to build a 
#' \code{\link[GenomicRanges]{GRanges-class}}. If not provided the default 
#' values given for the function definition will be used.
#' @param group A character string denoting the group representation for the 
#' given base or codon as shown in reference (1).
#' @param output Only used with \emph{get_coord} function. See \emph{details} 
#' section.
#' @param base_seq Logical. Only used with \emph{matrices} function. See 
#' \emph{details} section.
#' @param granges Only used with \emph{seqranges} function. See \emph{details}
#' section.
#' @details Symbols "-" and "N" usually found in DNA sequence alignments to 
#' denote gaps and missing/unknown bases are represented by the number: '-1' on
#' Z4 and '0' in Z5. In Z64 the symbol 'NA' will be returned for codons 
#' including symbols "-" and "N".
#' 
#' ## codon_coord:
#' This function returns a \code{\link[GenomicRanges]{GRanges-class}} object
#' carrying the codon sequence(s) and their respective coordinates in the 
#' requested Abelian group or simply, when \emph{group =  'R^3'} the coordinates
#' inserted in 3D space, which are derive from Z5 as indicated in reference (3).
#' Notice that the coordinates can be 3D or just one-dimension ("Z64" or
#' "Z125"). Hence, the pairwise alignment provided in argument 
#' \emph{\strong{codon}} must correspond to codon sequences. 
#' 
#' ## base_coord:
#' This function returns a \code{\link[GenomicRanges]{GRanges-class}} object 
#' carrying the DNA sequence(s) and their respective coordinates in the
#' requested Abelian group of base representation (one-dimension, "Z4" or "Z5").
#' Observe that to get coordinates in the set of of integer numbers ("Z") is
#' also possible but they are not defined to integrate a Abelian group. These
#' are just used for the further insertion the codon set in the 3D space (R^3).
#' 
#' ## get_coord:
#' Although the \code{\link[GenomicRanges]{GRanges-class}} object returned by
#' functions \code{\link{codon_coord}} and \code{\link{base_coord}} are useful
#' to store genomic information, the base and codon coordinates are not given on
#' them as numeric magnitudes. Function \code{\link{get_coord}} provides the way
#' to get the coordinates in a numeric object in object from and still to
#' preserve the base/codon sequence information. An object from
#' \code{\link{CodonSeq}} class is returned when \emph{output = "all"}. This has
#' two slots, the first one carrying a list of matrices and the second one
#' carrying the codon/base sequence information. That is, if \emph{x} is an
#' object from \code{\link{CodonSeq}} class, then a list of matrices of codon
#' coordinate can be retrieved as x@CoordList and the information on the codon
#' sequence as x@SeqRanges.
#' 
#' ## matrices & seqranges:
#' These are alternative ways to get the list of matrices of base/codon
#' coordinate and the information on the codon sequence from
#' \code{\link{CodonSeq}} and \code{\link{MatrixList}} class objects. These
#' functions can either take the output from functions \code{\link{base_coord}}
#' and \code{\link{codon_coord}}  or to operate directly on a
#' \code{\link[Biostrings]{DNAStringSet}} or to retrieve the a DNA sequence
#' alignment from a file.
#' 
#' \emph{\strong{base_seq}} parameter will determine whether to return the
#' matrices of coordinate for a DNA or codon sequence. While in function
#' \code{\link{seqranges}}, \emph{\strong{granges}} parameter will determine
#' whether to return a \code{\link[GenomicRanges]{GRanges-class}} object or a
#' \code{\link[S4Vectors]{DataFrame}}.
#' 
#' @seealso [Symmetric Group of the Genetic-Code Cubes.](
#' https://github.com/genomaths/GenomeAlgebra_SymmetricGroup)
#' @import Biostrings
#' @importFrom GenomicRanges makeGRangesFromDataFrame GRanges GRangesList
#' @importFrom S4Vectors mcols DataFrame
#' @importFrom Biostrings DNAStringSet
#' @importFrom methods new
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
#' 
#' @aliases codon_coord
setGeneric("codon_coord",
    function(
            codon = NULL,
            filepath = NULL,
            cube = "ACGT",
            group = "Z4",
            ...)
        standardGeneric("codon_coord"))


#' @aliases codon_coord
#' @rdname codon_coord
setMethod("codon_coord", signature(codon = "DNAStringSet_OR_NULL"),
    function(
            codon = NULL,
            filepath = NULL,
            cube = c(
                    "ACGT","AGCT","TCGA","TGCA","CATG",
                    "GTAC","CTAG","GATC","ACTG","ATCG",
                    "GTCA","GCTA","CAGT","TAGC","TGAC",
                    "CGAT","AGTC","ATGC","CGTA","CTGA",
                    "GACT","GCAT","TACG","TCAG"), 
            group = c("Z4","Z5", "Z64", "Z125", "Z4^3", "Z5^3"),
            start = NA,
            end = NA,
            chr = 1L,
            strand = "+") {
    
        if (!is.null(filepath) && is.character(filepath)) 
            codon <- readDNAMultipleAlignment(filepath = filepath)
        
        if (any(nchar(codon) %% 3 != 0)) {
            stop("*** 'codon' argument is not a base-triplet sequence.",
                 " A base-triplet sequence is multiple of 3.")
        }
        
        cube <- match.arg(cube)
        group <- match.arg(group)
        
        base_grp <- group
        if (is.element(group, "Z64")) 
            base_grp <- "Z4"
        if (is.element(group, "Z125")) 
            base_grp <- "Z5"
        if (is.element(group, "Z4^3")) 
            base_grp <- "Z4"
        if (is.element(group, "Z5^3")) 
            base_grp <- "Z5"
            
        codon <- base_coord(
                            base = codon, 
                            filepath = NULL,
                            cube = cube,
                            group = base_grp,
                            start = start,
                            end = end,
                            chr = chr,
                            strand = strand)
        
        if (length(codon) %% 3 != 0) {
            stop("*** 'codon' argument is not a base-triplet sequence.",
                 " A base-triplet sequence is multiple of 3.")
        }
        
        codon <- mcols(codon)
        idx_seq <- grep("seq", colnames(codon))
        idx_coord <- grep("coord", colnames(codon))
        
        sq <- data.frame(codon[, idx_seq])
        f <- factor(as.vector(sapply(seq_len(nrow(sq)/3), rep, 3)))
        sq <- split(sq, f)
        
        crd <- data.frame(codon[, idx_coord])
        crd <- split(crd, f)
        
        idx <- seq_along(sq)
        if (is.element(group, c("Z4","Z5","Z4^3","Z5^3"))) {
            codon <- lapply(idx, function(k) {
                c(apply(sq[[k]], 2, paste, collapse = ""),
                apply(crd[[k]], 2, paste, collapse = ","))
            })
        }
        else {
            fun <- switch(group,
                Z64 = CodonCoordZ4toZ64,
                Z125 = CodonCoordZ5toZ125,
            )
            
            codon <- lapply(idx, function(k) {
                c(apply(sq[[k]], 2, paste, collapse = ""),
                  apply(crd[[k]], 2, fun))
            })
        }

        rm(sq, crd); gc()
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
            cube = "ACGT",
            group = "Z4",
            ...)
        standardGeneric("base_coord"))


#' @aliases base_coord
#' @rdname codon_coord
#' @export
setMethod("base_coord", signature(base = "DNAStringSet_OR_NULL"),
    function(
        base = NULL, 
        filepath = NULL,
        cube = c(
            "ACGT","AGCT","TCGA","TGCA","CATG",
            "GTAC","CTAG","GATC","ACTG","ATCG",
            "GTCA","GCTA","CAGT","TAGC","TGAC",
            "CGAT","AGTC","ATGC","CGTA","CTGA",
            "GACT","GCAT","TACG","TCAG"), 
        group = c("Z4","Z5"),
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


#' @aliases get_coord
#' @rdname codon_coord
#' @export
setGeneric("get_coord",
    function(
            x,
            ...)
        standardGeneric("get_coord"))

#' @aliases get_coord
#' @rdname codon_coord
#' @export
setMethod("get_coord", signature(x = "GRanges"),
    function(
        x, 
        output = c("all", "matrix.list")) {
        
        output <- match.arg(output)
        
        gr <- x
        x <- mcols(x)
        nms <- colnames(x)
        idx <- grep("coord[[:digit:]]+", colnames(x))
        
        idx_seq <- grepl("seq[[:digit:]]+", colnames(x))
        if (any(idx_seq)) {
            gr <- gr[, which(idx_seq)]
            gr <- makeGRangesFromDataFrame(data.frame(gr), 
                                        keep.extra.columns = TRUE)
        }
        
        if (length(idx) == 0) 
            stop("Argument 'x' must carry at leat one named columns: coord1, ", 
                "coord2 & so on, as returned by function 'codon_coord'.")
        if (nchar(x[1, idx[1]]) >= 5) {
            m <- lapply(idx, function(k) {
                m <- do.call(rbind,strsplit(x[, k], ","))
                m <- apply(m,2, as.numeric)
            })
        }
        else {
            m <- lapply(idx, function(k) {
                m <- x[, k]
                as.numeric(m)
            })
        }
        names(m) <- nms[ idx ]
        res <- switch(output,
            "all" = new("CodonSeq", CoordList = m, SeqRanges = gr),
            "matrix.list" = new("MatrixList", matrices = m, names = nms)
        )
        
    }
)


#' @aliases get_coord
#' @rdname codon_coord
#' @export
setMethod("get_coord", signature(x = "DNAStringSet_OR_NULL"),
    function(
            x, 
            output = c("all", "matrix.list"),
            base_seq = TRUE,
            filepath = NULL,
            cube = "ACGT",
            group = "Z4",
            start = NA,
            end = NA,
            chr = 1L,
            strand = "+") {
        
        if (!is.null(filepath)) 
            x <- NULL
        if (base_seq) {
            x <- base_coord(
                    base = x,
                    filepath = filepath,
                    cube = cube,
                    group = group,
                    start = start,
                    end = end,
                    chr = chr,
                    strand = strand)
        }
        else {
            x <- codon_coord(
                    codon = x,
                    filepath = filepath,
                    cube = cube,
                    group = group,
                    start = start,
                    end = end,
                    chr = chr,
                    strand = strand)
        }
        x <- get_coord(x, output = output)
        return(x)
    }
)

#' @aliases matrices
#' @rdname codon_coord
setGeneric("matrices",
           function(
                x,
                ...)
        standardGeneric("matrices"))

#' @aliases matrices
#' @rdname codon_coord
#' @export
setMethod("matrices", signature(x = "CodonSeq_OR_MatrixList"),
    function(x) {
        if (inherits(x, "CodonSeq")) 
            x <- x@CoordList
        else
            x <- x@matrices
        return(x)
    }
)


#' @aliases matrices
#' @rdname codon_coord
#' @export
setMethod("matrices", signature(x = "DNAStringSet_OR_NULL"),
    function(
            x,
            base_seq = TRUE,
            filepath = NULL,
            cube = "ACGT",
            group = "Z4",
            start = NA,
            end = NA,
            chr = 1L,
            strand = "+") {
        
        if (!is.null(filepath)) 
            x <- NULL
        if (base_seq) {
            x <- base_coord(
                        base = x,
                        filepath = filepath,
                        cube = cube,
                        group = group,
                        start = start,
                        end = end,
                        chr = chr,
                        strand = strand)
        }
        else {
            x <- codon_coord(
                        codon = x,
                        filepath = filepath,
                        cube = cube,
                        group = group,
                        start = start,
                        end = end,
                        chr = chr,
                        strand = strand)
        }
        x <- get_coord(x, output = "matrix.list")
        return(matrices(x))
    }
)


#' @aliases seqranges
#' @rdname codon_coord
setGeneric("seqranges",
    function(
            x,
            ...)
        standardGeneric("seqranges"))

#' @aliases seqranges
#' @rdname codon_coord
#' @export
setMethod("seqranges", signature(x = "CodonSeq"),
    function(
        x,
        granges = TRUE) {
        x <- x@SeqRanges
        if (granges)
            x <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
        return(x)
    }
)


#' @aliases seqranges
#' @rdname codon_coord
#' @export
setMethod("seqranges", signature(x = "DNAStringSet_OR_NULL"),
    function(
            x,
            granges = TRUE,
            base_seq = TRUE,
            filepath = NULL,
            cube = "ACGT",
            group = "Z4",
            start = NA,
            end = NA,
            chr = 1L,
            strand = "+") {
        
        if (!is.null(filepath)) 
            x <- NULL
        if (base_seq) {
            x <- base_coord(
                        base = x,
                        filepath = filepath,
                        cube = cube,
                        group = group,
                        start = start,
                        end = end,
                        chr = chr,
                        strand = strand)
        }
        else {
            x <- codon_coord(
                        codon = x,
                        filepath = filepath,
                        cube = cube,
                        group = group,
                        start = start,
                        end = end,
                        chr = chr,
                        strand = strand)
        }
        x <- get_coord(x, output = "all")
        return(seqranges(x, granges = granges))
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
        base[ base == "-" ] <- NA
        base[ base == "N" ] <- NA
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
    
    if (is.data.frame(base)) 
        base <- apply(base, 2, as.numeric)
    else 
        base <- as.numeric(base)
    return(base)
}

CodonCoordZ4toZ64 <-  function(x) {
    if (any(is.na(x))) 
        res <- NA
    res <- 4 * x[1] + 16 * x[2] + x[3]
    return(res)
}
CodonCoordZ5toZ125 <-  function(x) 5 * x[1] + 25 * x[2] + x[3]


