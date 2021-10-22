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

#' @rdname automorphism
#' @title Compute the Automorphisms of Mutational Events Between two Codon 
#' Sequences Represented in a Given Abelian group.
#' @description 
#' ## automorphism:
#' Given two codon sequences represented in a given Abelian group,
#' this function computes the automorphisms describing codon mutational events.
#' Basically, this function is a wrapping to call the corresponding function 
#' for a specified Abelian group.
#' 
#' ## automorphismByRanges:
#' Herein, automorphisms are algebraic descriptions of mutational event
#' observed in codon sequences represented on different Abelian groups. In
#' particular, as described in references (3-4), for each representation of the
#' codon set on a defined Abelian group there are 24 possible isomorphic Abelian
#' groups. These Abelian groups can be labeled based on the DNA base-order used
#' to generate them. The set of 24 Abelian groups can be described as a group
#' isomorphic to the symmetric group of degree four (\eqn{S_4}, see reference
#' (4)). Function \code{\link{automorphismByRanges}} permits the classification
#' of the pairwise alignment of protein-coding sub-regions based on the
#' mutational events observed on it and on the genetic-code cubes that describe
#' them. 
#' 
#' @details 
#' ## automorphism:
#' Automorphisms in Z5, Z64 and Z125 are described as functions \eqn{f(x) =
#' k x mod 64} and \eqn{f(x) = k x mod 125}, where k and x are elements from the
#' set of integers modulo 64 or modulo 125, respectively. If an automorphism
#' cannot be found on any of the cubes provided in the argument \eqn{cube}, then
#' function \code{\link{automorphism}} will search for automorphisms in the
#' cubes provided in the argument \eqn{cube_alt}.
#' 
#' Automorphisms in Z5^3' are described as functions \eqn{f(x) = Ax mod Z5}, 
#' where A is diagonal matrix.
#' 
#' ## automorphismByRanges:
#' Arguments \emph{\strong{cube}} and \emph{\strong{cube_alt}} must be pairs of
#' dual cubes (see section 2.4 from reference 4).
#' 
#' 
#' @param seq An object from a \code{\link[Biostrings]{DNAStringSet}} or 
#' \code{\link[Biostrings]{DNAMultipleAlignment}} class carrying the DNA
#' pairwise alignment of two sequences. The pairwise alignment provided in
#' argument \emph{\strong{seq}} or the 'fasta' file \emph{\strong{filepath}}
#' must correspond to codon sequences.
#' @param filepath A character vector containing the path to a file in 
#' \emph{\strong{fasta}} format to be read. This argument must be given if 
#' \emph{codon & base} arguments are not provided.
#' @param cube,cube_alt A character string denoting pairs of the 24 Genetic-code
#' cubes, as given in references (2-3). That is, the base pairs from the given 
#' cubes must be complementary each other. Such a cube pair are call 
#' \eqn{dual cubes} and, as shown in reference (3), each pair integrates group.
#' @param start,end,chr,strand Optional parameters required to build a 
#' \code{\link[GenomicRanges]{GRanges-class}}. If not provided the default 
#' values given for the function definition will be used.
#' @param group A character string denoting the group representation for the 
#' given base or codon as shown in reference (1).
#' @return 
#' ## automorphism:
#' This function returns a \code{\link{Automorphism}} class object  with four
#' columns on its metacolumn named: \emph{seq1}, \emph{seq2}, \emph{autm}, and
#' \emph{cube}.
#' 
#' ## automorphismByRanges:
#' This function returns a \code{\link[GenomicRanges]{GRanges-class}} object.
#' Consecutive mutational events (on the codon sequence) described by 
#' automorphisms on a same cube are grouped in a range.
#' 
#' @seealso \code{\link{autZ64}}.
#' @importFrom numbers modlin
#' @export
#' @references 
#' \enumerate{
#'  \item Sanchez R, Morgado E, Grau R. Gene algebra from a genetic code 
#'  algebraic structure. J Math Biol. 2005 Oct;51(4):431-57. 
#'  doi: 10.1007/s00285-005-0332-8. Epub 2005 Jul 13. PMID: 16012800. (
#'  [PDF](https://arxiv.org/pdf/q-bio/0412033.pdf)).
#'  \item Robersy Sanchez, Jesús Barreto (2021) Genomic Abelian Finite 
#'   Groups.
#'  [doi: 10.1101/2021.06.01.446543](https://doi.org/10.1101/2021.06.01.446543).
#'  \item M. V José, E.R. Morgado, R. Sánchez, T. Govezensky, The 24 possible 
#'  algebraic representations of the standard genetic code in six or in three 
#'  dimensions, Adv. Stud. Biol. 4 (2012) 119–152.[PDF](https://is.gd/na9eap). 
#'  \item R. Sanchez. Symmetric Group of the Genetic–Code Cubes. Effect of the 
#'  Genetic–Code Architecture on the Evolutionary Process MATCH Commun. Math. 
#'  Comput. Chem. 79 (2018) 527-560. [PDF](https://bit.ly/2Z9mjM7)
#' }

#' @aliases automorphism
#' @export
setGeneric("automorphism",
    function(
            seq = NULL,
            filepath = NULL,
            group = "Z4",
            ...)
        standardGeneric("automorphism"))

#' @aliases automorphism
#' @rdname automorphism
#' @export
setMethod("automorphism", signature(seq = "DNAStringSet_OR_NULL"),
    function( 
            seq = NULL,
            filepath = NULL,
            group = c("Z4","Z5", "Z64", "Z125", "Z4^3", "Z5^3"),
            cube = c("ACGT", "TGCA"),
            cube_alt = c("CATG", "GTAC"),
            start = NA,
            end = NA,
            chr = 1L,
            strand = "+") {
        
        group <- match.arg(group)
        
        seq <- switch (group,
            "Z64" = autZ64(seq = seq,
                          filepath = filepath,
                          cube = cube,
                          cube_alt = cube_alt,
                          start = start,
                          end = end,
                          chr = chr,
                          strand = strand),
            "Z5^3" = aut3D(seq = seq,
                            filepath = filepath,
                            cube = cube,
                            group = "Z5^3",
                            cube_alt = cube_alt,
                            start = start,
                            end = end,
                            chr = chr,
                            strand = strand)
        )
        return(seq)
    }
)


#' @aliases automorphismByRanges
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @export
setGeneric("automorphismByRanges",
    function(
            autm,
            ...)
    standardGeneric("automorphismByRanges"))


#' @aliases automorphismByRanges
#' @rdname automorphism
#' @param autm An Automorphism-class object returned by function 
#' \code{\link{automorphism}}.
#' @importFrom data.table data.table
#' @export
setMethod("automorphismByRanges", signature(autm = "Automorphism"),
    function(autm) {
        i <- 1
        l <- length(autm)
        idx <- vector(mode = "numeric", length = length(autm))
        for (k in seq_len(l)) {
            if (autm$seq1[k] != autm$seq2[k])
                i <- i + 1
            idx[ k ] <- i
        }
        
        autm$idx <- factor(idx)
        autm <- data.table(data.frame(autm))
        autm <- autm[, list(
            seqnames = unique(seqnames), start = min(start),
            end = max(end), strand = unique(strand), 
            cube = unique(cube)), 
            by = idx ]
        autm <- autm[, c("seqnames", "start", "end", "strand", "cube")]
        autm <- makeGRangesFromDataFrame(autm, keep.extra.columns = TRUE)
        return(autm)
    }
)





