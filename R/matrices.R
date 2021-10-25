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

#' @rdname matrices
#' @title Codon coordinates on a given a given Abelian group representation.
#' @description Given a string denoting a codon or base from the DNA (or RNA) 
#' alphabet and a genetic-code Abelian group as given in reference (1).
#' @param x An object from a \code{\link[Biostrings]{DNAStringSet}} or 
#' \code{\link[Biostrings]{DNAMultipleAlignment}} class carrying the DNA
#' pairwise alignment of two sequences. 
#' @param base_seq Logical. Whether to return the base or codon coordinates on
#' the selected Abelian group. If codon coordinates are requested, then the 
#' number of the DNA bases in the given sequences must be multiple of 3.
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
#' @details 
#' These are alternative ways to get the list of matrices of base/codon
#' coordinate and the information on the codon sequence from
#' \code{\link{CodonSeq}} and \code{\link{MatrixList}} class objects. These
#' functions can either take the output from functions \code{\link{base_coord}}
#' and \code{\link{matrices}}  or to operate directly on a
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
#' @import GenomicRanges
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
#' @aliases matrices
setGeneric("matrices",
    function(
            x,
            ...)
    standardGeneric("matrices"))

#' @aliases matrices
#' @rdname matrices
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
#' @rdname matrices
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
            x <- matrices(
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
