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
#' @aliases codon_coord
#' @aliases base_coord
#' @title Represent codon (DNA base triplet) with its coordinates in a given 
#' Abelian group
#' @description Given a string denoting a codon or base from the DNA (or RNA) 
#' alphabet and a genetic-code Abelian group as given in reference [#1].
#' @param codon A character string of DNA/RNA base-triplets (i.e., with
#' letter from the DNA/RNA alphabet: "A","C","G","T", and "U").
#' @param cube A character string denoting one of the 24 Genetic-code cubes,
#' as given in references [#2] and [#3].
#' @param group A character string denoting the group representation for the 
#' given base or codon as shown in reference [#1].
#' @details Symbols "-" and "N" usually found in DNA sequence alignments to 
#' denote gaps and missing/unknown bases are represented by the number: '-1' on
#' Z4 and '0' in Z5. In Z64 the symbol 'NA' will be returned for codons 
#' including symbols "-" and "N".
#' @seealso [Symmetric Group of the Genetic-Code Cubes.](https://is.gd/pTl8Js)
#' @export
#' @author Robersy Sanchez <https://genomaths.com>
#' @references 
#' \enumerate{
#'  \item Robersy Sanchez, Jesús Barreto (2021) Genomic Abelian Finite 
#'   Groups. #1
#'  [doi: 10.1101/2021.06.01.446543](https://doi.org/10.1101/2021.06.01.446543).
#'  \item M. V José, E.R. Morgado, R. Sánchez, T. Govezensky, The 24 possible 
#'  algebraic representations of the standard genetic code in six or in three 
#'  dimensions, Adv. Stud. Biol. 4 (2012) 119–152.[PDF](https://is.gd/na9eap). 
#'  #2
#'  \item R. Sanchez. Symmetric Group of the Genetic–Code Cubes. Effect of the 
#'  Genetic–Code Architecture on the Evolutionary Process MATCH Commun. Math. 
#'  Comput. Chem. 79 (2018) 527-560.#3
#' }
#' @examples 
#' ## Codon coordinate in R^3 (the genetic-code cube inserted in the 3D space)
#' codon_coord("ACG", group = "R^3")
#' 
#' codon_coord("ACG", group = "R^3", cube = "CUAG") 
#' 
#' ## Codon coordinate in Z4
#' codon_coord("ACG", group = "Z4") 
#' 
#' codon_coord("ACG", group = "Z4", cube = "UGCA") 
#' 
#' ## Codon coordinate in Z5
#' codon_coord("ACG", group = "Z5") 
#' 
#' codon_coord("ACG", group = "Z5", cube = "UGCA")


#' @aliases codon_coord
setGeneric("codon_coord",
    function(
            codon,
            ...)
        standardGeneric("codon_coord"))


#' @aliases codon_coord
#' @rdname codon_coord
setMethod("codon_coord", signature(codon = "character"),
    function(
            codon,
            cube = "ACGT", 
            group = c("Z4","Z64", "GF(5)", "Z125", "R^3")) {
    
    group <- match.arg(group)
    cube <- toupper(cube)
    codon <- toupper(codon)
    alf <- strsplit(cube, "")[[1]]
    if (any(!is.element(alf, c("A","C","G","T","U","-", "N"))))
        stop("*** Some codon bases does not belong to ",
             "the DNA/RNA base alphabet.")
    
    codon <- switch(group,
        Z4 = base_repl(codon = codon, cube = cube, group = "Z4"),
        Z64 = CodonCoordZ4toZ64(
                base_repl(
                        codon = codon, 
                        cube = cube,
                        group = "Z4")),
        Z5 = base_repl(codon = codon, cube = cube, group = "Z5"),
        Z125 = CodonCoordZ5toZ125(
                            base_repl(
                                codon = codon, 
                                cube = cube,
                                group = "Z5")),
        "R^3" = Codon2CoordR(codon = codon, 
                             cube = cube)
    )
    return(codon)
    }
)


## --------------------------- Auxiliary functions --------------------------

base_repl <- function(codon, cube, group) {
    codon <- toupper(codon)
    alf <- strsplit(cube, "")[[1]]
    codon <- strsplit(codon, "")[[1]]
    
    if (group == "Z4") {
        codon[ codon == alf[1] ] <- 0
        codon[ codon == alf[2] ] <- 1
        codon[ codon == alf[3] ] <- 2
        codon[ codon == alf[4] ] <- 3
        codon[ codon == "-" ] <- -1
        codon[ codon == "N" ] <- -1
    }
    
    if (group == "Z5") {
        codon[ codon == alf[1] ] <- 1
        codon[ codon == alf[2] ] <- 2
        codon[ codon == alf[3] ] <- 3
        codon[ codon == alf[4] ] <- 4
        codon[ codon == "-" ] <- 0
        codon[ codon == "N" ] <- 0
    }
    return(as.numeric(codon))
}

CodonCoordZ4toZ64 <-  function(x) 4 * x[1] + 16 * x[2] + x[3]
CodonCoordZ5toZ125 <-  function(x) 5 * x[1] + 25 * x[2] + x[3]
Codon2CoordR <- function(codon, cube) {
    base <- strsplit(cube, "")[[1]]
    mc <- t(sapply(codon, function(x) match(strsplit(x, "")[[1]], base)))
    mc[is.na(mc)] <- 5
    #       1   2  3  4  5
    R <- c(-2, -1, 1, 2, 0)
    mc <- t(apply(mc, 1, function(k) c(R[k[1]], R[k[2]], R[k[3]])))
    return(as.vector(mc))
}

