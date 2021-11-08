## Copyright (C) 2021 Robersy Sanchez <https://genomaths.com/>
## Author: Robersy Sanchez This file is part of the R package
## 'GenomAutomorphism'.  'GenomAutomorphism' is a free
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

#' @rdname conserved_regions
#' @aliases conserved_regions
#' @title Conserved and Non-conserved Regions from a MSA
#' @description Returns the Conserved or the Non-conserved Regions from a MSA.
#' @param x A \code{\link{Automorphism}}, a \code{\link{AutomorphismList}}, 
#' a \code{\link{AtomorphismByCoef}} or a \code{\link{AtomorphismByCoefList}} 
#' class object.
#' @param conserved Logical, Whether to return the *conserved* or the 
#' *non-conserved regions*.
#' @return A \code{\link{AtomorphismByCoef}} class object containing the 
#' requested regions.
#' @importFrom S4Vectors mcols
#' @export
#' @examples 
#' ## Load dataset
#' data(autm, package = "GenomAutomorphism")
#' conserved_regions(autm[1:3])
setGeneric("conserved_regions",
           function(
               x,
               ...)
               standardGeneric("conserved_regions"))


#' @rdname conserved_regions
#' @aliases conserved_regions
#' @importFrom GenomicRanges GRanges
#' @export
setMethod("conserved_regions", signature = "Automorphism",
    function(x, conserved = TRUE) {
        x <- automorphismByCoef(x)
        if (conserved)
            x <- x[ x$autm == 1 ]
        else 
            x <- x[ x$autm != 1 ]
        x <- sortByChromAndStart(x)
        return(x)
    }    
)


#' @rdname conserved_regions
#' @aliases conserved_regions
#' @importFrom GenomicRanges GRanges
#' @export
setMethod("conserved_regions", signature = "AutomorphismList",
    function(x, conserved = TRUE) {
        x <- automorphismByCoef(x)
        x <-  unlist(x)
        if (conserved)
            x <- x[ x$autm == 1 ]
        else 
            x <- x[ x$autm != 1 ]
        x <- sortByChromAndStart(x)
        return(x)
    }    
)

#' @rdname conserved_regions
#' @aliases conserved_regions
#' @importFrom GenomicRanges GRanges
#' @export
setMethod("conserved_regions", signature = "AtomorphismByCoef",
    function(x, conserved = TRUE) {
        if (conserved)
            x <- x[ x$autm == 1 ]
        else 
            x <- x[ x$autm != 1 ]
        x <- sortByChromAndStart(x)
        return(x)
    }    
)

#' @rdname conserved_regions
#' @aliases conserved_regions
#' @importFrom GenomicRanges GRanges
#' @export
setMethod("conserved_regions", signature = "AtomorphismByCoefList",
    function(x, conserved = TRUE) {
        x <-  unlist(x)
        if (conserved)
            x <- x[ x$autm == 1 ]
        else 
            x <- x[ x$autm != 1 ]
        x <- sortByChromAndStart(x)
        return(x)
    }    
)


