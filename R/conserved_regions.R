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
#' a \code{\link{AutomorphismByCoef}} or a \code{\link{AutomorphismByCoefList}} 
#' class object.
#' @param conserved Logical, Whether to return the *conserved* or the 
#' *non-conserved regions*.
#' @return A \code{\link{AutomorphismByCoef}} class object containing the 
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
    function(
            x, 
            conserved = TRUE,
            output = c("all_pairs", "unique_pairs", "unique")) {
        
        output <- match.arg(output)
        
        x <- automorphismByCoef(x)
        x <- conserved_regions(
                                x, 
                                conserved = conserved, 
                                output = output)
        return(x)
    }    
)


#' @rdname conserved_regions
#' @aliases conserved_regions
#' @importFrom GenomicRanges GRanges
#' @export
setMethod("conserved_regions", signature = "AutomorphismList",
    function(
            x, 
            conserved = TRUE,
            output = c("all_pairs", "unique_pairs", "unique")) {
        
        output <- match.arg(output)
        
        x <- automorphismByCoef(x)
        x <-  unlist(x)
        x <- conserved_regions(
                                x, 
                                conserved = conserved, 
                                output = output)
        return(x)
    }    
)

#' @rdname conserved_regions
#' @aliases conserved_regions
#' @importFrom GenomicRanges GRanges
#' @export
setMethod("conserved_regions", signature = "AutomorphismByCoef",
    function(
            x, 
            conserved = TRUE,
            output = c("all_pairs", "unique_pairs", "unique")) {
        
        output <- match.arg(output)
        
        if (conserved)
            x <- x[ x$autm == 1 ]
        else 
            x <- x[ x$autm != 1 ]
        x <- sortByChromAndEnd(x)
        
        x <- switch (output,
                all_pairs = x,
                unique_pairs = unique(x),
                unique = {
                    x <- unique(x)
                    x <- data.table(data.frame(x))
                    x <- x[, list(seqnames = unique(seqnames), 
                                  end = min(end), 
                                  strand = unique(strand), 
                                  autm = unique(autm)), 
                           by = c("start", "cube") ]
                    x <-makeGRangesFromDataFrame(x, 
                                                keep.extra.columns = TRUE)
                    x <- x[, c("autm", "cube")]
                }
        )
        return(x)
    }    
)

#' @rdname conserved_regions
#' @aliases conserved_regions
#' @importFrom GenomicRanges GRanges
#' @export
setMethod("conserved_regions", signature = "AutomorphismByCoefList",
    function(
            x, 
            conserved = TRUE,
            output = c("all_pairs", "unique_pairs", "unique")) {
        
        output <- match.arg(output)
        
        x <-  unlist(x)
        x <- conserved_regions(
                                x, 
                                conserved = conserved, 
                                output = output)
        return(x)
    }    
)


