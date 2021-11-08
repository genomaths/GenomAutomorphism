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

#' @rdname getAutomorphisms
#' @aliases getAutomorphisms
#' @title Get Automorphisms
#' @description This function returns an \code{\link{AutomorphismList}}-class 
#' object as a list of \code{\link{Automorphism}}-class objects, which inherits 
#' from \code{\link[GenomicRanges]{GRanges-class}} objects. 
#' @details For the sake of saving memory, each \code{\link{Automorphism}}-class 
#' objects is stored in an \code{\link{AutomorphismList}}-class, which  does 
#' not inherits from a \code{\link[GenomicRanges]{GRanges-class}}. This
#' function just transform each \code{\link{Automorphism}}-class object into
#' an object from the same class but now inheriting from a 
#' \code{\link[GenomicRanges]{GRanges-class}}. 
#' 
#' @importFrom S4Vectors mcols
#' @export
#' @examples 
#' ## Load dataset
#' data(autm, package = "GenomAutomorphism")
#' 
#' x1 <- autm[1:2]
#' x1
#' 
#' ## A list of DataFrame objects
#' as.list(x1)  
#' 
#' ## Get automorphism on GRanges objects
#' x1 <- getAutomorphisms(autm[1:2])
#' x1
#' 
#' ## A list of GRanges objects
#' as.list(x1)
#' 
setGeneric("getAutomorphisms",
    function(
        x,
        ...)
    standardGeneric("getAutomorphisms"))


#' @rdname getAutomorphisms
#' @aliases getAutomorphisms
#' @importFrom GenomicRanges GRanges
#' @export
setMethod("getAutomorphisms", signature = "AutomorphismList",
    function(x) {
        gr <- x@SeqRanges
        x <- x@DataList  
        if (length(gr) > 0) {
            x <- lapply(x, function(y) {
                    mcols(gr) <- y
                    x <- as(y, "Automorphism")
                    return(gr)
            })
            x <- new("AutomorphismList", 
                    DataList = x,
                    SeqRanges = GRanges()
            )
        } 
        else {
            x <- lapply(x, as, "Automorphism")
            x <- new("AutomorphismList", 
                    DataList = x,
                    SeqRanges = GRanges()
            )
        }
        return(x)
    }
)


#' @rdname getAutomorphisms
#' @aliases getAutomorphisms
#' @importFrom GenomicRanges GRanges
#' @export
setMethod("getAutomorphisms", signature = "list",
    function(x) {
        x <- as.AutomorphismList(x)
        x <- getAutomorphisms(x)
        return(x)
    }
)


#' @rdname getAutomorphisms
#' @aliases getAutomorphisms
#' @importFrom GenomicRanges GRanges
#' @export
setMethod("getAutomorphisms", signature = "DataFrame_OR_data.frame",
    function(x) {
        as(x, "Automorphism")       
    }
)

setMethod("[", "AutomorphismList", 
    function(x, i, ...) {
        x@DataList <- x@DataList[ i ]
        return(x)
    }
)

setMethod("[[", "AutomorphismList", 
    function(x, i, ...) {
        x <- getAutomorphisms(x)
        x <- x@DataList[[ i ]]
        return(x)
    }
)


