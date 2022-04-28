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
## 
#' An S4 class to extract elements from AutomorphismList-class object.
#' @rdname extract-methods
#' @aliases [
#' @aliases AutomorphismList-methods
#' @param x An \code{\link{AutomorphismList-class}} object.
#' @param i,j indices specifying elements to extract.
#' @description First and second level subsetting of 'x'. Extraction using names
#' can be done as x$name.
#' @return An element of x, an \code{\link{AutomorphismList-class}} object.
#' @export
#' @author Robersy Sanchez <https://genomaths.com>
setMethod("[", signature(x = "AutomorphismList", i = "integer"), 
    function(x, i, ...) {
        x@DataList <- x@DataList[ i ]
        return(x)
    }
)

#' @rdname extract-methods
#' @aliases [[
#' @aliases AutomorphismList-methods
#' @param x An \code{\link{AutomorphismList-class}} object.
#' @param i An integer to subsetting 'x'.
#' @description Second level subsetting of 'x'.
#' @return An element of x, an \code{\link{Automorphism-class}} object.
#' @exportMethod "[["
#' @export
setMethod("[[", signature(x = "AutomorphismList", i = "ANY"), 
    function(x, i, ...) {
        x <- x[ i ]
        x <- getAutomorphisms(x)
        x <- as(x, "GRangesList")
        return(x[[1]])
    }
)


#' @rdname extract-methods
#' @aliases $
#' @aliases AutomorphismList-methods
#' @param x An \code{\link{AutomorphismList-class}} object
#' @param name A literal character string naming an element from 'x'.
#' @description Subsetting of 'x' by element name.
#' @return An element of x, an \code{\link{Automorphism-class}} object.
#' @exportMethod "$"
#' @export
setMethod("$", signature(x = "AutomorphismList"), 
    function(x, name) {
        x@DataList <- x@DataList[ match(name, names(x)) ]
        x <- getAutomorphisms(x)
        x <- as(x, "GRangesList")
        return(x[[1]])    
    }
)
