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

#' @rdname allClasses
#' @title A class definition to store codon coordinates given in the Abelian
#' group and the codon sequence.
#' @description An objects from 'CodonSeq' or 'MatrixList' class is returned by
#' function \code{\link{get_coord}}. This object will store the coordinate of
#' each sequence in a list of 3D-vectors or a list of vectors located in the
#' slot named 'CoordList'. The original codon sequence (if provided) will be
#' stored in the slot named 'SeqRanges'.
 

#' @aliases CodonSeq
setClass("CodonSeq",
        slots = c(
            CoordList = "list",
            SeqRanges = "DataFrame_OR_NULL"
        )
)

#' @aliases MatrixList
#' @rdname allClasses
setClass("MatrixList",
        slots = c(
            matrices = "list",
            names = "character"
        )
)

setClassUnion("CodonSeq_OR_MatrixList", c("CodonSeq", "MatrixList"))


# ======================= Show method =================================

#' @rdname allClasses
#' @aliases show-CodonSeq
#' @title Show method for 'CodonSeq' or 'MatrixList' class object
#' @param object An object from 'CodonSeq' or 'MatrixList' class 
#' @importFrom methods show
#' @keywords internal
#' @export
setMethod(
    "show",
    signature = "CodonSeq",
    definition = function(object) {
        nams <- names(object@CoordList)
        cat(class(object)," object of length: ",
            length(object@CoordList), "\n", sep = "")
        cat(paste0("names(", length(object@CoordList), "):"), nams, "\n")
        cat("------- \n")
        print(.showMatrix(object@CoordList[[1]]))
        cat("...\n")
        cat("<",
            length(object@CoordList) - 1,
            " more ", class(object@CoordList[[1]])[1], " element(s)>\n",
            sep = "")
        cat("Two slots: 'CoordList' & 'SeqRanges'\n")
        cat("------- \n")
        invisible(object)
    }
)

#' @rdname allClasses
#' @aliases show-MatrixList
#' @title Show method for 'MatrixList' class object
#' @param object An object from 'MatrixList' class
#' @importFrom methods show
#' @keywords internal
#' @export
setMethod(
    "show",
    signature = "MatrixList",
    definition = function(object) {
        nams <- names(object@matrices)
        cat(class(object)," object of length: ",
            length(object@matrices), "\n", sep = "")
        cat(paste0("names(", length(object@matrices), "):"), nams, "\n")
        cat("------- \n")
        print(.showMatrix(object@matrices[[1]]))
        cat("...\n")
        cat("<",
            length(object@matrices) - 1,
            " more ", class(object@matrices[[1]])[1], " element(s)>\n",
            sep = "")
        cat("Two slots: 'matrices' & 'names'\n")
        cat("------- \n")
        invisible(object)
    }
)

.showMatrix <- function(x) {
    d <- dim(x)
    cat("Matrix with", d[1], "rows and", d[2], "columns:\n" )
    if (d[1] > 10) {
        r <- c()
        for (k in c(1:5, (d[1] - 5):d[1])) {
            r <- rbind(r, x[k,])
        }
        r[6, ] <- "..."
    } 
    r <- data.frame(r)
    rown <- paste0(c(1:5, (d[1] - 5):d[1]), ":")
    rown[6] <- "..."
    rownames(r) <- rown
    return(r)
}


