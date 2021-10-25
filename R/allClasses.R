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

# ========================== BaseGroup ============================= #

#' @aliases BaseGroup
#' @rdname BaseGroup
#' @title A class definition to store codon automorphisms in given in the 
#' Abelian group representation.
#' @seealso \code{\link{automorphism}}
#' @export
setClass("BaseGroup",
         slots = c(
             seqnames = "Rle",
             ranges = "IRanges_OR_IPos",
             strand = "Rle",
             elementMetadata = "DataFrame",
             seqinfo = "Seqinfo",
             colnames = "character"
         ),
         contains = "GRanges"
)

# ====================  Validity BaseGroup ======================== #
#' @rdname BaseGroup
#' @title Valid BaseGroup mcols
#' @param x A 'BaseGroup' object
#' @keywords internal
valid.BaseGroup.mcols <- function(x) {
    if (length(x) > 0) {
        coln <- x@colnames
        if (any(!is.element(coln,
                            c("seq1", "seq2", "coord1", "coord2")))) {
            return("*** This is not a valid  BaseGroup-class object.",
                   "Columns from the matacolumn have the wrong names")
        }
        if (unique(nchar(x$seq1)) != 1 || unique(nchar(x$seq2)) != 1) 
            stop("*** This is not a valid  BaseGroup-class object.",
                 "seq1 or seq2 columns is not a base sequence") 
    }
    NULL
}

#' @rdname BaseGroup
#' @title Valid 'BaseGroup' inheritance from 'GRanges' class
#' @param x A 'BaseGroup object'
#' @keywords internal
valid.GRanges <- function(x) {
    if (length(x) > 0) {
        if (!inherits(x, "GRanges")) {
            return("*** This is not a valid  Automorphism-class object.")
        }
    }
    NULL
}

#' @rdname valid.BaseGroup
#' @title Valid BaseGroup
#' @param x A 'BaseGroup object'
#' @keywords internal
valid.BaseGroup <- function(x) 
    c(valid.GRanges(x), valid.BaseGroup.mcols(x))

S4Vectors:::setValidity2("BaseGroup", valid.BaseGroup)


# ========================== CodonGroup ============================= #

#' @aliases CodonGroup
#' @rdname CodonGroup
#' @title A class definition to store codon automorphisms in given in the 
#' Abelian group representation.
#' @seealso \code{\link{automorphism}}
#' @export
setClass("CodonGroup",
    slots = c(
            seqnames = "Rle",
            ranges = "IRanges_OR_IPos",
            strand = "Rle",
            elementMetadata = "DataFrame",
            seqinfo = "Seqinfo",
            colnames = "character"
    ),
    contains = "GRanges"
)

# ====================  Validity CodonGroup ======================== #
#' @rdname CodonGroup
#' @title Valid CodonGroup mcols
#' @param x A 'CodonGroup' object
#' @keywords internal
valid.CodonGroup.mcols <- function(x) {
    if (length(x) > 0) {
        coln <- x@colnames
        if (any(!is.element(coln,
                            c("seq1", "seq2", "coord1", "coord2")))) {
            return("*** This is not a valid  CodonGroup-class object.",
                "Columns from the matacolumn have the wrong names")
        }
        if (unique(nchar(x$seq1)) != 3 || unique(nchar(x$seq2)) != 3) 
            stop("*** This is not a valid  BaseGroup-class object.",
                "seq1 or seq2 columns is not a base-triplet sequence") 
    }
    NULL
}

#' @rdname valid.CodonGroup
#' @title Valid CodonGroup
#' @param x A 'CodonGroup object'
#' @keywords internal
valid.CodonGroup <- function(x) 
    c(valid.GRanges(x), valid.CodonGroup.mcols(x))

S4Vectors:::setValidity2("CodonGroup", valid.CodonGroup)


setClassUnion("BaseGroup_OR_CodonGroup", c("BaseGroup","CodonGroup"))


# ========================== CodonSeq ============================= #

#' @rdname CodonSeq
#' @title A class definition to store codon coordinates given in the Abelian
#' group and the codon sequence.
#' @description An objects from 'CodonSeq' or 'MatrixList' class is returned by
#' function \code{\link{get_coord}}. This object will store the coordinate of
#' each sequence in a list of 3D-vectors or a list of vectors located in the
#' slot named 'CoordList'. The original codon sequence (if provided) will be
#' stored in the slot named 'SeqRanges'.
#' @import GenomicRanges
#' @importFrom S4Vectors setValidity2
#' @importFrom methods validObject
 
#' @aliases CodonSeq
#' @export
setClass("CodonSeq",
        slots = c(
            CoordList = "list",
            SeqRanges = "GenomicRanges_OR_missing"
        )
)

#' @aliases MatrixList
#' @rdname allClasses
#' @export
setClass("MatrixList",
        slots = c(
            matrices = "list",
            names = "character"
        )
)

setClassUnion("CodonSeq_OR_MatrixList", c("CodonSeq", "MatrixList"))
setClassUnion("DNAStringSet_OR_NULL",
              c("DNAStringSet", "DNAMultipleAlignment", "NULL", "missing"))


#' @aliases Automorphism
#' @rdname Automorphism
#' @title A class definition to store codon automorphisms in given in the 
#' Abelian group representation.
#' @seealso \code{\link{automorphism}}
#' @export
setClass("Automorphism",
    slots = c(
        seqnames = "Rle",
        ranges = "IRanges_OR_IPos",
        strand = "Rle",
        elementMetadata = "DataFrame",
        seqinfo = "Seqinfo",
        colnames = "character"
    ),
    contains = "GRanges"
)

# ======================== Validity CodonSeq ======================= #
#' @rdname Automorphism
#' @title Valid Automorphism mcols
#' @param x A 'Automorphism object'
#' @keywords internal
valid.Automorphism.mcols <- function(x) {
    if (length(x) > 0) {
        coln <- x@colnames
        if (any(!is.element(coln,
                c("seq1", "seq2", "coord1", "coord2", "autm", "cube")))) {
            return("*** This is not a valid  Automorphism-class object.",
                "Columns from the matacolumn have the wrong names")
        }
    }
    NULL
}

#' @rdname valid.Automorphism
#' @title Valid Automorphism
#' @param x A 'Automorphism object'
#' @keywords internal
valid.Automorphism <- function(x) 
    c(valid.GRanges(x), valid.Automorphism.mcols(x))

S4Vectors:::setValidity2("Automorphism", valid.Automorphism)


# ======================= Show method =================================


#' @rdname allClasses
#' @aliases show-CodonSeq
#' @title Show method for 'CodonSeq' class object
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

# ========================== MatrixList ============================= #

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
    if (!is.null(d)) {
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
    }
    else {
        l <- length(x)
        cat("Vector of length:", l, "\n")
        if (l > 10) {
            r <- x[ c(1:5, (l - 5):l) ]
            r[6] <- "..."
            r <- data.frame(matrix(r, ncol = length(r)))
            nms <- paste0(c(1:5, (l - 5):l), ":")
            nms[6] <- r[6] 
            colnames(r) <- nms
        }
    }
    return(r)
}


