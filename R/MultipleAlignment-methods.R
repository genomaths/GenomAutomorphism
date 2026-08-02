## S4 methods for MultipleAlignment::DNAMultipleAlignment inputs.
##
## Bundled .rda objects may carry a legacy class definition that does not
## dispatch through the DNAStringSet_OR_NULL union after installation.
## These methods unwrap @unmasked and delegate to the union methods.

#' @importClassesFrom MultipleAlignment DNAMultipleAlignment
#' @importFrom methods callGeneric
NULL

.dma_unmasked <- function(x) {
    x@unmasked
}

#' @aliases base_coord
#' @rdname base_methods
setMethod(
    "base_coord", signature(base = "DNAMultipleAlignment"),
    function(base = NULL, filepath = NULL, ...) {
        callGeneric(base = .dma_unmasked(base), filepath = filepath, ...)
    }
)

#' @aliases seq2granges
#' @rdname base_methods
setMethod(
    "seq2granges", signature(base = "DNAMultipleAlignment"),
    function(base = NULL, filepath = NULL, ...) {
        callGeneric(base = .dma_unmasked(base), filepath = filepath, ...)
    }
)

#' @aliases base_matrix
#' @rdname base_methods
setMethod(
    "base_matrix", signature(base = "DNAMultipleAlignment"),
    function(base = NULL, ...) {
        callGeneric(base = .dma_unmasked(base), ...)
    }
)

#' @aliases codon_coord
#' @rdname codon_coord
setMethod(
    "codon_coord", signature(codon = "DNAMultipleAlignment"),
    function(codon = NULL, filepath = NULL, ...) {
        callGeneric(codon = .dma_unmasked(codon), filepath = filepath, ...)
    }
)

#' @aliases get_coord
#' @rdname get_coord
setMethod(
    "get_coord", signature(x = "DNAMultipleAlignment"),
    function(x, ...) {
        callGeneric(x = .dma_unmasked(x), ...)
    }
)

#' @aliases matrices
#' @rdname matrices
setMethod(
    "matrices", signature(x = "DNAMultipleAlignment"),
    function(x, ...) {
        callGeneric(x = .dma_unmasked(x), ...)
    }
)

#' @aliases seqranges
#' @rdname seqranges
setMethod(
    "seqranges", signature(x = "DNAMultipleAlignment"),
    function(x, ...) {
        callGeneric(x = .dma_unmasked(x), ...)
    }
)

#' @aliases automorphisms
#' @rdname automorphisms
setMethod(
    "automorphisms", signature(seqs = "DNAMultipleAlignment"),
    function(seqs = NULL, filepath = NULL, ...) {
        callGeneric(seqs = .dma_unmasked(seqs), filepath = filepath, ...)
    }
)
