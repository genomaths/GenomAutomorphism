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


#' @aliases automorphismByCoef
#' @rdname automorphismByCoef
#' @title Autmorphism Grouping by Coefficient
#' @description Automorphisms with the same automorphism's coefficients are 
#' grouped.
#' @return An \code{\link{AutomorphismByCoef}} class object.
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @export
#' @examples 
#' ## Load dataset
#' data(autm, package = "GenomAutomorphism")
#' 
#' automorphismByCoef(x = autm[1:2])
setGeneric("automorphismByCoef",
    function(
            x,
            ...)
            standardGeneric("automorphismByCoef"))

#' @aliases automorphismByCoef
#' @rdname automorphismByCoef
#' @param x An automorphism-class object returned by function 
#' \code{\link{automorphism}}.
#' @importFrom data.table data.table
#' @export
setMethod("automorphismByCoef", signature(x = "Automorphism"),
    function(x) {
        i <- 1
        l <- length(x)
        idx <- vector(mode = "numeric", length = length(x))
        coefs <- x$autm[1]
        for (k in seq_len(l)) {
            if ( x$autm[k] != coefs ) {
                i <- i + 1
                coefs <- x$autm[ k ]
            } 
            idx[ k ] <- i
        }
        
        x$idx <- factor(idx)
        x <- data.table(data.frame(x))
        x <- x[, list(
                    seqnames = unique(seqnames), start = min(start),
                    end = max(end), strand = unique(strand), 
                    autm = unique(autm),
                    cube = unique(cube)), 
                    by = idx ]
        x <- x[, c( "seqnames", "start", "end", 
                          "strand", "autm", "cube") ]
        x <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
        x <- sortByChromAndStart(x)
        return(as(x, "AutomorphismByCoef"))
    }
)



#' @aliases automorphismByCoef
#' @rdname automorphismByCoef
#' @param x An AutomorphismList-class object returned by function 
#' \code{\link{automorphism}}.
#' @param num.cores,tasks Integers. Argument \emph{num.cores} denotes the 
#' number of cores to use, i.e. at most how many child processes will be run
#' simultaneously (see \code{\link[BiocParallel]{bplapply}} function from
#' BiocParallel package). Argument \emph{tasks} denotes the number of tasks per
#' job. value must be a scalar integer >= 0L. In this documentation a job is
#' defined as a single call to a function, such as
#' \code{\link[BiocParallel]{bplapply}}. A task is the division of the \eqn{X}
#' argument into chunks. When tasks == 0 (default), \eqn{X} is divided as evenly
#' as possible over the number of workers (see
#' \code{\link[BiocParallel]{MulticoreParam}} from BiocParallel package).
#' @importFrom GenomicRanges GRangesList
#' @importFrom parallel detectCores
#' @importFrom BiocParallel MulticoreParam bplapply SnowParam
#' @importFrom data.table data.table
#' @export
setMethod("automorphismByCoef", signature(x = "AutomorphismList"),
    function(
            x,
            min.len = 1L,
            num.cores = detectCores(),
            tasks = 0L,
            verbose = TRUE) {
              
        gr <- try(x@SeqRanges, silent = TRUE)
        if (!inherits(gr, "try-error"))  
            x <- getAutomorphisms(x)
        
        x <- x@DataList
        
        ## ---------------- Setting parallel distribution --------- ##
        
        if (num.cores > 1) 
            num.cores <- num.cores - 1
        progressbar = FALSE
        if (verbose) progressbar = TRUE
        if (Sys.info()["sysname"] == "Linux")
            bpparam <- MulticoreParam(workers = num.cores, tasks = tasks,
                                      progressbar = progressbar)
        else bpparam <- SnowParam(workers = num.cores, type = "SOCK",
                                  progressbar = progressbar)
        
        ## ------------------------------------------------------- ##
        
        if (length(gr) > 0) {
            x <- lapply(x, function(x) {
                mcols(gr) <- x
                x <- automorphismByCoef(x)
                return(x)
            })
        }
        
        idx <- which(sapply(x, function(x) length(x) > min.len))
        x <- x[ idx ]
        
        return(as(x, "AutomorphismByCoefList"))
        return(x)
    }
)


