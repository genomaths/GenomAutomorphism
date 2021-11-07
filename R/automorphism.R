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

#' @rdname automorphism
#' @title Compute the Automorphisms of Mutational Events Between two Codon 
#' Sequences Represented in a Given Abelian group.
#' @description 
#' ## automorphism:
#' Given two codon sequences represented in a given Abelian group,
#' this function computes the automorphisms describing codon mutational events.
#' Basically, this function is a wrapping to call the corresponding function 
#' for a specified Abelian group.
#' 
#' ## automorphismByRanges:
#' Herein, automorphisms are algebraic descriptions of mutational event
#' observed in codon sequences represented on different Abelian groups. In
#' particular, as described in references (3-4), for each representation of the
#' codon set on a defined Abelian group there are 24 possible isomorphic Abelian
#' groups. These Abelian groups can be labeled based on the DNA base-order used
#' to generate them. The set of 24 Abelian groups can be described as a group
#' isomorphic to the symmetric group of degree four (\eqn{S_4}, see reference
#' (4)). Function \code{\link{automorphismByRanges}} permits the classification
#' of the pairwise alignment of protein-coding sub-regions based on the
#' mutational events observed on it and on the genetic-code cubes that describe
#' them. 
#' 
#' @details 
#' ## automorphism:
#' Automorphisms in Z5, Z64 and Z125 are described as functions \eqn{f(x) =
#' k x mod 64} and \eqn{f(x) = k x mod 125}, where k and x are elements from the
#' set of integers modulo 64 or modulo 125, respectively. If an automorphism
#' cannot be found on any of the cubes provided in the argument \eqn{cube}, then
#' function \code{\link{automorphism}} will search for automorphisms in the
#' cubes provided in the argument \eqn{cube_alt}.
#' 
#' Automorphisms in Z5^3' are described as functions \eqn{f(x) = Ax mod Z5}, 
#' where A is diagonal matrix.
#' 
#' ## automorphismByRanges:
#' Arguments \emph{\strong{cube}} and \emph{\strong{cube_alt}} must be pairs of
#' dual cubes (see section 2.4 from reference 4).
#' 
#' 
#' @param seq An object from a \code{\link[Biostrings]{DNAStringSet}} or 
#' \code{\link[Biostrings]{DNAMultipleAlignment}} class carrying the DNA
#' pairwise alignment of two sequences. The pairwise alignment provided in
#' argument \emph{\strong{seq}} or the 'fasta' file \emph{\strong{filepath}}
#' must correspond to codon sequences.
#' @param filepath A character vector containing the path to a file in 
#' \emph{\strong{fasta}} format to be read. This argument must be given if 
#' \emph{codon & base} arguments are not provided.
#' @param cube,cube_alt A character string denoting pairs of the 24 Genetic-code
#' cubes, as given in references (2-3). That is, the base pairs from the given 
#' cubes must be complementary each other. Such a cube pair are call 
#' \eqn{dual cubes} and, as shown in reference (3), each pair integrates group.
#' @param nms Optional. Only used if the DNA sequence alignment provided carries
#' more than two sequences. A character string giving short names for the
#' alignments to be compared. If not given then the automorphisms between 
#' pairwise alignment are named as: "aln_1", "aln_2", and so on.
#' @param start,end,chr,strand Optional parameters required to build a 
#' \code{\link[GenomicRanges]{GRanges-class}}. If not provided the default 
#' values given for the function definition will be used.
#' @param group A character string denoting the group representation for the 
#' given base or codon as shown in reference (1).
#' @return 
#' ## automorphism:
#' This function returns a \code{\link{Automorphism}} class object  with four
#' columns on its metacolumn named: \emph{seq1}, \emph{seq2}, \emph{autm}, and
#' \emph{cube}.
#' 
#' ## automorphismByRanges:
#' This function returns a \code{\link[GenomicRanges]{GRanges-class}} object.
#' Consecutive mutational events (on the codon sequence) described by 
#' automorphisms on a same cube are grouped in a range.
#' 
#' ## automorphismByCoef
#' This function returns a \code{\link[GenomicRanges]{GRanges-class}} object.
#' Consecutive mutational events (on the codon sequence) described by 
#' the same automorphism coefficients are grouped in a range.
#' 
#' ## getAutomorphisms
#' This function returns an AutomorphismList-class object as a list of
#' Automorphism-class objects, which inherits from
#' \code{\link[GenomicRanges]{GRanges-class}} objects.
#' 
#' @seealso \code{\link{autZ64}}.
#' @importFrom numbers modlin
#' @export
#' @references 
#' \enumerate{
#'  \item Sanchez R, Morgado E, Grau R. Gene algebra from a genetic code 
#'  algebraic structure. J Math Biol. 2005 Oct;51(4):431-57. 
#'  doi: 10.1007/s00285-005-0332-8. Epub 2005 Jul 13. PMID: 16012800. (
#'  [PDF](https://arxiv.org/pdf/q-bio/0412033.pdf)).
#'  \item Robersy Sanchez, Jesús Barreto (2021) Genomic Abelian Finite 
#'   Groups.
#'  [doi: 10.1101/2021.06.01.446543](https://doi.org/10.1101/2021.06.01.446543).
#'  \item M. V José, E.R. Morgado, R. Sánchez, T. Govezensky, The 24 possible 
#'  algebraic representations of the standard genetic code in six or in three 
#'  dimensions, Adv. Stud. Biol. 4 (2012) 119–152.[PDF](https://is.gd/na9eap). 
#'  \item R. Sanchez. Symmetric Group of the Genetic–Code Cubes. Effect of the 
#'  Genetic–Code Architecture on the Evolutionary Process MATCH Commun. Math. 
#'  Comput. Chem. 79 (2018) 527-560. [PDF](https://bit.ly/2Z9mjM7)
#' }
#' @examples 
#' ## Load a pairwise alignment
#' data(aln)
#' aln
#' 
#' ## Automorphism on "Z5^3"
#' autms <- automorphism(seq = aln, group = "Z5^3")
#' autms
#' 
#' ## Automorphism on "Z64"
#' autms <- automorphism(seq = aln, group = "Z64")
#' autms
#' 
#' ## Grouping into ranges the automorphisms by cubes 
#' automorphismByRanges(autms)
#' 
#' ## Automorphism on "Z64" from position 1 to 33
#' autms <- automorphism(seq = aln,
#'                       group = "Z64",
#'                       start = 1,
#'                       end = 33)
#' autms
#' 
#' ## Grouping into ranges the automorphisms by cubes 
#' automorphismByRanges(autms)
#' 
#' @aliases automorphism
#' @export
setGeneric("automorphism",
    function(
            seqs = NULL,
            filepath = NULL,
            group = "Z4",
            ...)
        standardGeneric("automorphism"))

#' @aliases automorphism
#' @rdname automorphism
#' @importFrom foreach foreach %dopar% %:%
#' @importFrom doParallel registerDoParallel
#' @importFrom stats setNames
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom BiocGenerics width 
#' @import Biostrings
#' @export
setMethod("automorphism", signature(seqs = "DNAStringSet_OR_NULL"),
    function( 
            seqs = NULL,
            filepath = NULL,
            group = c("Z5", "Z64", "Z125", "Z5^3"),
            cube = c("ACGT", "TGCA"),
            cube_alt = c("CATG", "GTAC"),
            nms = NULL,
            start = NA,
            end = NA,
            chr = 1L,
            strand = "+") {
        
        if (is.null(seqs) && is.null(filepath)) 
            stop("*** One of the arguments 'seqs' or 'filepath' must
                be provided.")
        
        group <- match.arg(group)
        
        if (!is.null(filepath) && is.character(filepath)) 
            seqs <- readDNAMultipleAlignment(filepath = filepath)
        
        if (inherits(seqs, "DNAStringSet"))
            nr <- length(seqs)
        
        if (inherits(seqs, "DNAMultipleAlignment"))
            nr <- nrow(seqs)
        
        if (nr < 3) {
            seqs <- selectAutomorphism(
                                    seq = seqs, 
                                    filepath = filepath, 
                                    group = group, 
                                    cube = cube, 
                                    cube_alt = cube_alt, 
                                    start = start, 
                                    end = end, 
                                    chr = chr, 
                                    strand = strand)
        } 
        else {
            
            gr <- seqs@unmasked[ c(1, 1) ]
            gr <- selectAutomorphism(
                seq = gr, 
                filepath = NULL, 
                group = group, 
                cube = cube, 
                cube_alt = cube_alt, 
                start = start, 
                end = end, 
                chr = chr, 
                strand = strand)
            mcols(gr) <- NULL

            ## ------------ Setting up parallel computation ------------ #
            no_cores <- detectCores() - 1  
            cl <- makeCluster(no_cores, type="FORK")  
            registerDoParallel(cl)  
            
            if (is.null(nms)) 
                nms <- paste0("aln_", seq_len(nr))
            nams <- function(x) {
                nms <- rev(nms[ (length(nms) - seq_along(x) + 1) ] )
                nms
            }
            
            ## -------------------------------------------------------- #
            
            seqs <- try(foreach(k = seq_len(nr - 1), 
                        .final = function(x) setNames(x, nms[ -nr ])) %:%
                        foreach(j = seq((k + 1), nr, 1), 
                            .final = function(x) setNames(x, nams(x))) %dopar%
                {
                    aln <- seqs@unmasked[ c(k, j) ]
                    
                    aln <- selectAutomorphism(
                                            seq = aln, 
                                            filepath = NULL, 
                                            group = group, 
                                            cube = cube, 
                                            cube_alt = cube_alt, 
                                            start = start, 
                                            end = end, 
                                            chr = chr, 
                                            strand = strand)
                    mcols(aln)
                }, silent = TRUE)
            
            if (inherits(seqs, "try-error"))  {   
                stopCluster(cl)
                stop("*** Automorphism cannot be computed from
                     the MSA.")
            }
            else {
                stopCluster(cl)
                seqs <- unlist(seqs, recursive = FALSE)
                seqs <- as.AutomorphismList(seqs, grs = gr)
            }
        }

        return(seqs)
    }
)


## ========================== automorphismByRanges ============================

#' @aliases automorphismByRanges
#' @rdname automorphism
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @export
setGeneric("automorphismByRanges",
    function(
            autm,
            ...)
    standardGeneric("automorphismByRanges"))


#' @aliases automorphismByRanges
#' @rdname automorphism
#' @param autm An Automorphism-class object returned by function 
#' \code{\link{automorphism}}.
#' @importFrom data.table data.table
#' @export
setMethod("automorphismByRanges", signature(autm = "Automorphism"),
    function(autm) {
        i <- 1
        l <- length(autm)
        
        if (!inherits(autm, "GRanges")) {
            gr <- try(autm@SeqRanges, silent = TRUE)
            if (!inherits(gr, "try-error")) { 
                autm <- getAutomorphisms(autm)
                autm <- autm@DataList[[1]]
            }
        } 

        idx <- vector(mode = "numeric", length = length(autm))
        cube <- autm$cube[1]
        for (k in seq_len(l)) {
            if ( autm$cube[k] != cube ) {
                i <- i + 1
                cube <- autm$cube[k]
            } 
            idx[ k ] <- i
        }
        
        autm$idx <- factor(idx)
        autm <- data.table(data.frame(autm)) 
        autm <- autm[, list(
            seqnames = unique(seqnames), start = min(start),
            end = max(end), strand = unique(strand), 
            cube = unique(cube)), 
            by = idx ]
        autm <- autm[, c("seqnames", "start", "end", "strand", "cube")]
        autm <- makeGRangesFromDataFrame(autm, keep.extra.columns = TRUE)
        return(autm)
    }
)


#' @aliases automorphismByRanges
#' @rdname automorphism
#' @param autm An AutomorphismList-class object returned by function 
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
setMethod("automorphismByRanges", signature(autm = "AutomorphismList"),
    function(
            autm,
            min.len = 1L,
            num.cores = detectCores(),
            tasks = 0L,
            verbose = TRUE) {
        
        gr <- try(autm@SeqRanges, silent = TRUE)
        if (!inherits(gr, "try-error"))  
            autm <- getAutomorphisms(autm)
        
        autm <- autm@DataList
        
        ## ---------------- Setting parallel distribution --------------- ##
        
        if (num.cores > 1) 
            num.cores <- num.cores - 1
        progressbar = FALSE
        if (verbose) progressbar = TRUE
        if (Sys.info()["sysname"] == "Linux")
            bpparam <- MulticoreParam(workers = num.cores, tasks = tasks,
                                    progressbar = progressbar)
        else bpparam <- SnowParam(workers = num.cores, type = "SOCK",
                                progressbar = progressbar)
        
        ## -------------------------------------------------------------- ##
        
        if (length(gr) > 0) {
            autm <- lapply(autm, function(x) {
                        mcols(gr) <- x
                        x <- automorphismByRanges(x)
                        return(x)
                    })
        }
        
        idx <- which(sapply(autm, function(x) length(x) > min.len))
        autm <- autm[ idx ]
        
        return(as(autm, "GRangesList"))
    }
)


## ========================== automorphismByCoef ============================

#' @aliases automorphismByCoef
#' @rdname automorphism
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @export
setGeneric("automorphismByCoef",
           function(
               autm,
               ...)
               standardGeneric("automorphismByCoef"))

#' @aliases automorphismByCoef
#' @rdname automorphism
#' @param autm An automorphism-class object returned by function 
#' \code{\link{automorphism}}.
#' @importFrom data.table data.table
#' @export
setMethod("automorphismByCoef", signature(autm = "Automorphism"),
    function(autm) {
        i <- 1
        l <- length(autm)
        idx <- vector(mode = "numeric", length = length(autm))
        coefs <- autm$autm[1]
        for (k in seq_len(l)) {
            if ( autm$autm[k] != coefs ) {
                i <- i + 1
                coefs <- autm$autm[k]
            } 
            idx[ k ] <- i
        }
        
        autm$idx <- factor(idx)
        autm <- data.table(data.frame(autm))
        autm <- autm[, list(
            seqnames = unique(seqnames), start = min(start),
            end = max(end), strand = unique(strand), 
            autm = unique(autm),
            cube = unique(cube)), 
            by = idx ]
        autm <- autm[, c( "seqnames", "start", "end", 
                          "strand", "autm", "cube") ]
        autm <- makeGRangesFromDataFrame(autm, keep.extra.columns = TRUE)
        return(as(autm, "AtomorphismByCoef"))
    }
)



#' @aliases automorphismByCoef
#' @rdname automorphism
#' @param autm An AutomorphismList-class object returned by function 
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
setMethod("automorphismByCoef", signature(autm = "AutomorphismList"),
    function(
            autm,
            min.len = 1L,
            num.cores = detectCores(),
            tasks = 0L,
            verbose = TRUE) {
              
        gr <- try(autm@SeqRanges, silent = TRUE)
        if (!inherits(gr, "try-error"))  
            autm <- getAutomorphisms(autm)
        
        autm <- autm@DataList
        
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
            autm <- lapply(autm, function(x) {
                mcols(gr) <- x
                x <- automorphismByCoef(x)
                return(x)
            })
        }
        
        idx <- which(sapply(autm, function(x) length(x) > min.len))
        autm <- autm[ idx ]
        
        # return(as(autm, "AtomorphismByCoefList"))
        return(autm)
    }
)


## ========================== getAutomorphisms ============================

#' @rdname automorphism
#' @aliases getAutomorphisms
#' @importFrom S4Vectors mcols
#' @export
setGeneric("getAutomorphisms",
    function(
        x,
        ...)
    standardGeneric("getAutomorphisms"))


#' @rdname automorphism
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


#' @rdname automorphism
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


#' @rdname automorphism
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

## ===================== Auxiliary functions ===========================

selectAutomorphism <- function(
    seq, 
    filepath, 
    group, 
    cube, 
    cube_alt, 
    start, 
    end, 
    chr, 
    strand) {
    
    seq <- switch(group,
                  "Z5" = autZ5(
                      seq = seq,
                      filepath = filepath,
                      cube = cube,
                      cube_alt = cube_alt,
                      start = start,
                      end = end,
                      chr = chr,
                      strand = strand),
                  "Z64" = autZ64(
                      seq = seq,
                      filepath = filepath,
                      cube = cube,
                      cube_alt = cube_alt,
                      start = start,
                      end = end,
                      chr = chr,
                      strand = strand),
                  "Z5^3" = aut3D(
                      seq = seq,
                      filepath = filepath,
                      cube = cube,
                      cube_alt = cube_alt,
                      start = start,
                      end = end,
                      chr = chr,
                      strand = strand),
                  "Z125" = autZ125(
                      seq = seq,
                      filepath = filepath,
                      cube = cube,
                      cube_alt = cube_alt,
                      start = start,
                      end = end,
                      chr = chr,
                      strand = strand)
    )
    return(seq)
}


















