# GenomAutomorphism

Robersy Sanchez  
Department of Biology. Eberly College of Science.  
Pennsylvania State University, University Park, PA 16802  
<rus547@psu.edu>  
[ORCID:
orcid.org/0000-0002-5246-1453](https://orcid.org/0000-0002-5246-1453)

## Overview

This is a R package to compute the autimorphisms between pairwise
aligned DNA sequences represented as elements from a Genomic Abelian
group as described in the paper [Genomic Abelian Finite
Groups](https://www.biorxiv.org/content/10.1101/2021.06.01.446543v2). In
a general scenario, whole chromosomes or genomic regions from a
population (from any species or close related species) can be
algebraically represented as a direct sum of cyclic groups or more
specifically Abelian *p*-groups. Basically, we propose the
representation of multiple sequence alignments (MSA) of length *N* as a
finite Abelian group created by the direct sum of Abelian group of
*prime-power order*:

$$
\\begin{aligned}
G = (\\mathbb{Z}\_{p^{\\alpha\_{1}}\_1})^{n\_1} \\oplus (\\mathbb{Z}\_{p^{\\alpha\_{2}}\_1})^{n\_2} \\oplus \\dots \\oplus (\\mathbb{Z}\_{p^{\\alpha\_{k}}\_k})^{n\_k} 
\\end{aligned}
$$

Where, the *p*<sub>*i*</sub>’s are prime numbers, *α*<sub>*i*</sub> ∈ ℕ
and ℤ<sub>*p*<sub>*i*</sub><sup>*α*<sub>*i*</sub></sup></sub> is the
group of integer modulo *p*<sub>*i*</sub><sup>*α*<sub>*i*</sub></sup>.

For the purpose of automorphism between two aligned DNA sequences,
*p*<sub>*i*</sub><sup>*α*<sub>*i*</sub></sup> ∈ {5, 2<sup>6</sup>, 5<sup>3</sup>}.

------------------------------------------------------------------------

## Status

This application is under development. Watch this repo or check for
updates.

------------------------------------------------------------------------

## Dependences

This package depends, so far, from: *Biostrings*, *GenomicRanges*,
*numbers*, and *S4Vectors*.

------------------------------------------------------------------------

## Installation of R dependencies:

        if (!requireNamespace("BiocManager")) install.packages("BiocManager")
        BiocManager::install()
        
        BiocManager::install(c("Biostrings", "GenomicRanges", "S4Vectors"))
        install.packages(c("numbers", "devtools"), dependencies=TRUE)

------------------------------------------------------------------------

## You can install **GenomAutomorphism** package from GitHub

       devtools::install_git("https://github.com/genomaths/GenomAutomorphism.git")

------------------------------------------------------------------------

# References

1.  Sanchez R, Morgado E, Grau R. Gene algebra from a genetic code
    algebraic structure. J Math Biol. 2005 Oct;51(4):431-57. doi:
    10.1007/s00285-005-0332-8. Epub 2005 Jul 13. PMID: 16012800. (
    [PDF](https://arxiv.org/pdf/q-bio/0412033.pdf)).

2.  Robersy Sanchez, Jesús Barreto (2021) Genomic Abelian Finite Groups.
    [doi:
    10.1101/2021.06.01.446543](https://doi.org/10.1101/2021.06.01.446543).

3.  M. V José, E.R. Morgado, R. Sánchez, T. Govezensky, The 24 possible
    algebraic representations of the standard genetic code in six or in
    three dimensions, Adv. Stud. Biol. 4 (2012)
    119–152.[PDF](https://is.gd/na9eap).

4.  R. Sanchez. Symmetric Group of the Genetic–Code Cubes. Effect of the
    Genetic–Code Architecture on the Evolutionary Process MATCH Commun.
    Math. Comput. Chem. 79 (2018) 527-560.
    [PDF](https://bit.ly/2Z9mjM7).
