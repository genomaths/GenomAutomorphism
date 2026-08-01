# GenomAutomorphism 1.15.0

* Fix a bug in the internal function 'weighted_manhattan' (used by
  'codon_dist', 'codon_dist_matrix', and 'aminoacid_dist') where, for
  \eqn{group = "Z4"}, the second coordinate weight ('w[2]') was silently
  omitted and 'w[3]' was applied twice instead, contradicting the documented
  weighted Manhattan distance formula. This changes the numeric distances
  returned for \emph{group = "Z4"} (the default). The bundled 'cdm_z64'
  dataset was regenerated accordingly.

* Documentation updated 

# GenomAutomorphism 1.8.1

* Add new a new function: 'automorphism_prob', which applies a 
  Dirichlet-Multinomial Modelling (in a Bayesian framework) to compute the 
  posterior probability of each type of mutational event.

# GenomAutomorphism 1.5.1

* Introducing new functions for DNA and aminoacid sequence representations
  with physicochemical properties of DNA and aminoacids, which would be 
  useful for further downstream statistical analysis in R.

# GenomAutomorphism 1.0.2

* Fix error of parallel computation on Windows (12/08/2022)

# GenomAutomorphism 1.0.1

* Expanding analyses by including aminoacid similarity based on codon
  distances. Three new functions are added: codon_dist, codon_dist_matrix,
  and aminoacid_dist. See a tutorial applying these functions at 
  https://is.gd/oYLDK4.

# GenomAutomorphism 1.0.0

* Release in Bioconductor (version 3.16). 
  https://doi.org/doi:10.18129/B9.bioc.GenomAutomorphism

# GenomAutomorphism 0.99.4

* Expanding analyses by including amino acid similarity and statistical 
  protein contact potentials matrices from from Amino Acid Index Database 
  https://www.genome.jp/aaindex/.
* Improving documentation.
* Fixed a bug which introduced a change of protein coding frameshift.

# GenomAutomorphism 0.99.3

* Documentation improvement.      
* Package accepted on Bioconductor (07/18/22)

# GenomAutomorphism 0.99.2

* Updating several details after the review process in 
  Bioconductor (https://github.com/Bioconductor/Contributions/issues/2678)

# GenomAutomorphism 0.99.0

* Available at https://github.com/genomaths/GenomAutomorphism
* Release:
    February 28, 2022
    * Initial development.
