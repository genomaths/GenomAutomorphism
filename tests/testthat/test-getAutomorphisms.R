test_that("GenomAutomorphism function test", {
    data(autm, package = "GenomAutomorphism")
    aut <- autm@elementMetadata
    test1 <- length(ranges(getAutomorphisms(aut))) == 7100
    test2 <- inherits(getAutomorphisms(list(aut, aut)), "AutomorphismList")
    expect_true(all(test1,test2))
})
