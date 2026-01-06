test_that(".aamass", {
    expect_error(PTMods:::.aamass(1:3))
    expect_error(PTMods:::.aamass("ACE", "FOO"))
    expect_equal(PTMods:::.aamass("ACE"), PTMods:::.aamass("ACE", "MonoMass"))
    expect_equal(PTMods:::.aamass("ACE"), 303.088892)
    expect_equal(PTMods:::.aamass("ACE", "AvgMass"), 303.3348)
    expect_equal(PTMods:::.aamass(c(foo="ACE", bar="CEA")),
                 c(foo=303.088892, bar=303.088892))
})

test_that(".unimodMass", {
    expect_error(PTMods:::.unimodMass(1:3))
    expect_error(PTMods:::.unimodMass("ACE"))
    expect_error(PTMods:::.unimodMass("ACE", "FOO"))
    expect_error(PTMods:::.unimodMass("ACE", c("FOO", "BAR")))
    expect_equal(PTMods:::.unimodMass("ACE", "Met-loss:P-M"), 0)
    expect_equal(PTMods:::.unimodMass("MACE", "Met-loss:P-M"),
                 PTMods:::.unimodMass("MACE", "Met-loss:P-M", "MonoMass"))
    expect_equal(PTMods:::.unimodMass("MACE", "Met-loss:P-M", "MonoMass"),
                 -131.040485)
    expect_equal(PTMods:::.unimodMass("MACE", "Met-loss:P-M", "AvgMass"),
                 -131.1961)
    expect_equal(PTMods:::.unimodMass("MCCE", "Met-loss:P-M"),
                 PTMods:::.unimodMass("MCCE", "Met-loss+Acetyl:P-M"))
    expect_equal(PTMods:::.unimodMass("MACE", "Met-loss+Acetyl:P-M"),
                 -89.02992)
    expect_equal(PTMods:::.unimodMass("ACE", "Acetyl:N-term"), 42.010565)
    expect_equal(PTMods:::.unimodMass("ACE", "Acetyl:C"), 42.010565)
    expect_equal(PTMods:::.unimodMass("ACE", "Acetyl:H"), 0)
    expect_equal(PTMods:::.unimodMass(c("ABE", "ACE"), "Acetyl:C"),
                 c(0, 42.010565))
    expect_message(PTMods:::.unimodMass("ACE", "Unknown:420:N-term"),
                   "Applying the default rule for .* create an issue")
    expect_silent(PTMods:::.unimodMass("ACE", "Unknown:420:N-term", msg=FALSE))
})

test_that(".unimodSequence", {
    expect_error(PTMods:::.unimodSequence(1:3))
    expect_error(PTMods:::.unimodSequence("ACE"))
    expect_error(PTMods:::.unimodSequence("ACE", "FOO"))
    expect_error(PTMods:::.unimodSequence("ACE", c("FOO", "BAR")))
    expect_equal(PTMods:::.unimodSequence("ACE", "Met-loss:P-M"), "ACE")
    expect_equal(PTMods:::.unimodSequence("MACE", "Met-loss:P-M"), "ACE")
    expect_equal(PTMods:::.unimodSequence("MACE", "Met-loss+Acetyl:P-M"),
                 "ACE")
    expect_equal(PTMods:::.unimodSequence(c("MBCE", "MACE"), "Met-loss:P-M"),
                 c("MBCE", "ACE"))
})

test_that(".mass", {
    expect_error(PTMods:::.mass(1:3))
    expect_error(PTMods:::.mass("ACE", "FOO"))
    expect_error(PTMods:::.mass("ACE", fixedModifications="FOO"), "is not part")
    expect_error(PTMods:::.mass("ACE",
                                fixedModifications=c("FOO", "Met-loss:P-M",
                                                     "BAR")), "are not part")
    expect_error(PTMods:::.mass("ACE", fixedModifications=1:3),
                 "must be a `character`")
    expect_error(PTMods:::.mass("ACE",
                                fixedModifications=c("Acetyl:K",
                                                     "Carbamidomethyl:K"),
                       "Duplicated fixed modification sites are not allowed"))
    r <- 303.088892
    attr(r, "sequence") <- "ACE"
    expect_equal(PTMods:::.mass("ACE"), r)
    expect_equal(PTMods:::.mass("ACE", fixedModifications="Acetyl:K"), r)
    r <- 345.099457
    attr(r, "sequence") <- "ACE"
    expect_equal(PTMods:::.mass("ACE", fixedModifications="Acetyl:N-term"), r)
    expect_equal(PTMods:::.mass("MACE",
                                fixedModifications="Met-loss+Acetyl:P-M"), r)
    r <- 402.120921
    attr(r, "sequence") <- "ACE"
    expect_equal(PTMods:::.mass("ACE",
                                fixedModifications=c("Acetyl:N-term",
                                                     "Carbamidomethyl:C")), r)
    r <- 478.119206
    attr(r, "sequence") <- "MDCE"
    expect_equal(PTMods:::.mass("MDCE", fixedModifications="Met-loss+Acetyl:P-M"), r)
})

test_that(".countSite", {
    expect_error(PTMods:::.countSite(1:3, "C"))
    expect_error(PTMods:::.countSite("ACE", 1:3))
    expect_error(PTMods:::.countSite("ACE", c("A", "C")))
    expect_equal(PTMods:::.countSite("ACE", "Nterm"), 1)
    expect_equal(PTMods:::.countSite("ACE", "Cterm"), 1)
    expect_equal(PTMods:::.countSite("ACE", "C"), 1)
    expect_equal(PTMods:::.countSite(c("ACCE", "ACE"), "Nterm"), c(1, 1))
    expect_equal(PTMods:::.countSite(c("ACCE", "ACE"), "Cterm"), c(1, 1))
    expect_equal(PTMods:::.countSite(c("ACCE", "ACE"), "C"), 2:1)
    expect_equal(PTMods:::.countSite(c("ACCE", "ACE"), "H"), c(0, 0))
})
