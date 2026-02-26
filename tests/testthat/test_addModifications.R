# Test addModifications

test_that("addModifications with fixed modifications only", {
    result <- addModifications(
        c("ATCK", "PEPCK"),
        fixedModifications = c(C = 57.02146)
    )
    expect_equal(result, c("ATC[+57.02146]K", "PEPC[+57.02146]K"))
})

test_that("addModifications with variable modifications only", {
    result <- addModifications(
        "ATSK",
        variableModifications = c(S = 79.966331, T = 79.966331),
        maxMods = 1
    )
    ## Modifiable: T, S → 1+2 = 3 combos
    expect_length(result, 3)
    expect_true("ATSK" %in% result)
    expect_true("AT[+79.966331]SK" %in% result)
    expect_true("ATS[+79.966331]K" %in% result)
})

test_that("addModifications with both fixed and variable modifications", {
    result <- addModifications(
        "ATCSK",
        fixedModifications = c(C = 57.02146),
        variableModifications = c(S = 79.966331, T = 79.966331),
        maxMods = 1
    )
    ## Fixed first: ATC[+57.02146]SK
    ## Then variable on T and S: 1+2 = 3 combos
    expect_length(result, 3)
    expect_true("ATC[+57.02146]SK" %in% result)
    expect_true("AT[+79.966331]C[+57.02146]SK" %in% result)
    expect_true("ATC[+57.02146]S[+79.966331]K" %in% result)
})

test_that("addModifications returns annotations in name style", {
    result <- addModifications(
        "ATSK",
        variableModifications = c(S = 79.966331, T = 79.966331),
        convertToStyle = "name",
        maxMods = 1
    )
    expect_true("ATSK" %in% result)
    expect_true("AT[Phospho]SK" %in% result)
    expect_true("ATS[Phospho]K" %in% result)
})

test_that("addModifications returns annotations in unimodId style", {
    result <- addModifications(
        "ATSK",
        variableModifications = c(S = 79.966331, T = 79.966331),
        convertToStyle = "unimodId",
        maxMods = 1
    )
    expect_true("ATSK" %in% result)
    expect_true("AT[UNIMOD:21]SK" %in% result)
    expect_true("ATS[UNIMOD:21]K" %in% result)
})

test_that("addModifications with N-terminal fixed modification", {
    result <- addModifications(
        "ATSK",
        fixedModifications = c(Nterm = 304),
        variableModifications = c(S = 79.966331)
    )
    ## Fixed: [+304]-ATSK
    ## Variable on S: 2 combos
    expect_length(result, 2)
    expect_true("[+304]-ATSK" %in% result)
    expect_true("[+304]-ATS[+79.966331]K" %in% result)
})

test_that("addModifications with NULL fixed and variable returns unchanged", {
    result <- addModifications("PEPTIDE",
        fixedModifications = NULL,
        variableModifications = NULL
    )
    expect_equal(result, "PEPTIDE")
})

test_that("addModifications passes maxMods through dots", {
    result <- addModifications(
        "ATSK",
        variableModifications = c(A = 4, T = 5, S = 8),
        maxMods = 2
    )
    ## Modifiable:
    ## 0: 1, 1-mod: 3, 2-mod: 3 → total 7
    expect_length(result, 7)
})
