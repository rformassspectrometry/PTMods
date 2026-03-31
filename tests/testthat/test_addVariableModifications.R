# Test addVariableModifications

test_that("addVariableModifications generates all combos with maxMods", {
    result <- addVariableModifications("APGHKA",
        variableModifications = c(A = 4, K = 5, S = 8),
        maxMods = 2
    )
    ## Modifiable positions:
    ## 0 mods: APGHKA
    ## 1 mod:  A[+4]PGHKA, APGHK[+5]A, APGHKA[+4]
    ## 2 mods: A[+4]PGHK[+5]A, A[+4]PGHKA[+4], APGHK[+5]A[+4]
    expect_length(result, 7)
    expect_true("APGHKA" %in% result)
    expect_true("A[+4]PGHKA" %in% result)
    expect_true("APGHK[+5]A" %in% result)
    expect_true("APGHKA[+4]" %in% result)
    expect_true("A[+4]PGHK[+5]A" %in% result)
    expect_true("A[+4]PGHKA[+4]" %in% result)
    expect_true("APGHK[+5]A[+4]" %in% result)
})

test_that("addVariableModifications generates all combos with maxMods = 3", {
    result <- addVariableModifications("APGHKA",
        variableModifications = c(A = 4, K = 5, S = 8),
        maxMods = 3
    )
    ## Same as maxMods=2 plus the 3-mod combo:
    ## A[+4]PGHK[+5]A[+4]
    expect_length(result, 8)
    expect_true("A[+4]PGHK[+5]A[+4]" %in% result)
})

test_that("addVariableModifications handles character modifications preserving style", {
    result <- addVariableModifications("ATK",
        variableModifications = c(T = "Phospho", A = "UNIMOD:1")
    )
    ## Annotation style of the parameter is preserved in the output:
    ## 0 mods: ATK
    ## 1 mod:  A[UNIMOD:1]TK, AT[Phospho]K
    ## 2 mods: A[UNIMOD:1]T[Phospho]K
    expect_length(result, 4)
    expect_true("ATK" %in% result)
    expect_true("A[UNIMOD:1]TK" %in% result)
    expect_true("AT[Phospho]K" %in% result)
    expect_true("A[UNIMOD:1]T[Phospho]K" %in% result)
})

test_that("addVariableModifications preserves existing annotation style while adding character mods", {
    result <- addVariableModifications("A[+42]TK",
        variableModifications = c(T = "Phospho")
    )
    ## Existing [+42] stays in deltaMass; new mod uses name style as supplied
    expect_true("A[+42]TK" %in% result)
    expect_true("A[+42]T[Phospho]K" %in% result)
    expect_length(result, 2)
})

test_that("addVariableModifications handles vectorised input", {
    result <- addVariableModifications(c("ATK", "PAK"),
        variableModifications = c(T = 10, A = 5)
    )
    ## ATK: modifiable A(1), T(2) → 4 combos
    ## PAK: modifiable A(2) → 2 combos
    expect_length(result, 6)
})

test_that("addVariableModifications preserves existing modifications", {
    result <- addVariableModifications("A[+10]KT[-5]",
        variableModifications = c(A = -5, T = 8)
    )
    ## Modifiable: position 1 (A), position 3 (T)
    expect_true("A[+10]KT[-5]" %in% result)
    expect_length(result, 4)
})

test_that("addVariableModifications respects maxMods = 1", {
    result <- addVariableModifications("ATSK",
        variableModifications = c(A = 4, T = 5, S = 8),
        maxMods = 1
    )
    expect_length(result, 4)
    expect_true("ATSK" %in% result)
    expect_true("A[+4]TSK" %in% result)
    expect_true("AT[+5]SK" %in% result)
    expect_true("ATS[+8]K" %in% result)
})

test_that("addVariableModifications with no matching residues returns original", {
    result <- addVariableModifications("GGG",
        variableModifications = c(A = 4, S = 8)
    )
    expect_equal(result, "GGG")
})

test_that("addVariableModifications preserves terminal modifications", {
    result <- addVariableModifications("[+304]-ATK",
        variableModifications = c(A = 5)
    )
    expect_true("[+304]-ATK" %in% result)
    expect_true("[+304]-A[+5]TK" %in% result)
    expect_length(result, 2)
})

test_that("addVariableModifications returns original if maxMods = 0", {
    result <- addVariableModifications("ATSK",
        variableModifications = c(A = 4, T = 5),
        maxMods = 0
    )
    expect_equal(result, "ATSK")
})

test_that("addVariableModifications applies two mods to same residue consecutively", {
    result <- "ATK" |>
        addVariableModifications(c(T = "+57.024")) |>
        addVariableModifications(c(T = "Phospho"))
    ## All four combinations: unmodified, first mod only, second mod only,
    ## and both mods stacked on T
    expect_length(result, 4)
    expect_true("ATK" %in% result)
    expect_true("AT[+57.024]K" %in% result)
    expect_true("AT[Phospho]K" %in% result)
    expect_true("AT[+57.024][Phospho]K" %in% result)
})

test_that("addVariableModifications returns original with NULL mods", {
    result <- addVariableModifications("ATSK",
        variableModifications = NULL
    )
    expect_equal(result, "ATSK")
})

test_that("addVariableModifications preserves NAs in input", {
    result <- addVariableModifications(c(NA, "ATK"))
    expect_true(is.na(result[[1]]))
    expect_equal(result[[2]], "ATK")
})
