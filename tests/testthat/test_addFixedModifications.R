# Test addFixedModifications

test_that("addFixedModifications applies numeric fixed modification", {
    result <- addFixedModifications("ATCK", fixedModifications = c(C = 57.02146))
    expect_equal(result, "ATC[+57.02146]K")
})

test_that("addFixedModifications applies multiple numeric modifications", {
    result <- addFixedModifications("ATSK",
        fixedModifications = c(
            A = -42,
            S = 79.966331
        )
    )
    expect_equal(result, "A[-42]TS[+79.966331]K")
})

test_that("addFixedModifications applies character name modifications preserving style", {
    result <- addFixedModifications("ATK",
        fixedModifications = c(
            T = "Phospho",
            A = "UNIMOD:1"
        )
    )
    ## Annotation style of the parameter is preserved in the output
    expect_equal(result, "A[UNIMOD:1]T[Phospho]K")
})

test_that("addFixedModifications preserves existing annotation style while adding character mods", {
    result <- addFixedModifications("A[+42]TK",
        fixedModifications = c(T = "Phospho")
    )
    ## Existing [+42] stays in deltaMass, new mod uses name style as supplied
    expect_equal(result, "A[+42]T[Phospho]K")
})

test_that("addFixedModifications preserves existing modifications", {
    result <- addFixedModifications("AT[+79.966]AK",
        fixedModifications = c(A = 42)
    )
    expect_equal(result, "A[+42]T[+79.966]A[+42]K")
})

test_that("addFixedModifications stacks on already-modified residues", {
    result <- addFixedModifications("ATA[+79.966]K",
        fixedModifications = c(A = 42)
    )
    expect_equal(result, "A[+42]TA[+79.966][+42]K")
})

test_that("addFixedModifications handles vectorised input", {
    result <- addFixedModifications(c("ATCK", "PEPCK"),
        fixedModifications = c(C = 57.02146)
    )
    expect_length(result, 2)
    expect_equal(result, c("ATC[+57.02146]K", "PEPC[+57.02146]K"))
    expect_null(names(result))
})

test_that("addFixedModifications handles sequences without target amino acids", {
    result <- addFixedModifications("ATGK", fixedModifications = c(C = 57.02146))
    expect_equal(result, "ATGK")
})


test_that("addFixedModifications adds + prefix for numeric-as-character values", {
    result <- addFixedModifications("ATCK",
        fixedModifications = c(Nterm = "Phospho", C = 57.02146)
    )
    expect_equal(result, "[Phospho]-ATC[+57.02146]K")
})

test_that("addFixedModifications preserves - prefix for negative numeric-as-character values", {
    result <- addFixedModifications("ATK",
        fixedModifications = c(A = "-42", T = "79.966331")
    )
    expect_equal(result, "A[-42]T[+79.966331]K")
})

test_that("addFixedModifications applies two mods to same residue consecutively", {
    result <- "ATK" |>
        addFixedModifications(c(T = "+57.024")) |>
        addFixedModifications(c(T = "Phospho"))
    expect_equal(result, "AT[+57.024][Phospho]K")
})

# Test internal .addFixedModifications

test_that(".addFixedModifications applies N-terminal modification", {
    result <- PTMods:::.addFixedModifications("ATCK",
        fixedModifications = c(Nterm = 304, C = 57.02146)
    )
    expect_equal(result, "[+304]-ATC[+57.02146]K")
})

test_that(".addFixedModifications applies C-terminal modification", {
    result <- PTMods:::.addFixedModifications("ATCK",
        fixedModifications = c(Cterm = 42)
    )
    expect_equal(result, "ATCK-[+42]")
})

test_that(".addFixedModifications applies both terminal modifications", {
    result <- PTMods:::.addFixedModifications("ATK",
        fixedModifications = c(Nterm = 304, Cterm = 42)
    )
    expect_equal(result, "[+304]-ATK-[+42]")
})

test_that(".addFixedModifications preserves existing terminal modifications", {
    result <- PTMods:::.addFixedModifications("[+10]-ATK",
        fixedModifications = c(A = 42)
    )
    expect_equal(result, "[+10]-A[+42]TK")
})

test_that(".addFixedModifications handles negative modifications", {
    result <- PTMods:::.addFixedModifications("ATK",
        fixedModifications = c(A = -42)
    )
    expect_equal(result, "A[-42]TK")
})


# Test internal .addTerminalModifications

test_that(".addTerminalModifications adds N-terminal modification", {
    result <- PTMods:::.addTerminalModifications("ATK", c(Nterm = 67))
    expect_equal(result, "[+67]-ATK")
})

test_that(".addTerminalModifications adds C-terminal modification", {
    result <- PTMods:::.addTerminalModifications("ATK", c(Cterm = 42))
    expect_equal(result, "ATK-[+42]")
})

test_that(".addTerminalModifications adds both terminal modifications", {
    result <- PTMods:::.addTerminalModifications(
        "ATK",
        c(Nterm = 304, Cterm = 42)
    )
    expect_equal(result, "[+304]-ATK-[+42]")
})

test_that(".addTerminalModifications skips if no terminal mods", {
    result <- PTMods:::.addTerminalModifications("ATK", c(A = 42))
    expect_equal(result, "ATK")
})
