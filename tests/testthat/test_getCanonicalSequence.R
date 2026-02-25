## Test getCanonicalSequence

test_that("getCanonicalSequence removes delta mass modifications", {
    expect_equal(
        getCanonicalSequence("EM[-15.9949]EVEES[+79.9663]PEK"),
        "EMEVEESPEK"
    )
    expect_equal(
        getCanonicalSequence("APEPT[+79.966]IDEK"),
        "APEPTIDEK"
    )
})

test_that("getCanonicalSequence removes name modifications", {
    expect_equal(
        getCanonicalSequence("EK[Acetyl]EVEES[Phospho]PEK"),
        "EKEVEESPEK"
    )
    expect_equal(
        getCanonicalSequence("M[Oxidation]PEPTIDE"),
        "MPEPTIDE"
    )
})

test_that("getCanonicalSequence removes unimod ID modifications", {
    expect_equal(
        getCanonicalSequence("EK[UNIMOD:1]EVEES[UNIMOD:21]PEK"),
        "EKEVEESPEK"
    )
})

test_that("getCanonicalSequence removes terminal modifications", {
    expect_equal(getCanonicalSequence("[+304]-ATCK"), "ATCK")
    expect_equal(getCanonicalSequence("ATCK-[+42]"), "ATCK")
    expect_equal(getCanonicalSequence("[+304]-ATCK-[+42]"), "ATCK")
})

test_that("getCanonicalSequence handles unmodified sequences", {
    expect_equal(getCanonicalSequence("PEPTIDE"), "PEPTIDE")
    expect_equal(getCanonicalSequence("ATK"), "ATK")
})

test_that("getCanonicalSequence handles vectorised input", {
    result <- getCanonicalSequence(c(
        "EM[-15.9949]EVEES[+79.9663]PEK",
        "EK[Acetyl]EVEES[UNIMOD:21]PEK",
        "PEPTIDE"
    ))
    expect_equal(result, c("EMEVEESPEK", "EKEVEESPEK", "PEPTIDE"))
})

test_that("getCanonicalSequence removes stacked modifications", {
    expect_equal(getCanonicalSequence("AT[+10][-5]K"), "ATK")
    expect_equal(getCanonicalSequence("A[Phospho][UNIMOD:1]TK"), "ATK")
})


# Test internal helpers

test_that(".extractTerminalMods parses N-terminal modifications", {
    result <- PTMods:::.extractTerminalMods("[+304]-ATK")
    expect_equal(result$sequence, "ATK")
    expect_equal(result$nterm, "[+304]-")
    expect_equal(result$cterm, "")
})

test_that(".extractTerminalMods parses C-terminal modifications", {
    result <- PTMods:::.extractTerminalMods("ATK-[+42]")
    expect_equal(result$sequence, "ATK")
    expect_equal(result$nterm, "")
    expect_equal(result$cterm, "-[+42]")
})

test_that(".extractTerminalMods parses both terminal modifications", {
    result <- PTMods:::.extractTerminalMods("[+304]-ATK-[+42]")
    expect_equal(result$sequence, "ATK")
    expect_equal(result$nterm, "[+304]-")
    expect_equal(result$cterm, "-[+42]")
})

test_that(".extractTerminalMods handles no terminal modifications", {
    result <- PTMods:::.extractTerminalMods("ATK")
    expect_equal(result$sequence, "ATK")
    expect_equal(result$nterm, "")
    expect_equal(result$cterm, "")
})

test_that(".extractTerminalMods preserves internal modifications", {
    result <- PTMods:::.extractTerminalMods("[+304]-AT[+79.966]K-[+42]")
    expect_equal(result$sequence, "AT[+79.966]K")
    expect_equal(result$nterm, "[+304]-")
    expect_equal(result$cterm, "-[+42]")
})

test_that(".reattachTerminalMods prepends and appends correctly", {
    expect_equal(
        PTMods:::.reattachTerminalMods("ATK", "[+304]-", ""),
        "[+304]-ATK"
    )
    expect_equal(
        PTMods:::.reattachTerminalMods("ATK", "", "-[+42]"),
        "ATK-[+42]"
    )
    expect_equal(
        PTMods:::.reattachTerminalMods("ATK", "[+304]-", "-[+42]"),
        "[+304]-ATK-[+42]"
    )
    expect_equal(
        PTMods:::.reattachTerminalMods("ATK", "", ""),
        "ATK"
    )
})
