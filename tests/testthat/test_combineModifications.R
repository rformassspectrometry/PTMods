# Test combineModifications

test_that("combineModifications sums stacked delta mass modifications", {
    expect_equal(combineModifications("AT[+10][-5]K"), "AT[+5]K")
})

test_that("combineModifications combines stacked N-terminal modifications", {
    expect_equal(combineModifications("[+10]-[-5]-ATK"), "[+5]-ATK")
})

test_that("combineModifications converts name and unimod to delta mass", {
    result <- combineModifications("AT[Phospho][UNIMOD:1]K")
    expect_equal(result, "AT[+121.976896]K")
})

test_that("combineModifications handles multiple modified residues", {
    expect_equal(
        combineModifications("EM[+15.9949][-10]EK[+42][-8]"),
        "EM[+5.9949]EK[+34]"
    )
})

test_that("combineModifications handles vectorised input", {
    result <- combineModifications(c("AT[+10][-5]K", "EM[+15.9949][-10]EK"))
    expect_equal(result, c("AT[+5]K", "EM[+5.9949]EK"))
})

test_that("combineModifications handles unmodified sequences", {
    expect_equal(combineModifications("PEPTIDE"), "PEPTIDE")
})

test_that("combineModifications with single modification is unchanged", {
    expect_equal(combineModifications("AT[+10]K"), "AT[+10]K")
})

test_that("combineModifications combines stacked C-terminal mods", {
    expect_equal(combineModifications("ATK-[+42]-[-8]"), "ATK-[+34]")
})

test_that("combineModifications combines both terminal stacked mods", {
    expect_equal(
        combineModifications("[+304]-[-54]-ATK-[+42]-[-8]"),
        "[+250]-ATK-[+34]"
    )
})


# Test internal .parseModifiedSequence

test_that(".parseModifiedSequence parses single modification", {
    result <- PTMods:::.parseModifiedSequence("APEPT[+79.966]IDEK")
    expect_equal(
        as.numeric(result),
        c(0, 0, 0, 0, 79.966, 0, 0, 0, 0)
    )
})

test_that(".parseModifiedSequence parses stacked modifications", {
    result <- PTMods:::.parseModifiedSequence("APEPT[+79.966][+10]IDEK")
    expect_equal(
        as.numeric(result),
        c(0, 0, 0, 0, 89.966, 0, 0, 0, 0)
    )
})

test_that(".parseModifiedSequence parses mixed sign stacked modifications", {
    result <- PTMods:::.parseModifiedSequence("EM[+15.9949][-10]EK")
    expect_equal(as.numeric(result), c(0, 5.9949, 0, 0))
})

test_that(".parseModifiedSequence sets terminal attributes", {
    result <- PTMods:::.parseModifiedSequence("[+304]-APEPTIDEK")
    expect_equal(attr(result, "Nterm"), 304)
    expect_true(is.na(attr(result, "Cterm")))

    result2 <- PTMods:::.parseModifiedSequence("ATK-[+42]")
    expect_true(is.na(attr(result2, "Nterm")))
    expect_equal(attr(result2, "Cterm"), 42)
})

test_that(".parseModifiedSequence handles unmodified sequence", {
    result <- PTMods:::.parseModifiedSequence("PEPTIDE")
    expect_equal(as.numeric(result), rep(0, 7))
})


# Test internal .combineTerminalMods

test_that(".combineTerminalMods combines stacked N-terminal mods", {
    expect_equal(
        PTMods:::.combineTerminalMods("[+304]-[-54]-ATK"),
        "[+250]-ATK"
    )
})

test_that(".combineTerminalMods combines stacked C-terminal mods", {
    expect_equal(
        PTMods:::.combineTerminalMods("ATK-[+42]-[-8]"),
        "ATK-[+34]"
    )
})

test_that(".combineTerminalMods combines both terminal stacks", {
    expect_equal(
        PTMods:::.combineTerminalMods("[+304]-[-54]-ATK-[+42]-[-8]"),
        "[+250]-ATK-[+34]"
    )
})

test_that(".combineTerminalMods handles no stacked mods", {
    expect_equal(PTMods:::.combineTerminalMods("ATK"), "ATK")
    expect_equal(PTMods:::.combineTerminalMods("[+304]-ATK"), "[+304]-ATK")
    expect_equal(PTMods:::.combineTerminalMods("ATK-[+42]"), "ATK-[+42]")
})
