context("convertAnnotation")

# Test main convertAnnotation function ====

test_that("convertAnnotation converts name to delta_mass", {
    # Single modification
    expect_equal(convertAnnotation("M[Oxidation]PEPTIDE", convertToStyle = "deltaMass"),
                 "M[+15.994915]PEPTIDE")

    # Multiple modifications
    expect_equal(convertAnnotation("M[Oxidation]EVNES[Phospho]PEK",
                                   convertToStyle = "deltaMass"),
                 "M[+15.994915]EVNES[+79.966331]PEK")

    # No modifications
    expect_equal(convertAnnotation("PEPTIDE", convertToStyle = "deltaMass"),
                 "PEPTIDE")
})

test_that("convertAnnotation converts name to unimod_id", {
    # Single modification
    expect_equal(convertAnnotation("M[Oxidation]PEPTIDE", convertToStyle = "unimodId"),
                 "M[UNIMOD:35]PEPTIDE")

    # Multiple modifications
    expect_equal(convertAnnotation("C[Carbamidomethyl]PEPTIDE",
                                   convertToStyle = "unimodId"),
                 "C[UNIMOD:4]PEPTIDE")
})

test_that("convertAnnotation converts delta_mass to name", {
    # Single modification
    expect_equal(convertAnnotation("M[+15.995]PEPTIDE", convertToStyle = "name"),
                 "M[Oxidation]PEPTIDE")

    # Multiple modifications
    expect_equal(convertAnnotation("S[+79.966]EQUENCE", convertToStyle = "name"),
                 "S[Phospho]EQUENCE")
})

test_that("convertAnnotation converts delta_mass to unimod_id", {
    expect_equal(convertAnnotation("M[+15.995]PEPTIDE", convertToStyle = "unimodId"),
                 "M[UNIMOD:35]PEPTIDE")

    expect_equal(convertAnnotation("S[+79.966]PEPTIDE", convertToStyle = "unimodId"),
                 "S[UNIMOD:21]PEPTIDE")
})

test_that("convertAnnotation converts unimod_id to name", {
    expect_equal(convertAnnotation("M[UNIMOD:35]PEPTIDE", convertToStyle = "name"),
                 "M[Oxidation]PEPTIDE")

    expect_equal(convertAnnotation("C[UNIMOD:4]PEPTIDE", convertToStyle = "name"),
                 "C[Carbamidomethyl]PEPTIDE")
})

test_that("convertAnnotation converts unimod_id to delta_mass", {
    expect_equal(convertAnnotation("M[UNIMOD:35]PEPTIDE",
                                   convertToStyle = "deltaMass"),
                 "M[+15.994915]PEPTIDE")

    expect_equal(convertAnnotation("S[UNIMOD:21]PEPTIDE",
                                   convertToStyle = "deltaMass"),
                 "S[+79.966331]PEPTIDE")
})

test_that("convertAnnotation handles vectorized input", {
    sequences <- c("M[Oxidation]PEPTIDE", "EVNES[Phospho]PEK", "PEPTIDE")
    result <- convertAnnotation(sequences, convertToStyle = "deltaMass")

    expect_length(result, 3)
    expect_equal(result[1], "M[+15.994915]PEPTIDE")
    expect_equal(result[2], "EVNES[+79.966331]PEK")
    expect_equal(result[3], "PEPTIDE")

    # No names should be preserved
    expect_null(names(result))
})

test_that("convertAnnotation respects mass_tolerance", {
    # With default tolerance (0.01), should match
    expect_equal(convertAnnotation("M[+15.995]PEPTIDE", convertToStyle = "name"),
                 "M[Oxidation]PEPTIDE")

    # With very strict tolerance, may not match
    expect_warning(
        result <- convertAnnotation("M[+16.1]PEPTIDE",
                                   convertToStyle = "name",
                                   massTolerance = 0.01),
        "Could not find Unimod entry"
    )
    expect_equal(result, "M[+16.1]PEPTIDE")
})

test_that("convertAnnotation handles sequences already in target format", {
    # Already in delta_mass format, converting to delta_mass
    expect_equal(convertAnnotation("M[+15.994915]PEPTIDE",
                                   convertToStyle = "deltaMass"),
                 "M[+15.994915]PEPTIDE")

    # Already in name format, converting to name
    expect_equal(convertAnnotation("M[Oxidation]PEPTIDE", convertToStyle = "name"),
                 "M[Oxidation]PEPTIDE")

    # Already in unimod_id format, converting to unimod_id
    expect_equal(convertAnnotation("M[UNIMOD:35]PEPTIDE",
                                   convertToStyle = "unimodId"),
                 "M[UNIMOD:35]PEPTIDE")
})

test_that("convertAnnotation input validation", {
    # Non-character input
    expect_error(convertAnnotation(123, convertToStyle = "deltaMass"),
                 "x must be a character vector")

    expect_error(convertAnnotation(list("M[Oxidation]PEPTIDE"),
                                   convertToStyle = "deltaMass"),
                 "x must be a character vector")

    # Invalid convert_to argument
    expect_error(convertAnnotation("M[Oxidation]PEPTIDE", convertToStyle = "invalid"),
                 "'arg' should be one of")
})

test_that("convertAnnotation warns on unknown modifications", {
    expect_warning(
        result <- convertAnnotation("M[UnknownMod]PEPTIDE",
                                   convertToStyle = "deltaMass"),
        "Could not find Unimod entry"
    )
    # Should return original when conversion fails
    expect_equal(result, "M[UnknownMod]PEPTIDE")
})


# Test internal helper functions ====

test_that(".detectModificationType correctly identifies modification types", {
    # Unimod ID formats
    expect_equal(unimod:::.detectModificationType("UNIMOD:35"), "unimodId")
    expect_equal(unimod:::.detectModificationType("U:35"), "unimodId")
    expect_equal(unimod:::.detectModificationType("unimod:35"), "unimodId")

    # Delta mass formats
    expect_equal(unimod:::.detectModificationType("+15.995"), "deltaMass")
    expect_equal(unimod:::.detectModificationType("-17.026"), "deltaMass")
    expect_equal(unimod:::.detectModificationType("+79.966331"), "deltaMass")
    expect_equal(unimod:::.detectModificationType("+15"), "deltaMass")

    # Name format
    expect_equal(unimod:::.detectModificationType("Oxidation"), "name")
    expect_equal(unimod:::.detectModificationType("Phospho"), "name")
    expect_equal(unimod:::.detectModificationType("Carbamidomethyl"), "name")
})

test_that(".extractUnimodId extracts numeric ID", {
    expect_equal(unimod:::.extractUnimodId("UNIMOD:35"), 35)
    expect_equal(unimod:::.extractUnimodId("U:35"), 35)
    expect_equal(unimod:::.extractUnimodId("unimod:35"), 35)
    expect_equal(unimod:::.extractUnimodId("UNIMOD:4"), 4)
})

test_that(".formatMass formats mass values correctly", {
    # Positive mass
    expect_equal(unimod:::.formatMass(15.994915), "+15.994915")
    expect_equal(unimod:::.formatMass(79.966331), "+79.966331")

    # Negative mass
    expect_equal(unimod:::.formatMass(-17.026549), "-17.026549")

    # Zero
    expect_equal(unimod:::.formatMass(0), "+0")

    # Custom precision
    expect_equal(unimod:::.formatMass(15.994915, digits = 3), "+15.995")
    expect_equal(unimod:::.formatMass(15.994915, digits = 2), "+15.99")
})

test_that(".lookupByName finds modifications correctly", {
    unimod_data <- unimod::modifications[!unimod::modifications$NeutralLoss,
                                         c("UnimodId", "Name", "MonoMass")]
    unimod_data <- unimod_data[!duplicated(unimod_data$Name), ]

    # Valid lookup
    result <- unimod:::.lookupByName("Oxidation", unimod_data)
    expect_type(result, "list")
    expect_equal(result$unimodId, 35)
    expect_equal(result$name, "Oxidation")
    expect_equal(result$mass, 15.994915, tolerance = 0.000001)

    # Invalid lookup
    result <- unimod:::.lookupByName("NonExistentMod", unimod_data)
    expect_null(result)
})

test_that(".lookupByMass finds modifications correctly", {
    unimod_data <- unimod::modifications[!unimod::modifications$NeutralLoss,
                                         c("UnimodId", "Name", "MonoMass")]
    unimod_data <- unimod_data[!duplicated(unimod_data$Name), ]

    # Valid lookup with tolerance
    result <- unimod:::.lookupByMass(15.995, unimod_data, 0.01)
    expect_type(result, "list")
    expect_equal(result$unimodId, 35)
    expect_equal(result$name, "Oxidation")

    # Exact match
    result <- unimod:::.lookupByMass(15.994915, unimod_data, 0.01)
    expect_type(result, "list")
    expect_equal(result$unimodId, 35)

    # No match with strict tolerance
    result <- unimod:::.lookupByMass(999.999, unimod_data, 0.01)
    expect_null(result)
})

test_that(".lookupByUnimodId finds modifications correctly", {
    unimod_data <- unimod::modifications[!unimod::modifications$NeutralLoss,
                                         c("UnimodId", "Name", "MonoMass")]
    unimod_data <- unimod_data[!duplicated(unimod_data$Name), ]

    # Valid lookup
    result <- unimod:::.lookupByUnimodId(35, unimod_data)
    expect_type(result, "list")
    expect_equal(result$unimodId, 35)
    expect_equal(result$name, "Oxidation")
    expect_equal(result$mass, 15.994915, tolerance = 0.000001)

    # Another valid lookup
    result <- unimod:::.lookupByUnimodId(21, unimod_data)
    expect_type(result, "list")
    expect_equal(result$unimodId, 21)
    expect_equal(result$name, "Phospho")

    # Invalid lookup
    result <- unimod:::.lookupByUnimodId(99999, unimod_data)
    expect_null(result)
})

test_that(".convertAnnotation handles single sequences", {
    # Valid conversion
    expect_equal(unimod:::.convertAnnotation("M[Oxidation]PEPTIDE",
                                            convertToStyle = "deltaMass"),
                 "M[+15.994915]PEPTIDE")

    # Input validation
    expect_error(unimod:::.convertAnnotation(c("A", "B"),
                                            convertToStyle = "deltaMass"),
                 "must be a single character string")

    expect_error(unimod:::.convertAnnotation(123, convertToStyle = "deltaMass"),
                 "must be a single character string")
})
