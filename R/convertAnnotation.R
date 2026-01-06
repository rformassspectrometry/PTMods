utils::globalVariables("modifications")

#' @name convertAnnotation
#'
#' @aliases convertAnnotation .convertAnnotation
#'
#' @title Convert sequences from one annotation style to another
#'
#' @description
#' Converts modifications between different annotation formats for multiple
#' sequences at once. See the details and examples sections for more
#' information.  The annotation styles are inferred from the `modifications`
#' dataframe (see `?modifications`).
#'
#' @param x Character vector with peptide sequences in ProForma format
#'
#' @param convertToStyle Character string specifying target format. Options:
#'   "deltaMass", "unimodId", "name"
#'
#' @param massTolerance Numeric mass tolerance in Daltons for matching
#'   modifications (default: 0.01). Used when converting from deltaMass.
#'
#' @return Character vector with the sequences in the target annotation format
#'
#' @author Guillaume Deflandre <guillaume.deflandre@uclouvain.be>
#'
#' @details
#' The function handles three main conversion scenarios:
#'
#' - Name to deltaMass: "M[Oxidation]PEPTIDE" -> "M[+15.994915]PEPTIDE"
#' - Name to unimodId: "M[Oxidation]PEPTIDE" -> "M[UNIMOD:35]PEPTIDE"
#' - deltaMass to name: "M[+15.995]PEPTIDE" -> "M[Oxidation]PEPTIDE"
#' - deltaMass to unimodId: "M[+15.995]PEPTIDE" -> "M[UNIMOD:35]PEPTIDE"
#' - unimodId to name: "M[UNIMOD:35]PEPTIDE" -> "M[Oxidation]PEPTIDE"
#' - unimodId to deltaMass: "M[UNIMOD:35]PEPTIDE" -> "M[+15.994915]PEPTIDE"
#'
#' @examples
#' # Convert sequence from name to delta mass
#' convertAnnotation("M[Oxidation]PEPTIDE", convertToStyle = "deltaMass")
#' # Result: "M[+15.994915]PEPTIDE"
#'
#' # Name to Unimod ID
#' convertAnnotation("M[Oxidation]PEPTIDE", convertToStyle = "unimodId")
#' # Result: "M[UNIMOD:35]PEPTIDE"
#'
#' # Delta mass to name
#' convertAnnotation("M[+15.995]PEPTIDE", convertToStyle = "name")
#' # Result: "M[Oxidation]PEPTIDE"
#'
#' # Multiple modifications
#' convertAnnotation("M[Oxidation]EVNES[Phospho]PEK", convertToStyle = "deltaMass")
#' # Result: "M[+15.994915]EVNES[+79.966331]PEK"
#'
#' # Convert multiple sequences from name to delta mass
#' sequences <- c("M[Oxidation]PEPTIDE", "EVNES[Phospho]PEK", "PEPTIDE")
#' convertAnnotation(sequences, convertToStyle = "deltaMass")
#' # Result: c("M[+15.994915]PEPTIDE", "EVNES[+79.966331]PEK", "PEPTIDE")
#'
#' # Convert from delta mass to name
#' sequences <- c("M[+15.995]PEPTIDE", "S[+79.966]EQUENCE")
#' convertAnnotation(sequences, convertToStyle = "name")
#' # Result: c("M[Oxidation]PEPTIDE", "S[Phospho]EQUENCE")
#'
#' # Convert to Unimod IDs
#' sequences <- c("M[Oxidation]PEPTIDE", "C[Carbamidomethyl]PEPTIDE")
#' convertAnnotation(sequences, convertToStyle = "unimodId")
#' # Result: c("M[UNIMOD:35]PEPTIDE", "C[UNIMOD:4]PEPTIDE")
#'
#' @export
convertAnnotation <- function(x,
                              convertToStyle = c("deltaMass", "unimodId", "name"),
                              massTolerance = 0.01) {

    # Validate inputs
    if (!is.character(x)) {
        stop("x must be a character vector")
    }

    convertToStyle <- match.arg(convertToStyle)

    # Apply .convertAnnotation to each element
    vapply(x,
           function(seq) .convertAnnotation(seq,
                                           convertToStyle = convertToStyle,
                                           massTolerance = massTolerance),
           character(1),
           USE.NAMES = FALSE)
}


#' Converter for a character of length(1L)
#'
#' @noRd
.convertAnnotation <- function(x,
                              convertToStyle = c("deltaMass", "unimodId", "name"),
                              massTolerance = 0.01) {
    if (!is.character(x) || length(x) != 1L) {
        stop("x must be a single character string")
    }

    convertToStyle <- match.arg(convertToStyle)

    # Get Unimod data - simplified without priority rules
    unimodData <- modifications[!modifications$NeutralLoss,
                                         c("UnimodId", "Name", "MonoMass")]
    unimodData <- unimodData[!duplicated(unimodData$Name), ]

    # Find all modifications in the sequence
    # Pattern matches: [content] where content is not empty
    pattern <- "\\[([^]]+)\\]"

    # Extract all modifications
    matches <- gregexpr(pattern, x)
    mod_strings <- regmatches(x, matches)[[1]]

    if (length(mod_strings) == 0L) {
        # No modifications found, return unchanged
        return(x)
    }

    # Extract content within brackets
    modContents <- gsub("\\[|\\]", "", mod_strings)

    # Convert each modification
    result <- x
    for (i in seq_along(modContents)) {
        modContent <- modContents[i]

        # Determine the input type
        inputType <- .detectModificationType(modContent)

        # Skip if already in target format
        if (inputType == convertToStyle) {
            next
        }

        # Perform conversion
        converted <- .convertModificationFormat(
            modContent = modContent,
            inputType = inputType,
            outputType = convertToStyle,
            unimodData = unimodData,
            massTolerance = massTolerance
        )

        # Replace in the result string if conversion successful
        if (!is.null(converted)) {
            old_string <- mod_strings[i]
            new_string <- paste0("[", converted, "]")
            result <- sub(old_string, new_string, result, fixed = TRUE)
        }
    }

    return(result)
}


#' Detect modification type from string
#'
#' Determines whether a modification string represents a name, mass shift,
#' or Unimod ID.
#'
#' @param mod_string Character string with modification (without brackets)
#'
#' @return Character string: "name", "deltaMass", or "unimodId"
#'
#' @examples
#' # Detect Unimod ID format
#' .detectModificationType("UNIMOD:35")
#' # Result: "unimodId"
#'
#' # Detect delta mass format
#' .detectModificationType("+15.995")
#' # Result: "deltaMass"
#'
#' # Detect name format
#' .detectModificationType("Oxidation")
#' # Result: "name"
#'
#' @noRd
.detectModificationType <- function(mod_string) {
    # Check for Unimod ID (i.e. UNIMOD:35 or U:35)
    if (grepl("^(UNIMOD|U):[0-9]+$", mod_string, ignore.case = TRUE)) {
        return("unimodId")
    }

    # Check for mass shift (+15.995 or -17.026)
    if (grepl("^[+-][0-9]+(\\.[0-9]+)?$", mod_string)) {
        return("deltaMass")
    }

    # Otherwise assume it's a name
    return("name")
}


#' Convert Modification Between Formats
#'
#' Converts a single modification between different annotation formats using
#' the Unimod database.
#'
#' @param modContent Character string with modification content (without brackets)
#'
#' @param inputType Character string: "name", "deltaMass", or "unimodId"
#'
#' @param outputType Character string: "name", "deltaMass", or "unimodId"
#'
#' @param unimodData Data.frame with Unimod data
#'
#' @param massTolerance Numeric mass tolerance in Daltons
#'
#' @return Character string with converted modification, or NULL if conversion failed
#'
#' @examples
#' # Load Unimod data
#' unimodData <- PTMods::modifications[!PTMods::modifications$NeutralLoss,
#'                                      c("UnimodId", "Name", "MonoMass")]
#' unimodData <- unimodData[!duplicated(unimodData$Name), ]
#'
#' # Convert from name to delta mass
#' .convertModificationFormat("Oxidation", "name", "deltaMass",
#'                           unimodData, 0.01)
#' # Result: "+15.994915"
#'
#' # Convert from delta mass to name
#' .convertModificationFormat("+15.995", "deltaMass", "name",
#'                           unimodData, 0.01)
#' # Result: "Oxidation"
#'
#' @noRd
.convertModificationFormat <- function(modContent, inputType, outputType,
                                       unimodData, massTolerance) {

    # Strategy: always convert to intermediate representation first
    # inputType → intermediate (UnimodId + Mass + Name) → outputType

    intermediate <- NULL

    # Step 1: Convert input to intermediate representation
    if (inputType == "name") {
        # Look up by name
        intermediate <- .lookupByName(modContent, unimodData)

    } else if (inputType == "deltaMass") {
        # Extract mass value
        mass_value <- as.numeric(modContent)
        # Look up by mass
        intermediate <- .lookupByMass(mass_value, unimodData, massTolerance)

    } else if (inputType == "unimodId") {
        # Extract Unimod ID
        unimodId <- .extractUnimodId(modContent)
        # Look up by ID
        intermediate <- .lookupByUnimodId(unimodId, unimodData)
    }

    # Check if lookup was successful
    if (is.null(intermediate)) {
        warning(paste0(
            "Could not find Unimod entry for modification ",
            modContent, ", see `?modifications`"))
        return(NULL)
    }

    # Step 2: Convert intermediate to output format
    if (outputType == "name") {
        return(intermediate$name)

    } else if (outputType == "deltaMass") {
        # Format mass with appropriate precision (by default 6 decimals)
        return(.formatMass(intermediate$mass))

    } else if (outputType == "unimodId") {
        return(paste0("UNIMOD:", intermediate$unimodId))
    }

    return(NULL)
}


#' Look Up Modification by Name
#'
#' Finds a modification in the Unimod database by its name.
#'
#' @param name Character string with modification name
#'
#' @param unimodData Data.frame with Unimod data
#'
#' @return List with unimodId, name, and mass, or NULL if not found
#'
#' @examples
#' # Load Unimod data
#' unimodData <- PTMods::modifications[!PTMods::modifications$NeutralLoss,
#'                                      c("UnimodId", "Name", "MonoMass")]
#' unimodData <- unimodData[!duplicated(unimodData$Name), ]
#'
#' # Look up by name
#' .lookupByName("Oxidation", unimodData)
#' # Result: list(unimodId = 35, name = "Oxidation", mass = 15.994915)
#'
#' # Look up non-existent modification
#' .lookupByName("NonExistent", unimodData)
#' # Result: NULL
#'
#' @noRd
.lookupByName <- function(name, unimodData) {
    # Match by name (exact match)
    match_idx <- which(unimodData$Name == name)

    if (length(match_idx) == 0L) {
        return(NULL)
    }

    # Return the first match
    list(
        unimodId = unimodData$UnimodId[match_idx[1]],
        name = unimodData$Name[match_idx[1]],
        mass = unimodData$MonoMass[match_idx[1]]
    )
}


#' Look Up Modification by Mass
#'
#' Finds a modification in the Unimod database by its monoisotopic mass.
#'
#' @param mass Numeric mass value
#'
#' @param unimodData Data.frame with Unimod data
#'
#' @param tolerance Numeric mass tolerance in Daltons
#'
#' @return List with unimodId, name, and mass, or NULL if not found
#'
#' @examples
#' # Load Unimod data
#' unimodData <- PTMods::modifications[!PTMods::modifications$NeutralLoss,
#'                                      c("UnimodId", "Name", "MonoMass")]
#' unimodData <- unimodData[!duplicated(unimodData$Name), ]
#'
#' # Look up by mass with tolerance
#' .lookupByMass(15.995, unimodData, 0.01)
#' # Result: list(unimodId = 35, name = "Oxidation", mass = 15.994915)
#'
#' # Look up with no match
#' .lookupByMass(999.999, unimodData, 0.01)
#' # Result: NULL
#'
#' @noRd
.lookupByMass <- function(mass, unimodData, tolerance) {
    # Match by mass within tolerance
    mass_diff <- abs(unimodData$MonoMass - mass)
    match_idx <- which(mass_diff <= tolerance)

    if (length(match_idx) == 0L) {
        return(NULL)
    }

    # If multiple matches, return the closest one
    if (length(match_idx) > 1L) {
        closest_idx <- match_idx[which.min(mass_diff[match_idx])]
        match_idx <- closest_idx
    }

    # Return the best match
    list(
        unimodId = unimodData$UnimodId[match_idx[1]],
        name = unimodData$Name[match_idx[1]],
        mass = unimodData$MonoMass[match_idx[1]]
    )
}


#' Look Up Modification by Unimod ID
#'
#' Finds a modification in the Unimod database by its Unimod ID.
#'
#' @param unimodId Numeric Unimod ID
#'
#' @param unimodData Data.frame with Unimod data
#'
#' @return List with unimodId, name, and mass, or NULL if not found
#'
#' @examples
#' # Load Unimod data
#' unimodData <- PTMods::modifications[!PTMods::modifications$NeutralLoss,
#'                                      c("UnimodId", "Name", "MonoMass")]
#' unimodData <- unimodData[!duplicated(unimodData$Name), ]
#'
#' # Look up by Unimod ID
#' .lookupByUnimodId(35, unimodData)
#' # Result: list(unimodId = 35, name = "Oxidation", mass = 15.994915)
#'
#' # Look up non-existent ID
#' .lookupByUnimodId(99999, unimodData)
#' # Result: NULL
#'
#' @noRd
.lookupByUnimodId <- function(unimodId, unimodData) {
    # Match by UnimodId
    match_idx <- which(unimodData$UnimodId == unimodId)

    if (length(match_idx) == 0L) {
        return(NULL)
    }

    # Return the match
    list(
        unimodId = unimodData$UnimodId[match_idx[1]],
        name = unimodData$Name[match_idx[1]],
        mass = unimodData$MonoMass[match_idx[1]]
    )
}


#' Extract Unimod ID from String
#'
#' Extracts the numeric Unimod ID from strings like "UNIMOD:35" or "U:35".
#'
#' @param unimod_string Character string with Unimod ID
#'
#' @return Numeric Unimod ID
#'
#' @examples
#' # Extract from full format
#' .extractUnimodId("UNIMOD:35")
#' # Result: 35
#'
#' # Extract from short format
#' .extractUnimodId("U:35")
#' # Result: 35
#'
#' # Case insensitive
#' .extractUnimodId("unimod:35")
#' # Result: 35
#'
#' @noRd
.extractUnimodId <- function(unimod_string) {
    # Extract the number after "UNIMOD:" or "U:"
    id_str <- sub("^(UNIMOD|U):", "", unimod_string, ignore.case = TRUE)
    as.integer(id_str)
}


#' Format Mass Value
#'
#' Formats a mass value as a string with appropriate sign and precision.
#'
#' @param mass Numeric mass value
#'
#' @param digits Integer number of decimal places (default: 6)
#'
#' @return Character string with formatted mass (e.g., "+15.994915")
#'
#' @examples
#' # Format positive mass
#' .formatMass(15.994915)
#' # Result: "+15.994915"
#'
#' # Format negative mass
#' .formatMass(-17.026549)
#' # Result: "-17.026549"
#'
#' # Custom precision
#' .formatMass(15.994915, digits = 3)
#' # Result: "+15.995"
#'
#' @noRd
.formatMass <- function(mass, digits = 6) {
    # Format with sign
    if (mass >= 0) {
        paste0("+", round(mass, digits))
    } else {
        as.character(round(mass, digits))
    }
}