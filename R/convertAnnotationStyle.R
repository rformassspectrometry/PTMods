#' @name convertAnnotation
#' 
#' @title Convert sequences from one annotation style to another
#'
#' @description
#' Converts modifications between different annotation formats for multiple 
#' sequences at once. See the details and examples sections for more information.
#' The annotation styles are inferred from the `modifications` dataframe. 
#'
#' @param x Character vector with peptide sequences in ProForma format
#'
#' @param convert_to Character string specifying target format. Options:
#'   "delta_mass", "unimod_id", "name"
#'
#' @param mass_tolerance Numeric mass tolerance in Daltons for matching 
#'   modifications (default: 0.01). Used when converting from delta_mass.
#'
#' @return Character vector with the sequences in the target annotation format
#'
#' @author Guillaume Deflandre <guillaume.deflandre@uclouvain.be>
#' 
#' @details
#' The function handles three main conversion scenarios:
#' \itemize{
#'   \item Name to delta_mass: "M[Oxidation]PEPTIDE" -> "M[+15.994915]PEPTIDE"
#'   \item Name to unimod_id: "M[Oxidation]PEPTIDE" -> "M[UNIMOD:35]PEPTIDE"
#'   \item Delta_mass to name: "M[+15.995]PEPTIDE" -> "M[Oxidation]PEPTIDE"
#'   \item Delta_mass to unimod_id: "M[+15.995]PEPTIDE" -> "M[UNIMOD:35]PEPTIDE"
#'   \item Unimod_id to name: "M[UNIMOD:35]PEPTIDE" -> "M[Oxidation]PEPTIDE"
#'   \item Unimod_id to delta_mass: "M[UNIMOD:35]PEPTIDE" -> "M[+15.994915]PEPTIDE"
#' }
#' 
#' @examples
#' # Convert sequence from name to delta mass
#' convertAnnotation("M[Oxidation]PEPTIDE", convert_to = "delta_mass")
#' # Result: "M[+15.994915]PEPTIDE"
#' 
#' # Name to Unimod ID
#' convertAnnotation("M[Oxidation]PEPTIDE", convert_to = "unimod_id")
#' # Result: "M[UNIMOD:35]PEPTIDE"
#' 
#' # Delta mass to name
#' convertAnnotation("M[+15.995]PEPTIDE", convert_to = "name")
#' # Result: "M[Oxidation]PEPTIDE"
#' 
#' # Multiple modifications
#' convertAnnotation("M[Oxidation]EVNES[Phospho]PEK", convert_to = "delta_mass")
#' # Result: "M[+15.994915]EVNES[+79.966331]PEK"
#' 
#' # Convert multiple sequences from name to delta mass
#' sequences <- c("M[Oxidation]PEPTIDE", "EVNES[Phospho]PEK", "PEPTIDE")
#' convertAnnotation(sequences, convert_to = "delta_mass")
#' # Result: c("M[+15.994915]PEPTIDE", "EVNES[+79.966331]PEK", "PEPTIDE")
#' 
#' # Convert from delta mass to name
#' sequences <- c("M[+15.995]PEPTIDE", "S[+79.966]EQUENCE")
#' convertAnnotation(sequences, convert_to = "name")
#' # Result: c("M[Oxidation]PEPTIDE", "S[Phospho]EQUENCE")
#' 
#' # Convert to Unimod IDs
#' sequences <- c("M[Oxidation]PEPTIDE", "C[Carbamidomethyl]PEPTIDE")
#' convertAnnotation(sequences, convert_to = "unimod_id")
#' # Result: c("M[UNIMOD:35]PEPTIDE", "C[UNIMOD:4]PEPTIDE")
#'
#' @export
convertAnnotation <- function(x, 
                              convert_to = c("delta_mass", "unimod_id", "name"),
                              mass_tolerance = 0.01) {
    
    # Validate inputs
    if (!is.character(x)) {
        stop("x must be a character vector")
    }
    
    convert_to <- match.arg(convert_to)
    
    # Apply .convertAnnotation to each element
    vapply(x, 
           function(seq) .convertAnnotation(seq, 
                                           convert_to = convert_to, 
                                           mass_tolerance = mass_tolerance),
           character(1),
           USE.NAMES = FALSE)
}


#' Converter for a character of length(1L)
#' 
#' @noRd
.convertAnnotation <- function(x, 
                              convert_to = c("delta_mass", "unimod_id", "name"),
                              mass_tolerance = 0.01) {
    
    if (!is.character(x) || length(x) != 1L) {
        stop("x must be a single character string")
    }
    
    convert_to <- match.arg(convert_to)
    
    # Get Unimod data - simplified without priority rules
    unimod_data <- unimod::modifications[!unimod::modifications$NeutralLoss, 
                                         c("UnimodId", "Name", "MonoMass")]
    unimod_data <- unimod_data[!duplicated(unimod_data$Name), ]
    
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
    mod_contents <- gsub("\\[|\\]", "", mod_strings)
    
    # Convert each modification
    result <- x
    for (i in seq_along(mod_contents)) {
        mod_content <- mod_contents[i]
        
        # Determine the input type
        input_type <- .detectModificationType(mod_content)
        
        # Skip if already in target format
        if (input_type == convert_to) {
            next
        }
        
        # Perform conversion
        converted <- .convertModificationFormat(
            mod_content = mod_content,
            input_type = input_type,
            output_type = convert_to,
            unimod_data = unimod_data,
            mass_tolerance = mass_tolerance
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
#' @return Character string: "name", "delta_mass", or "unimod_id"
#'
#' @examples
#' # Detect Unimod ID format
#' .detectModificationType("UNIMOD:35")
#' # Result: "unimod_id"
#' 
#' # Detect delta mass format
#' .detectModificationType("+15.995")
#' # Result: "delta_mass"
#' 
#' # Detect name format
#' .detectModificationType("Oxidation")
#' # Result: "name"
#'
#' @noRd
.detectModificationType <- function(mod_string) {
    # Check for Unimod ID (i.e. UNIMOD:35 or U:35)
    if (grepl("^(UNIMOD|U):[0-9]+$", mod_string, ignore.case = TRUE)) {
        return("unimod_id")
    }
    
    # Check for mass shift (+15.995 or -17.026)
    if (grepl("^[+-][0-9]+(\\.[0-9]+)?$", mod_string)) {
        return("delta_mass")
    }
    
    # Otherwise assume it's a name
    return("name")
}


#' Convert Modification Between Formats
#'
#' Converts a single modification between different annotation formats using
#' the Unimod database.
#'
#' @param mod_content Character string with modification content (without brackets)
#'
#' @param input_type Character string: "name", "delta_mass", or "unimod_id"
#'
#' @param output_type Character string: "name", "delta_mass", or "unimod_id"
#'
#' @param unimod_data Data.frame with Unimod data
#'
#' @param mass_tolerance Numeric mass tolerance in Daltons
#'
#' @return Character string with converted modification, or NULL if conversion failed
#'
#' @examples
#' # Load Unimod data
#' unimod_data <- unimod::modifications[!unimod::modifications$NeutralLoss, 
#'                                      c("UnimodId", "Name", "MonoMass")]
#' unimod_data <- unimod_data[!duplicated(unimod_data$Name), ]
#' 
#' # Convert from name to delta mass
#' .convertModificationFormat("Oxidation", "name", "delta_mass", 
#'                           unimod_data, 0.01)
#' # Result: "+15.994915"
#' 
#' # Convert from delta mass to name
#' .convertModificationFormat("+15.995", "delta_mass", "name", 
#'                           unimod_data, 0.01)
#' # Result: "Oxidation"
#'
#' @noRd
.convertModificationFormat <- function(mod_content, input_type, output_type, 
                                       unimod_data, mass_tolerance) {
    
    # Strategy: always convert to intermediate representation first
    # input_type → intermediate (UnimodId + Mass + Name) → output_type
    
    intermediate <- NULL
    
    # Step 1: Convert input to intermediate representation
    if (input_type == "name") {
        # Look up by name
        intermediate <- .lookupByName(mod_content, unimod_data)
        
    } else if (input_type == "delta_mass") {
        # Extract mass value
        mass_value <- as.numeric(mod_content)
        # Look up by mass
        intermediate <- .lookupByMass(mass_value, unimod_data, mass_tolerance)
        
    } else if (input_type == "unimod_id") {
        # Extract Unimod ID
        unimod_id <- .extractUnimodId(mod_content)
        # Look up by ID
        intermediate <- .lookupByUnimodId(unimod_id, unimod_data)
    }
    
    # Check if lookup was successful
    if (is.null(intermediate)) {
        warning(paste0(
            "Could not find Unimod entry for modification ", mod_content))
        return(NULL)
    }
    
    # Step 2: Convert intermediate to output format
    if (output_type == "name") {
        return(intermediate$name)
        
    } else if (output_type == "delta_mass") {
        # Format mass with appropriate precision (by default 6 decimals)
        return(.formatMass(intermediate$mass))
        
    } else if (output_type == "unimod_id") {
        return(paste0("UNIMOD:", intermediate$unimod_id))
    }
    
    return(NULL)
}


#' Look Up Modification by Name
#'
#' Finds a modification in the Unimod database by its name.
#'
#' @param name Character string with modification name
#'
#' @param unimod_data Data.frame with Unimod data
#'
#' @return List with unimod_id, name, and mass, or NULL if not found
#'
#' @examples
#' # Load Unimod data
#' unimod_data <- unimod::modifications[!unimod::modifications$NeutralLoss, 
#'                                      c("UnimodId", "Name", "MonoMass")]
#' unimod_data <- unimod_data[!duplicated(unimod_data$Name), ]
#' 
#' # Look up by name
#' .lookupByName("Oxidation", unimod_data)
#' # Result: list(unimod_id = 35, name = "Oxidation", mass = 15.994915)
#' 
#' # Look up non-existent modification
#' .lookupByName("NonExistent", unimod_data)
#' # Result: NULL
#'
#' @noRd
.lookupByName <- function(name, unimod_data) {
    # Match by name (exact match)
    match_idx <- which(unimod_data$Name == name)
    
    if (length(match_idx) == 0L) {
        return(NULL)
    }
    
    # Return the first match
    list(
        unimod_id = unimod_data$UnimodId[match_idx[1]],
        name = unimod_data$Name[match_idx[1]],
        mass = unimod_data$MonoMass[match_idx[1]]
    )
}


#' Look Up Modification by Mass
#'
#' Finds a modification in the Unimod database by its monoisotopic mass.
#'
#' @param mass Numeric mass value
#'
#' @param unimod_data Data.frame with Unimod data
#'
#' @param tolerance Numeric mass tolerance in Daltons
#'
#' @return List with unimod_id, name, and mass, or NULL if not found
#'
#' @examples
#' # Load Unimod data
#' unimod_data <- unimod::modifications[!unimod::modifications$NeutralLoss, 
#'                                      c("UnimodId", "Name", "MonoMass")]
#' unimod_data <- unimod_data[!duplicated(unimod_data$Name), ]
#' 
#' # Look up by mass with tolerance
#' .lookupByMass(15.995, unimod_data, 0.01)
#' # Result: list(unimod_id = 35, name = "Oxidation", mass = 15.994915)
#' 
#' # Look up with no match
#' .lookupByMass(999.999, unimod_data, 0.01)
#' # Result: NULL
#'
#' @noRd
.lookupByMass <- function(mass, unimod_data, tolerance) {
    # Match by mass within tolerance
    mass_diff <- abs(unimod_data$MonoMass - mass)
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
        unimod_id = unimod_data$UnimodId[match_idx[1]],
        name = unimod_data$Name[match_idx[1]],
        mass = unimod_data$MonoMass[match_idx[1]]
    )
}


#' Look Up Modification by Unimod ID
#'
#' Finds a modification in the Unimod database by its Unimod ID.
#'
#' @param unimod_id Numeric Unimod ID
#'
#' @param unimod_data Data.frame with Unimod data
#'
#' @return List with unimod_id, name, and mass, or NULL if not found
#'
#' @examples
#' # Load Unimod data
#' unimod_data <- unimod::modifications[!unimod::modifications$NeutralLoss, 
#'                                      c("UnimodId", "Name", "MonoMass")]
#' unimod_data <- unimod_data[!duplicated(unimod_data$Name), ]
#' 
#' # Look up by Unimod ID
#' .lookupByUnimodId(35, unimod_data)
#' # Result: list(unimod_id = 35, name = "Oxidation", mass = 15.994915)
#' 
#' # Look up non-existent ID
#' .lookupByUnimodId(99999, unimod_data)
#' # Result: NULL
#'
#' @noRd
.lookupByUnimodId <- function(unimod_id, unimod_data) {
    # Match by UnimodId
    match_idx <- which(unimod_data$UnimodId == unimod_id)
    
    if (length(match_idx) == 0L) {
        return(NULL)
    }
    
    # Return the match
    list(
        unimod_id = unimod_data$UnimodId[match_idx[1]],
        name = unimod_data$Name[match_idx[1]],
        mass = unimod_data$MonoMass[match_idx[1]]
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
