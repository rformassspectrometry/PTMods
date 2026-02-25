#' @title Combines consecutive modifications on amino acids in peptide
#' sequences.
#'
#' @param sequences `character`. A character vector of peptide sequences with
#' modifications in any format accepted by `convertAnnotation`.
#' Multiple modifications per amino acid are allowed.
#'
#' @return A character vector of the same length as `sequences` with
#' consecutive modifications on each amino acid summed into a single
#' modification. This summed modifications is written in delta mass format.
#'
#' @author Guillaume Deflandre <guillaume.deflandre@uclouvain.be>
#'
#' @examples
#' combineModifications("AT[+10][-5]K")
#' combineModifications("[+10]-[-5]-ATK")
#' combineModifications("AT[Phospho][UNIMOD:1]K")
#' combineModifications("EM[+15.9949][-10]EK[+42][-8]")
#' combineModifications(c("AT[+10][-5]K", "EM[+15.9949][-10]EK"))
#'
#' @export
combineModifications <- function(sequences) {
    sequences <- convertAnnotation(sequences)

    sapply(sequences, function(sequence) {
        sequence <- .combineTerminalMods(sequence)
        terms <- .extractTerminalMods(sequence)
        sequence <- terms$sequence

        masses <- .parseModifiedSequence(sequence)
        rx_residue <- "[A-Z](?:\\[(?:[GMURX]:)?[+-][0-9.]+\\])*"
        tokens <- regmatches(
            sequence,
            gregexpr(rx_residue, sequence, perl = TRUE)
        )[[1L]]
        aa <- substring(tokens, 1L, 1L)

        inner <- paste(
            ifelse(masses != 0,
                paste0(aa, "[", ifelse(masses > 0, "+", ""), masses, "]"),
                aa
            ),
            collapse = ""
        )

        .reattachTerminalMods(inner, terms$nterm, terms$cterm)
    }, USE.NAMES = FALSE)
}

#' @title Generates a numeric with positioned modifications
#'
#' @param sequence `character(1L)`. A peptide sequence that may have
#' modifications in delta mass format. Multiple modifications per
#' amino acid are allowed.
#'
#' @return `numeric()` of length equal to canonical sequence with specified
#' modifications.
#'
#' @author Guillaume Deflandre <guillaume.deflandre@uclouvain.be>
#'
#' @noRd
#'
#' @examples
#' PTMods:::.parseModifiedSequence("APEPT[+79.966]IDEK")
#' PTMods:::.parseModifiedSequence("[+304]-APEPTIDEK")
#' PTMods:::.parseModifiedSequence("EM[+15.9949]EVEES[+79.9663]PEK")
#' PTMods:::.parseModifiedSequence("APEPT[+79.966][+10]IDEK")
#' PTMods:::.parseModifiedSequence("EM[+15.9949][-10]EK")
.parseModifiedSequence <- function(sequence) {
    ## Match each amino acid followed by zero or more modification brackets
    rx_residue <- "[A-Z](?:\\[(?:[GMURX]:)?[+-][0-9.]+\\])*"
    ## Match delta masses: sign followed by digits, inside brackets
    rx_mass <- "[+-][0-9.]+"

    tokens <- regmatches(
        sequence,
        gregexpr(rx_residue, sequence, perl = TRUE)
    )[[1L]]

    masses <- vapply(tokens,
        function(token) {
            ## Remove leading amino acid character before extracting masses
            bracket_section <- substring(token, 2L)
            masses <- regmatches(
                bracket_section,
                gregexpr(rx_mass, bracket_section, perl = TRUE)
            )[[1L]]

            if (length(masses) == 0L) {
                return(0)
            }
            sum(as.double(masses))
        },
        numeric(1L),
        USE.NAMES = FALSE
    )

    ## Extract terminal modification values as attributes
    terms <- .extractTerminalMods(sequence)

    .extractTerminalMass <- function(term) {
        if (!nchar(term)) {
            return(NA_real_)
        }
        masses <- regmatches(term, gregexpr(rx_mass, term, perl = TRUE))[[1L]]
        sum(as.numeric(masses))
    }

    attr(masses, "Nterm") <- .extractTerminalMass(terms$nterm)
    attr(masses, "Cterm") <- .extractTerminalMass(terms$cterm)

    masses
}

#' @title Combine stacked terminal modifications into a single bracket.
#'
#' @param sequence `character(1L)`. A peptide sequence possibly containing
#' stacked N-terminal or C-terminal modifications.
#'
#' @return `character(1L)` with each terminal stack summed into a single
#' bracketed delta mass.
#'
#' @noRd
#'
#' @examples
#' .combineTerminalMods("[+304]-[-54]-ATK")
#' .combineTerminalMods("ATK-[+42]-[-8]")
#' .combineTerminalMods("[+304]-[-54]-ATK-[+42]-[-8]")
#' .combineTerminalMods("ATK")
.combineTerminalMods <- function(sequence) {
    rx_mass <- "[+-][0-9.]+"

    ## N-terminal: one or more [±x]- prefixes
    nterm_match <- regmatches(
        sequence,
        regexpr("^(\\[[^\\]]*\\]-)+", sequence, perl = TRUE)
    )

    if (length(nterm_match) && nchar(nterm_match) > 0L) {
        masses <- as.numeric(regmatches(
            nterm_match,
            gregexpr(rx_mass, nterm_match, perl = TRUE)
        )[[1L]])
        total <- sum(masses)
        sequence <- paste0(
            "[", ifelse(total > 0, "+", ""), total, "]-",
            substring(sequence, nchar(nterm_match) + 1L)
        )
    }

    ## C-terminal: one or more -[±x] suffixes
    cterm_match <- regmatches(
        sequence,
        regexpr("(-\\[[^\\]]*\\])+$", sequence, perl = TRUE)
    )

    if (length(cterm_match) && nchar(cterm_match) > 0L) {
        masses <- as.numeric(regmatches(
            cterm_match,
            gregexpr(rx_mass, cterm_match, perl = TRUE)
        )[[1L]])
        total <- sum(masses)
        sequence <- paste0(
            substring(sequence, 1L, nchar(sequence) - nchar(cterm_match)),
            "-[", ifelse(total > 0, "+", ""), total, "]"
        )
    }

    sequence
}
