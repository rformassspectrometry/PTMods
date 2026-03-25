#' @title Add fixed modifications
#'
#' @description Applies fixed modifications to peptide sequences.
#'
#' The annotation style of modifications already present in `sequences` is
#' preserved as-is.  The annotation style of newly applied modifications
#' matches the style supplied in `fixedModifications`: numeric values produce
#' deltaMass notation (e.g. `[+79.966331]`), while character values are used
#' verbatim (e.g. `[Phospho]` or `[UNIMOD:21]`).
#'
#' @param sequences Character vector. Peptide sequences that may have
#' modifications or not.
#'
#' @param fixedModifications Named `numeric` or `character`. If a
#' `character` is given, values must be in UniMod name or UniMod ID format
#' (e.g. `"Phospho"`, `"UNIMOD:21"`). The annotation style of the values is
#' preserved in the output.
#' Specifies which fixed modifications are applied to which amino acids. By
#' default, carbamidomethylation is applied (+57.021464 on cysteines).
#' Supplying multiple modifications for the same amino acid is not supported;
#' call the function consecutively instead:
#' \preformatted{
#'   "ATK" |>
#'     addFixedModifications(c(T = "+57.024")) |>
#'     addFixedModifications(c(T = "Phospho"))
#' }
#'
#' @return A named list of character strings with all fixed modifications
#' applied to the corresponding amino acids, one element per input sequence.
#'
#' @author Guillaume Deflandre <guillaume.deflandre@uclouvain.be>
#'
#' @export
#'
#' @examples
#' addFixedModifications("ATCK")
#' addFixedModifications("ATCK", fixedModifications = c(Nterm = 304, C = 57.02146))
#' addFixedModifications("ATSK", fixedModifications = c(A = -42, S = 79.966331))
#' addFixedModifications("ATK", fixedModifications = c(T = "Phospho", A = 42))
#' addFixedModifications("A[+42]TK", fixedModifications = c(T = "Phospho"))
#' addFixedModifications("AT[+79.966]AK", fixedModifications = c(A = 42))
#' addFixedModifications("ATA[+79.966]K", fixedModifications = c(A = 42))
#'
#' ## Apply two modifications to the same residue consecutively
#' "ATK" |>
#'     addFixedModifications(c(T = "+57.024")) |>
#'     addFixedModifications(c(T = "Phospho"))
addFixedModifications <- function(sequences,
                                  fixedModifications = c(C = 57.021464)) {
    unname(sapply(sequences, .addFixedModifications,
        fixedModifications = fixedModifications
    ))
}

#' @title Applies fixed modifications to a peptide sequence
#'
#' @param sequence Character. A peptide sequence that may have modifications or
#' not.
#'
#' @param fixedModifications Named `numeric` or `character`. If a
#' `character` is given, it needs to be in unimodId or name format as
#' specified in `convertAnnotation`.
#' Specifies which fixed modifications are applied to which amino acids.
#' Use Nterm or Cterm as names for modifications that should be added to the
#' amino/carboxyl-termini. Default is set to carbamidomethylation
#' (C = 57.02146). Supplying multiple modifications for the same amino acid
#' is not supported; call `addFixedModifications()` consecutively instead:
#' \preformatted{
#'   "ATK" |>
#'     addFixedModifications(c(T = "+57.024")) |>
#'     addFixedModifications(c(T = "Phospho"))
#' }
#'
#' @return A character string with all fixed modifications applied to the
#' corresponding amino acids.
#'
#' @author Guillaume Deflandre <guillaume.deflandre@uclouvain.be>
#'
#' @noRd
#'
#' @examples
#' .addFixedModifications("ATCK", fixedModifications = c(Nterm = 304, C = 57.02146))
#' .addFixedModifications("ATSK", fixedModifications = c(A = -42, S = 79.966331))
#' .addFixedModifications("ATK", fixedModifications = c(T = "Phospho", A = "UNIMOD:1"))
#' .addFixedModifications("AT[+79.966]AK", fixedModifications = c(A = 42))
#' .addFixedModifications("ATA[+79.966]K", fixedModifications = c(A = 42))
#'
#' ## Apply two modifications to the same residue consecutively
#' "ATK" |>
#'     addFixedModifications(c(T = "+57.024")) |>
#'     addFixedModifications(c(T = "Phospho"))
.addFixedModifications <- function(sequence,
                                   fixedModifications = c(C = 57.02146)) {
    if (!is.null(fixedModifications)) {
        mod_table <- table(names(fixedModifications))

        if (!all(mod_table == 1)) {
            stop(
                "Applying multiple modifications to the same amino acid is not ",
                "supported.\n",
                "Please use addFixedModifications() consecutively instead:\n\n",
                "  'ATK' |>\n",
                "    addFixedModifications(c(T = '+57.024')) |>\n",
                "    addFixedModifications(c(T = 'Phospho'))"
            )
        }
    }

    terms <- .extractTerminalMods(sequence)
    sequence <- terms$sequence

    sequence <- regmatches(
        sequence,
        gregexpr("[A-Z](\\[[^\\]]*\\])*", sequence, perl = TRUE)
    )[[1]]
    seqAA <- sub("\\[.*", "", sequence)
    modifiedAA <- names(fixedModifications)

    fixedMods <- .resolveModsToStrings(fixedModifications)
    names(fixedMods) <- modifiedAA

    mStr <- setNames(rep(NA_character_, length(sequence)), sequence)
    applicable <- which(seqAA %in% modifiedAA)
    if (length(applicable)) {
        mStr[applicable] <- fixedMods[seqAA[applicable]]
    }

    sequence <- paste(
        ifelse(!is.na(mStr),
            paste0(names(mStr), "[", mStr, "]"),
            names(mStr)
        ),
        collapse = ""
    )

    sequence <- .addTerminalModifications(sequence, fixedMods)
    .reattachTerminalMods(sequence, terms$nterm, terms$cterm)
}

#' Adds nterm/cterm modifications to sequence
#'
#' @param sequence A `character(1L)` with a peptide sequence with or without
#' modifications.
#'
#' @param modifications A named `numeric` or `character` vector with all
#' modifications to be applied on the peptide sequence. Numeric values are
#' formatted as deltaMass (e.g. `[+67]`); character values are used verbatim
#' (e.g. `[Phospho]` or `[UNIMOD:21]`).
#'
#' @return The same sequence with additional Nterm or Cterm modifications if
#' present in `modifications`.
#'
#' @noRd
#'
#' @examples
#' PTMods:::.addTerminalModifications("ATK", c(Nterm = 67))
#' PTMods:::.addTerminalModifications("ATK", c(Nterm = "Acetyl"))
.addTerminalModifications <- function(sequence, modifications) {
    modStrings <- .resolveModsToStrings(modifications)

    if ("Nterm" %in% names(modStrings)) {
        val <- modStrings[["Nterm"]]
        sequence <- paste0("[", val, "]-", sequence)
    }

    if ("Cterm" %in% names(modStrings)) {
        val <- modStrings[["Cterm"]]
        sequence <- paste0(sequence, "-[", val, "]")
    }

    sequence
}
