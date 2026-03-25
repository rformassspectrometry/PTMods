#' @title Add variable modifications
#'
#' @description Generates all possible combinations of variable modifications
#' for peptide sequences.
#'
#' The annotation style of modifications already present in `sequences` is
#' preserved as-is.  The annotation style of newly applied modifications
#' matches the style supplied in `variableModifications`: numeric values
#' produce deltaMass notation (e.g. `[+79.966331]`), while character values
#' are used verbatim (e.g. `[Phospho]` or `[UNIMOD:21]`).
#'
#' @param sequences Character vector. Peptide sequences that may have
#' modifications or not.
#'
#' @param variableModifications Named `numeric` or `character`. If a
#' `character` is given, values must be in UniMod name or UniMod ID format
#' (e.g. `"Phospho"`, `"UNIMOD:21"`). The annotation style of the values is
#' preserved in the output.
#' Specifies which variable modifications are used on which amino acids.
#' Supplying multiple modifications for the same amino acid is not supported;
#' call the function consecutively instead:
#' \preformatted{
#'   seq |>
#'     addVariableModifications(c(T = "+57.024")) |>
#'     addVariableModifications(c(T = "Phospho"))
#' }
#'
#' @param maxMods Numeric. Indicates how many modifications can be applied
#' at once.
#'
#' @return A named list of character vectors with all possible combinations
#' of modifications, one element per input sequence.
#'
#' @author Guillaume Deflandre <guillaume.deflandre@uclouvain.be>
#'
#' @importFrom utils combn
#'
#' @export
#'
#' @examples
#' addVariableModifications(
#'     "APGHKA",
#'     variableModifications = c(A = 4, K = 5, S = 8),
#'     maxMods = 2
#' )
#' addVariableModifications(
#'     "APGHKA",
#'     variableModifications = c(A = 4, K = 5, S = 8),
#'     maxMods = 3
#' )
#' addVariableModifications(
#'     c("ATK", "PAK"),
#'     variableModifications = c(T = "Phospho", A = "UNIMOD:1")
#' )
#' addVariableModifications(
#'     "A[+10]KT[-5]",
#'     variableModifications = c(A = -20, T = "Phospho")
#' )
#'
#' ## Apply two modifications to the same residue consecutively
#' "ATK" |>
#'     addVariableModifications(c(T = "+57.024")) |>
#'     addVariableModifications(c(T = "Phospho"))
addVariableModifications <- function(sequences,
                                     variableModifications = NULL,
                                     maxMods = Inf) {
    unlist(lapply(sequences, .addVariableModifications,
        variableModifications = variableModifications,
        maxMods = maxMods
    ), use.names = FALSE)
}

#' @title Generates list of possible combinations of modifications
#'
#' @param sequence Character. A peptide sequence that may have modifications or
#' not.
#'
#' @param variableModifications Named `numeric` or `character`. If a
#' `character` is given, it needs to be in unimodId or name format as
#' specified in `convertAnnotation`.
#' Specifies which variable modifications are used on which amino acids.
#' Supplying multiple modifications for the same amino acid is not supported;
#' call `addVariableModifications()` consecutively instead:
#' \preformatted{
#'   seq |>
#'     addVariableModifications(c(T = "+57.024")) |>
#'     addVariableModifications(c(T = "Phospho"))
#' }
#'
#' @param maxMods Numeric. Indicates how many modifications can be applied at
#' once.
#'
#' @return A `character()` vector with all possible combinations of
#' modifications.
#'
#' @author Guillaume Deflandre <guillaume.deflandre@uclouvain.be>
#'
#' @importFrom utils combn
#'
#' @noRd
#'
#' @examples
#' .addVariableModifications(
#'     "APGHKA",
#'     variableModifications = c(A = 4, K = 5, S = 8),
#'     maxMods = 2
#' )
#' .addVariableModifications(
#'     "APGHKA",
#'     variableModifications = c(A = 4, K = 5, S = 8),
#'     maxMods = 3
#' )
#' .addVariableModifications(
#'     "ATK",
#'     variableModifications = c(T = "Phospho", A = "UNIMOD:1")
#' )
#' .addVariableModifications(
#'     "A[+10]KT[-5]",
#'     variableModifications = c(A = -5, T = "Phospho")
#' )
#'
#' ## Apply two modifications to the same residue consecutively
#' "ATK" |>
#'     addVariableModifications(c(T = "+57.024")) |>
#'     addVariableModifications(c(T = "Phospho"))
.addVariableModifications <- function(sequence,
                                      variableModifications = NULL,
                                      maxMods = Inf) {
    if (!is.null(variableModifications)) {
        mod_table <- table(names(variableModifications))

        if (!all(mod_table == 1)) {
            stop(
                "Applying multiple modifications to the same amino acid is not ",
                "supported.\n",
                "Please use addVariableModifications() consecutively instead:\n\n",
                "  'ATK' |>\n",
                "    addVariableModifications(c(T = '+57.024')) |>\n",
                "    addVariableModifications(c(T = 'Phospho'))"
            )
        }
    }

    if (!length(variableModifications) || maxMods <= 0) {
        return(sequence)
    }

    terms <- .extractTerminalMods(sequence)
    sequence <- terms$sequence

    sequence <- regmatches(
        sequence,
        gregexpr("[A-Z](\\[[^\\]]*\\])*", sequence, perl = TRUE)
    )[[1]]
    seqAA <- sub("\\[.*", "", sequence)

    modifiedAA <- names(variableModifications)

    varMods <- .resolveModsToStrings(variableModifications)
    names(varMods) <- modifiedAA

    modifiable_positions_var <- which(seqAA %in% modifiedAA)
    l <- length(modifiable_positions_var)
    max_mods <- min(maxMods, l)

    .mod <- function(cmb,
                     seq_split = sequence,
                     seq_aa = seqAA,
                     var_mods = varMods) {
        mStr <- setNames(rep(NA_character_, length(seq_split)), seq_split)
        vals <- unname(var_mods[seq_aa[cmb]])
        valid <- !is.na(vals)
        mStr[cmb[valid]] <- vals[valid]
        mStr
    }

    combos <- c(
        list(setNames(rep(NA_character_, length(sequence)), sequence)),
        unlist(
            lapply(seq_len(max_mods), function(n) {
                combn(modifiable_positions_var, n, FUN = .mod, simplify = FALSE)
            }),
            recursive = FALSE
        )
    )

    inner_sequences <- unique(unname(sapply(combos, function(x) {
        paste(
            ifelse(!is.na(x),
                paste0(names(x), "[", x, "]"),
                names(x)
            ),
            collapse = ""
        )
    })))

    vapply(inner_sequences, .reattachTerminalMods,
        nterm = terms$nterm, cterm = terms$cterm,
        FUN.VALUE = character(1L), USE.NAMES = FALSE
    )
}
