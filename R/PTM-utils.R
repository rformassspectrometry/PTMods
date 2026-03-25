#' @title Get canonical sequence from any peptidoform
#'
#' @description Remove any modifications annotation on a sequence. Annotation
#' from any type (so long as they are specified within brackets "[..]") are
#' removed.
#'
#' @param x `character`, ProForma sequence.
#'
#' @return `character`, a `character` vector with sequences cleaned of all
#' modifications.
#'
#' @export
#'
#' @examples
#' getCanonicalSequence(
#'     c(
#'         "EM[-15.9949]EVEES[+79.9663]PEK",
#'         "EK[Acetyl]EVEES[UNIMOD:21]PEK"
#'     )
#' )
getCanonicalSequence <- function(x) {
    x <- gsub("\\[.*?\\]", "", x)
    gsub("^-|-$", "", x) ## Remove any hyphens due to terminal mods
}

#' @title Extract modification counts from one or more peptide sequences
#'
#' @description
#' Given a character vector of peptide sequences containing modifications in
#' any annotation style (deltaMass, UniMod ID, or name), returns a named
#' integer vector with the count of each unique modification. Regardless of
#' the input annotation style, all modifications are converted to their
#' UniMod name in the output (e.g. `[+79.966]` and `[UNIMOD:21]` both
#' appear as `Phospho`).
#'
#' @param sequences A `character()` vector of peptide sequences in ProForma
#' format.
#'
#' @return A named integer vector where names are UniMod modification names
#' and values are their occurrence counts across all input sequences.
#' Returns an empty integer vector if no modifications are found.
#'
#' @author Guillaume Deflandre <guillaume.deflandre@uclouvain.be>
#'
#' @examples
#' seqs <- c(
#'     "[+304]-AT[Phospho]K",
#'     "AC[Carbamidomethyl]T[+79.966]S[Phospho]K",
#'     "PEPTIDE"
#' )
#' getModificationsCounts(seqs)
#'
#' @export
getModificationsCounts <- function(sequences) {
    matches <- .getModifications(sequences)

    if (!length(matches)) {
        return(setNames(integer(0L), character(0L)))
    }

    matches <- convertAnnotation(matches, "name", verbose = FALSE)

    mods <- gsub("^-?\\[|\\]$", "", matches)
    table_counts <- table(mods)
    setNames(as.integer(table_counts), names(table_counts))
}

#' Extract raw modification strings from peptide sequences
#'
#' Collapses all input sequences into one string and extracts every
#' bracket-enclosed modification token. The strings are returned in their
#' original annotation style (deltaMass, UniMod ID, or name) without any
#' conversion. Conversion to UniMod name occurs in the calling function
#' `getModificationsCounts`.
#'
#' @param sequences A `character()` vector of peptide sequences in ProForma
#'   format.
#'
#' @return A `character()` vector of raw modification tokens including
#'   their enclosing brackets (e.g. `"[+79.966]"`, `"[Phospho]"`,
#'   `"[UNIMOD:21]"`). Returns `character(0)` if no modifications are
#'   found.
#'
#' @noRd
.getModifications <- function(sequences) {
    stopifnot(is.character(sequences))

    collapsed <- paste(sequences, collapse = "")
    regmatches(
        collapsed,
        gregexpr("\\[([^\\[\\]]+)\\]", collapsed, perl = TRUE)
    )[[1L]]
}

#' @title Parse and strip terminal modifications from a sequence
#'
#' @param sequence `character(1L)`. A peptide sequence possibly containing
#' N-terminal (e.g. `[+304]-`) and/or C-terminal (e.g. `-[+42]`) modifications.
#'
#' @return A named `list` with elements:
#' \describe{
#'   \item{`sequence`}{The inner sequence with terminal mods removed.}
#'   \item{`nterm`}{`character(1L)` N-terminal prefix (e.g. `"[+304]-"`), or `""`.}
#'   \item{`cterm`}{`character(1L)` C-terminal suffix (e.g. `"-[+42]"`), or `""`.}
#' }
#'
#' @noRd
#'
#' @examples
#' PTMods:::.extractTerminalMods("[+304]-ATK")
#' PTMods:::.extractTerminalMods("ATK-[+42]")
#' PTMods:::.extractTerminalMods("[+304]-ATK-[+42]")
#' PTMods:::.extractTerminalMods("ATK")
.extractTerminalMods <- function(sequence) {
    nterm <- ""
    cterm <- ""

    ## N-terminal: optional bracket expression followed by "-"
    nterm_match <- regmatches(sequence, regexpr("^(\\[[^\\]]*\\])-",
        sequence,
        perl = TRUE
    ))
    if (length(nterm_match) && nchar(nterm_match) > 0L) {
        nterm <- nterm_match
        sequence <- substring(sequence, nchar(nterm_match) + 1L)
    }

    ## C-terminal: "-" followed by optional bracket expression at end
    cterm_match <- regmatches(sequence, regexpr("-(\\[[^\\]]*\\])$",
        sequence,
        perl = TRUE
    ))
    if (length(cterm_match) && nchar(cterm_match) > 0L) {
        cterm <- cterm_match
        sequence <- substring(
            sequence, 1L,
            nchar(sequence) - nchar(cterm_match)
        )
    }

    list(sequence = sequence, nterm = nterm, cterm = cterm)
}

#' @title Reattach terminal modifications to a sequence
#'
#' @param sequence `character(1L)`. Inner peptide sequence without terminal mods.
#' @param nterm `character(1L)`. N-terminal prefix to prepend (e.g. `"[+304]-"`).
#' @param cterm `character(1L)`. C-terminal suffix to append (e.g. `"-[+42]"`).
#'
#' @return `character(1L)` with terminal mods reattached.
#'
#' @noRd
#'
#' @examples
#' PTMods:::.reattachTerminalMods("ATK", "[+304]-", "")
#' PTMods:::.reattachTerminalMods("ATK", "", "-[+42]")
.reattachTerminalMods <- function(sequence, nterm, cterm) {
    paste0(nterm, sequence, cterm)
}

#' @title Convert modifications to bracket-content strings, preserving style
#'
#' @description Converts a named `numeric` or `character` modifications vector
#' to a named `character` vector of bracket contents:
#' \itemize{
#'   \item `numeric` values are formatted as `"+N"` or `"-N"` (deltaMass style).
#'   \item `character` values are returned unchanged, preserving `name` or
#'   `unimodId` annotation style.
#' }
#'
#' @param modifications Named `numeric` or `character`. Modifications to
#' resolve.
#'
#' @return Named `character` vector of bracket-content strings.
#'
#' @noRd
#'
#' @examples
#' PTMods:::.resolveModsToStrings(c(A = 42.010565, T = -15.994915))
#' PTMods:::.resolveModsToStrings(c(T = "Phospho", A = "UNIMOD:1"))
.resolveModsToStrings <- function(modifications) {
    if (is.numeric(modifications)) {
        result <- vapply(modifications, function(m) {
            if (m >= 0) paste0("+", m) else as.character(m)
        }, character(1L))
        names(result) <- names(modifications)
        result
    } else {
        ## character vector: elements that look like numbers (e.g. "57.02146")
        ## get the same +/- prefix treatment as numeric values; other elements
        ## (e.g. "Phospho", "UNIMOD:21") are returned verbatim.
        result <- vapply(modifications, function(m) {
            num <- suppressWarnings(as.numeric(m))
            if (!is.na(num) && num >= 0 && !startsWith(m, "+")) {
                paste0("+", m)
            } else {
                m
            }
        }, character(1L))
        names(result) <- names(modifications)
        result
    }
}
