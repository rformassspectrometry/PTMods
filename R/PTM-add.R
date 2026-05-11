#' @title Add fixed and variable modifications to peptide sequences
#'
#' @description
#' A convenience wrapper around [addFixedModifications()] and
#' [addVariableModifications()] that applies both fixed and variable
#' modifications in a single call. Fixed modifications are applied first,
#' followed by the enumeration of all possible combinations of variable
#' modifications. Optionally, the annotation style of
#' the output can be converted via [convertAnnotation()].
#'
#' @param sequences `character`. Peptide sequences with or without existing
#'   modification annotations. Modifications must be enclosed in square
#'   brackets (e.g. `"AT[+79.966]K"`, `"[+304]-ATK"`).
#'
#' @param fixedModifications Named `numeric` or `character`, or `NULL`.
#'   Fixed modifications to be applied unconditionally to the specified amino
#'   acids or termini. Names correspond to single-letter amino acid codes,
#'   `"Nterm"`, or `"Cterm"`. If `character`, values must be in UniMod name
#'   or UniMod ID format (e.g. `"Phospho"`, `"UNIMOD:21"`). See
#'   [addFixedModifications()] for full details. Set to `NULL` to skip.
#'
#' @param pos `integer` or `NULL`. When supplied, fixed modifications are
#'   applied at these 1-based canonical positions rather than by amino-acid
#'   name. Must have the same length as `sequences`. Forwarded to
#'   [addFixedModifications()]. Ignored when `fixedModifications` is `NULL`.
#'
#' @param variableModifications Named `numeric` or `character`, or `NULL`.
#'   Variable modifications to enumerate over the specified amino acids.
#'   Names correspond to single-letter amino acid codes. If `character`,
#'   values must be in UniMod name or UniMod ID format. See
#'   [addVariableModifications()] for full details. Set to `NULL` to skip.
#'
#' @param convertToStyle `character(1)`. Controls the annotation format of
#'   modifications in the returned sequences. One of:
#'   \describe{
#'     \item{`NULL`}{(default) Same as input modification style.}
#'     \item{`"deltaMass"`}{Delta mass notation, e.g. `[+79.966]`.}
#'     \item{`"unimodId"`}{UniMod ID notation, e.g. `[UNIMOD:21]`.}
#'     \item{`"name"`}{UniMod name notation, e.g. `[Phospho]`.}
#'   }
#'   Conversion is performed by [convertAnnotation()].
#'
#' @param ... Additional arguments passed to the underlying modification
#'   functions. Currently supported:
#'   \describe{
#'     \item{`maxMods`}{Forwarded to [addVariableModifications()]. A
#'       `numeric(1)` specifying the maximum number of variable modifications
#'       that may be applied simultaneously. Defaults to `Inf`.}
#'   }
#'
#' @return A `character` vector of peptidoforms. The length is greater than or
#'   equal to `length(sequences)`, as each input sequence may expand into
#'   multiple peptidoforms corresponding to all combinations of variable
#'   modifications.
#'
#' @seealso
#' * \code{\link{addFixedModifications}} for applying fixed modifications.
#'
#' * \code{\link{addVariableModifications}} for enumerating variable modification
#'   combinations.
#'
#' * \code{\link{convertAnnotation}} for converting between modification annotation
#'   styles.
#'
#' @author Guillaume Deflandre <guillaume.deflandre@uclouvain.be>
#'
#' @export
#'
#' @examples
#' ## Fixed modifications only
#' addModifications(
#'     c("ATCK", "PEPCK"),
#'     fixedModifications = c(C = 57.02146)
#' )
#'
#' ## Variable modifications only
#' addModifications(
#'     c("ATSK", "PEPTSK"),
#'     variableModifications = c(S = 79.966331, T = 79.966331),
#'     maxMods = 2
#' )
#'
#' ## Both fixed and variable modifications
#' addModifications(
#'     "ATCSK",
#'     fixedModifications = c(C = 57.02146),
#'     variableModifications = c(S = 79.966331, T = 79.966331),
#'     maxMods = 1
#' )
#'
#' ## Return annotations in UniMod name style
#' addModifications(
#'     "ATSK",
#'     variableModifications = c(S = 79.966331, T = 79.966331),
#'     convertToStyle = "name"
#' )
#'
#' ## N-terminal fixed modification with variable modifications
#' addModifications(
#'     "ATSK",
#'     fixedModifications = c(Nterm = 304),
#'     variableModifications = c(S = 79.966331)
#' )
#'
#' ## Positional fixed modification
#' addModifications(
#'     c("ATK", "PQTR"),
#'     fixedModifications = c("Phospho", 79.966),
#'     pos = c(2L, 3L)
#' )
addModifications <- function(sequences,
                             fixedModifications = NULL,
                             variableModifications = NULL,
                             convertToStyle = NULL,
                             pos = NULL,
                             ...) {
    dots <- list(...)

    ## --- Fixed modifications --- ##
    if (!is.null(fixedModifications)) {
        sequences <- addFixedModifications(
            sequences,
            fixedModifications = fixedModifications,
            pos = pos
        )
    }

    ## --- Variable modifications ---##
    if (!is.null(variableModifications)) {
        maxMods <- if (!is.null(dots[["maxMods"]])) dots[["maxMods"]] else Inf
        sequences <- addVariableModifications(
            sequences,
            variableModifications = variableModifications,
            maxMods = maxMods
        )
    }

    ## --- Convert annotation style --- ##
    if (!is.null(convertToStyle)) {
        sequences <- convertAnnotation(sequences,
            convertToStyle = convertToStyle
        )
    }

    sequences
}
