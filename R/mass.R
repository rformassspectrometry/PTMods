utils::globalVariables("aminoacids")


#' Calculate mass for a given aminoacid sequence
#'
#' @param x character, sequence
#' @param type character, mono or average mass
#' @return list of numeric (mass values)
#' @noRd
.aamass <- function(x, type=c("MonoMass", "AvgMass")) {
    m <- aminoacids[, match.arg(type)]
    names(m) <- aminoacids$OneLetter
    vapply(.string2character(x), function(xx)sum(m[xx]), NA_real_)
}

#' Calculate a single modification for a given sequence
#'
#' Throws a message if the default rule is applied and not verified (added to
#' the switch statement) by a human
#'
#' @param x character, sequence
#' @param id character, modification id
#' @param type character, mono or average mass
#' @param msg logical, show message if default was used?
#' @return numeric (mass values)
#' @noRd
.unimodMass <- function(x, id, type=c("MonoMass", "AvgMass"), msg=TRUE) {
    stopifnot(is.character(x) && is.character(id) && length(id) == 1L)
    stopifnot(id %in% modifications$Id)
    type <- match.arg(type)
    switch(id,
        ## 765: Met-loss
        "Met-loss:P-M" =
            # N-terminal initiator methionine is removed by a methionine
            # aminopeptidase from proteins where the residue following the
            # methionine is Ala, Cys, Gly, Pro, Ser, Thr or Val. This is
            # generally the final N-terminal state for proteins where the
            # following residue was a Cys, Pro or Val.
            grepl("^M([ACGPSTV])", x) * modifications[id, type],
        ## 766: Met-loss+Acetyl
        "Met-loss+Acetyl:P-M" =
            # The N-terminal initiator methionine is removed by a methionine
            # aminopeptidase from proteins whose residue following the
            # methionine is Ala, Cys, Gly, Pro, Ser, Thr or Val. Proteins
            # whose following residue was Ala, Gly, Ser or Thr are then
            # acetylated by an N(alpha)-acetyltransferase on the new
            # N-terminus.
            grepl("^M([AGST])", x) * modifications[id, type] +
            grepl("^M([CPV])", x) * modifications["Met-loss:P-M", type],
        ## default verified:
        ## 1: Acetyl
        "Acetyl:K" =, "Acetyl:N-term" =, "Acetyl:C" =, "Acetyl:S" =,
        "Acetyl:P-N-term" =, "Acetyl:T" =, "Acetyl:Y" =, "Acetyl:H" =,
        "Acetyl:R" =,
        ## 4: Carbamidomethyl
        "Carbamidomethyl:C" =, "Carbamidomethyl:K" =,
        "Carbamidomethyl:N-term" =, "Carbamidomethyl:H" =,
        "Carbamidomethyl:D" =, "Carbamidomethyl:E" =, "Carbamidomethyl:S" =,
        "Carbamidomethyl:T" =, "Carbamidomethyl:Y" =, "Carbamidomethyl:U" =,
        "Carbamidomethyl:M" =, "Carbamidomethyl:M:NL" =,
        ## 21: Phospho
        "Phospho:T" =, "Phospho:S" =, "Phospho:Y" =, "Phospho:D" =,
        "Phospho:H" =, "Phospho:C" =, "Phospho:R" =, "Phospho:K" =,
        "Phospho:T:NL" =, "Phospho:S:NL" =
        {
            .countSite(x, modifications[id, "Site"]) * modifications[id, type]
        },
        {
            if (msg) {
                message(
                    "Applying the default rule for the modification: ", id,
                    "\nPlease create an issue on: https://github.com/",
                    "ComputationalProteomicsUnit/unimod/issues/new ",
                    "\nto let us implement the correct rule or if the default ",
                    "one is already correct we could remove this message."
                )
            }
            .countSite(x, modifications[id, "Site"]) * modifications[id, type]
        }
    )
}

#' Modify the sequence for a single modification
#'
#' @param x character, sequence
#' @param id character, modification id
#' @return character (sequence)
#' @noRd
.unimodSequence <- function(x, id) {
    stopifnot(is.character(x) && is.character(id) && length(id) == 1L)
    stopifnot(id %in% modifications$Id)
    switch(id,
        ## 765: Met-loss
        "Met-loss:P-M" =,
        ## 766: Met-loss+Acetyl
        "Met-loss+Acetyl:P-M" =
            gsub("^M([ACGPSTV])", "\\1", x),
        x
    )
}

#' Find site occurrence
#'
#' @param x character, sequence
#' @param site character
#' @return numeric
#' @noRd
.countSite <- function(x, site) {
    stopifnot(is.character(x) && is.character(site) && length(site) == 1L)
    if (endsWith(site, "term")) {
        rep.int(1L, length(x))
    } else {
        nchar(x) - nchar(gsub(site, "", x, fixed=TRUE))
        ## while looking a little bit odd it is much faster than:
        # vapply(
        #     gregexpr(pattern=site, text=x, fixed=fixed),
        #     function(rx)sum(rx > 0L), double(1L),
        #     USE.NAMES=FALSE
        # )
    }
}

#' Calculate mass for a given aminoacid sequence and apply modifications
#'
#' @param x character, sequence
#' @param type character, mono or average mass
#' @param fixedModifications character/data.frame
#' @param variableModifications data.frame
#' @return list of numeric (mass values)
#' @noRd
.mass <- function(x, type=c("MonoMass", "AvgMass"),
                  fixedModifications=NULL,
                  variableModifications=NULL) {
    type <- match.arg(type)

    if (is.character(fixedModifications)) {
        i <- match(fixedModifications, modifications$Id)
        if (anyNA(i)) {
            stop(
                paste0(fixedModifications[is.na(i)], collapse=", "),
                ifelse(sum(is.na(i)) == 1L, " is ", " are "),
                "not part of the PTMods::modifications data set!"
            )
        }
        isUnimod <- TRUE
        fixedModifications <- modifications[i,,drop=FALSE]
    } else if (is.data.frame(fixedModifications)) {
        if (!all(c("Id", type, "Site") %in% colnames(fixedModifications))) {
            stop("The 'fixedModifications' `data.frame` has to have at least ",
                 "the columns 'Id', '", type, "', and 'Site'!")
        }
        isUnimod <- all(fixedModifications$Id %in% modifications$Id)
    } else if (!is.null(fixedModifications)) {
        stop("'fixedModifications' must be a `character`, `data.frame`, or ",
             "`NULL`!")
    }

    if (!is.null(variableModifications)) {
        stop("Variable Modifications are not supported yet!")
    }

    m <- .aamass(x, type)

    if (anyDuplicated(fixedModifications$Site)) {
        stop("Duplicated fixed modification sites are not allowed!")
    }

    if (!is.null(fixedModifications)) {
        for (i in seq_len(nrow(fixedModifications))) {
            if (isUnimod) {
                m <- m + .unimodMass(x, fixedModifications$Id[i], type=type)
                x <- .unimodSequence(x, fixedModifications$Id[i])
            } else {
                m <- m + .countSite(x, fixedModifications$Site[i]) *
                    fixedModifications[i, type]
            }
        }
    }
    attr(m, "sequence") <- x
    m
}
