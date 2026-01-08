#' Internal function to create the unimod datasets.
#'
#' Used for its sideeffects
#'
#' @param path path to save the datasets
#'
#' @noRd
.createDataSets <- function(path="data") {
    xml <- .unimodDb()
    elements <- .elements(xml)
    aminoacids <- .aminoacids(xml)
    modifications <- .modifications(xml)
    attr(modifications, "Created") <-
        attr(elements, "Created") <-
        attr(aminoacids, "Created") <-
        date()
    save(elements, file=file.path(path, "elements.RData"),
         compress="xz")
    save(aminoacids, file=file.path(path, "aminoacids.RData"),
         compress="xz")
    save(modifications, file=file.path(path, "modifications.RData"),
         compress="xz")
}

#' Elements data set.
#'
#' `data.frame` of chemical elements from the unimod database.
#'
#' It was created as
#' follows:
#'
#' ```
#' PTMods:::.createDataSets()
#' ```
#'
#' @usage data("elements")
#' @format A `data.frame` with 4 columns (Name, FullName, AvgMass, MonoMass) for
#' the chemical elements.
#' @source Taken from the unimod database: http://www.unimod.org/xml/unimod.xml.
#' @examples
#' data(elements)
#' head(elements)
"elements"

#' Amminoacids data set.
#'
#' `data.frame` of aminoacids from the unimod database.
#'
#' It was created as
#' follows:
#'
#' ```
#' PTMods:::.createDataSets()
#' ```
#'
#' @usage data("aminoacids")
#' @format A `data.frame` with 11 columns (OneLetter, ThreeLetter, FullName,
#' AvgMass, MonoMass, H, C, N, O S, Se) for the aminoacids. The H/C/N/O/S/Se
#' columns contain the number of elements that build the aminoacid.
#' @source Taken from the unimod database: http://www.unimod.org/xml/unimod.xml.
#' @examples
#' data(aminoacids)
#' head(aminoacids)
"aminoacids"

#' Modifications data set.
#'
#' `data.frame` of modifications from the unimod database.
#'
#' It was created as
#' follows:
#'
#' ```
#' PTMods:::.createDataSets()
#' ```
#'
#' @usage data("modifications")
#' @format A `data.frame` with 15 columns (Id (created by
#' Name:(Position-)Site(:NeutralLoss) because unimod id is not unique),
#' UnimodId, Name, Description, AvgMass, MonoMass, Site, Position,
#' Classification, SpecGroup, NeutralLoss, LastModification, Approved,
#' Hidden) for the modifications.
#' @source Taken from the unimod database: http://www.unimod.org/xml/unimod.xml.
#' @examples
#' data(modifications)
#' head(modifications)
"modifications"