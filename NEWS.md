# PTMods 1.1

## PTMods 1.1.1

- Nothing yet.

## PTMods 1.1.0

- New Bioconductor devel.

# PTMods 0.99

## PTMods 0.99.6 

- Adjusted documentation in `convertAnnotation()`.

## PTMods 0.99.5

- Import `stats::setNames`.
- Specify package name when loading data (fixes an
  [https://github.com/lgatto/MSnbase/issues/613](issue in MSnbase)).

## PTMods 0.99.4

- Added `getModificationsCounts()`
- `convertAnnotation()` now also has a `verbose` parameter to suppress warnings
  if need be.

## PTMods 0.99.3

- Added `addModifications()`, `addFixedModifications()`,
`addVariableModifications()`, `getCanonicalSequence()` and useful PTM utils.

## PTMods 0.99.2

- Update PTMods datasets [2026-02-18].
- Improved performance of `convertAnnotation()` by calling the `modifications`
dataset before `.convertAnnotation()`.

## PTMods 0.99.1

- Bump for re-build in Bioconductor submission (due to missing
  maintainer e-mail address)

## PTMods 0.99.0

- Bump version for Bioconductor submission

## PTMods 0.3

- Update PTMods datasets [2026-01-02 Fri].
- Rename the package to `PTMods` and add vignette.
- Add `convertAnnotation()` to parse and convert between 3 PTM annotation
  styles. Currently supports `unimodId`, unimod `name` and `deltaMass`.

## PTMods 0.2

- Update datasets to latest unimod version.
- In the new unimod version some IDs are duplicated using the current
  "Name:Site:Position:NeutralLoss" ID scheme. For these we add increasing
  numbers to the end ("Name:Site:Position:NeutralLoss:DuplicateNumber").
- Now depends on R 3.5 (serialisation of RData files changed).

## PTMods 0.1

- First public version.
