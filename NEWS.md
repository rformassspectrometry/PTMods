# PTMods 0.99

## PTMods 0.99.0

- Bump version for Bioconductor submission

## PTMods 0.3

- Update PTMods datasets [2026-01-02 Fri].
- Rename the package to `PTMods` and add vignette.
- Add `convertAnnotation()` to parse and convert between 3 PTM annotation
  styles. Currently supports `unimodId`, unimod `name` and `deltaMass`.

## unimod 0.2 [2023-02-05]

- Update datasets to latest unimod version.
- In the new unimod version some IDs are duplicated using the current
  "Name:Site:Position:NeutralLoss" ID scheme. For these we add increasing
  numbers to the end ("Name:Site:Position:NeutralLoss:DuplicateNumber").
- Now depends on R 3.5 (serialisation of RData files changed).

## unimod 0.1 [2014-04-14]

- First public version.