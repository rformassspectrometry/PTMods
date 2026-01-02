# unimod 0.x

## unimod 0.3 

- Add `convertAnnotation()` to parse and convert between 3 PTM annotation 
  styles. Currently only supports `unimodId`, unimod `name` and `deltaMass`.

## unimod 0.2 [2023-02-05]

- Update datasets to latest unimod version.
- In the new unimod version some IDs are duplicated using the current
  "Name:Site:Position:NeutralLoss" ID scheme. For these we add increasing
  numbers to the end ("Name:Site:Position:NeutralLoss:DuplicateNumber").
- Now depends on R 3.5 (serialisation of RData files changed).

## unimod 0.1 [2014-04-14]

- First public version.
