# Add fixed and variable modifications to peptide sequences

A convenience wrapper around \[addFixedModifications()\] and
\[addVariableModifications()\] that applies both fixed and variable
modifications in a single call. Fixed modifications are applied first,
followed by the enumeration of all possible combinations of variable
modifications. Optionally, the annotation style of the output can be
converted via \[convertAnnotation()\].

## Usage

``` r
addModifications(
  sequences,
  fixedModifications = NULL,
  variableModifications = NULL,
  convertToStyle = NULL,
  pos = NULL,
  ...
)
```

## Arguments

- sequences:

  \`character\`. Peptide sequences with or without existing modification
  annotations. Modifications must be enclosed in square brackets (e.g.
  \`"AT\[+79.966\]K"\`, \`"\[+304\]-ATK"\`).

- fixedModifications:

  Named \`numeric\` or \`character\`, or \`NULL\`. Fixed modifications
  to be applied unconditionally to the specified amino acids or termini.
  Names correspond to single-letter amino acid codes, \`"Nterm"\`, or
  \`"Cterm"\`. If \`character\`, values must be in UniMod name or UniMod
  ID format (e.g. \`"Phospho"\`, \`"UNIMOD:21"\`). See
  \[addFixedModifications()\] for full details. Set to \`NULL\` to skip.

- variableModifications:

  Named \`numeric\` or \`character\`, or \`NULL\`. Variable
  modifications to enumerate over the specified amino acids. Names
  correspond to single-letter amino acid codes. If \`character\`, values
  must be in UniMod name or UniMod ID format. See
  \[addVariableModifications()\] for full details. Set to \`NULL\` to
  skip.

- convertToStyle:

  \`character(1)\`. Controls the annotation format of modifications in
  the returned sequences. One of:

  \`NULL\`

  : (default) Same as input modification style.

  \`"deltaMass"\`

  : Delta mass notation, e.g. \`\[+79.966\]\`.

  \`"unimodId"\`

  : UniMod ID notation, e.g. \`\[UNIMOD:21\]\`.

  \`"name"\`

  : UniMod name notation, e.g. \`\[Phospho\]\`.

  Conversion is performed by \[convertAnnotation()\].

- pos:

  \`integer\` or \`NULL\`. When supplied, fixed modifications are
  applied at these 1-based canonical positions rather than by amino-acid
  name. Must have the same length as \`sequences\`. Forwarded to
  \[addFixedModifications()\]. Ignored when \`fixedModifications\` is
  \`NULL\`.

- ...:

  Additional arguments passed to the underlying modification functions.
  Currently supported:

  \`maxMods\`

  : Forwarded to \[addVariableModifications()\]. A \`numeric(1)\`
    specifying the maximum number of variable modifications that may be
    applied simultaneously. Defaults to \`Inf\`.

## Value

A \`character\` vector of peptidoforms. The length is greater than or
equal to \`length(sequences)\`, as each input sequence may expand into
multiple peptidoforms corresponding to all combinations of variable
modifications.

## See also

\*
[`addFixedModifications`](https://rformassspectrometry.github.io/PTMods/reference/addFixedModifications.md)
for applying fixed modifications.

\*
[`addVariableModifications`](https://rformassspectrometry.github.io/PTMods/reference/addVariableModifications.md)
for enumerating variable modification combinations.

\*
[`convertAnnotation`](https://rformassspectrometry.github.io/PTMods/reference/convertAnnotation.md)
for converting between modification annotation styles.

## Author

Guillaume Deflandre \<guillaume.deflandre@uclouvain.be\>

## Examples

``` r
## Fixed modifications only
addModifications(
    c("ATCK", "PEPCK"),
    fixedModifications = c(C = 57.02146)
)
#> [1] "ATC[+57.02146]K"  "PEPC[+57.02146]K"

## Variable modifications only
addModifications(
    c("ATSK", "PEPTSK"),
    variableModifications = c(S = 79.966331, T = 79.966331),
    maxMods = 2
)
#> [1] "ATSK"                           "AT[+79.966331]SK"              
#> [3] "ATS[+79.966331]K"               "AT[+79.966331]S[+79.966331]K"  
#> [5] "PEPTSK"                         "PEPT[+79.966331]SK"            
#> [7] "PEPTS[+79.966331]K"             "PEPT[+79.966331]S[+79.966331]K"

## Both fixed and variable modifications
addModifications(
    "ATCSK",
    fixedModifications = c(C = 57.02146),
    variableModifications = c(S = 79.966331, T = 79.966331),
    maxMods = 1
)
#> [1] "ATC[+57.02146]SK"             "AT[+79.966331]C[+57.02146]SK"
#> [3] "ATC[+57.02146]S[+79.966331]K"

## Return annotations in UniMod name style
addModifications(
    "ATSK",
    variableModifications = c(S = 79.966331, T = 79.966331),
    convertToStyle = "name"
)
#> [1] "ATSK"                   "AT[Phospho]SK"          "ATS[Phospho]K"         
#> [4] "AT[Phospho]S[Phospho]K"

## N-terminal fixed modification with variable modifications
addModifications(
    "ATSK",
    fixedModifications = c(Nterm = 304),
    variableModifications = c(S = 79.966331)
)
#> [1] "[+304]-ATSK"             "[+304]-ATS[+79.966331]K"

## Positional fixed modification
addModifications(
    c("ATK", "PQTR"),
    fixedModifications = c("Phospho", 79.966),
    pos = c(2L, 3L)
)
#> [1] "AT[Phospho]K"  "PQT[+79.966]R"
```
