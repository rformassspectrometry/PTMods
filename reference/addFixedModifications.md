# Add fixed modifications

Applies fixed modifications to peptide sequences.

The annotation style of modifications already present in \`sequences\`
is preserved as-is. The annotation style of newly applied modifications
matches the style supplied in \`fixedModifications\`: numeric values
produce deltaMass notation (e.g. \`\[+79.966331\]\`), while character
values are used verbatim (e.g. \`\[Phospho\]\` or \`\[UNIMOD:21\]\`).

## Usage

``` r
addFixedModifications(
  sequences,
  fixedModifications = c(C = 57.021464),
  pos = NULL
)
```

## Arguments

- sequences:

  Character vector. Peptide sequences that may have modifications or
  not.

- fixedModifications:

  Named \`numeric\` or \`character\`. If a \`character\` is given,
  values must be in UniMod name or UniMod ID format (e.g. \`"Phospho"\`,
  \`"UNIMOD:21"\`). The annotation style of the values is preserved in
  the output. Specifies which fixed modifications are applied to which
  amino acids. By default, carbamidomethylation is applied (+57.021464
  on cysteines). Supplying multiple modifications for the same amino
  acid is not supported; call the function consecutively instead:

        "ATK" |>
          addFixedModifications(c(T = "+57.024")) |>
          addFixedModifications(c(T = "Phospho"))

  When \`pos\` is provided, \`fixedModifications\` must have the same
  length as \`sequences\`, and names are ignored.

- pos:

  \`integer(1L)\` or \`NULL\`. When supplied, the modification is
  inserted at the specified positional residue. Names in
  \`fixedModifications\` are ignored. Must have the same length as
  \`sequences\`. Only one positional modification per sequence is
  allowed, call the function consecutively to add more. When \`pos\` is
  \`NULL\` (default), modifications are applied to all occurrences of
  the named amino acids in \`fixedModifications\`. An \`NA\` at position
  \`i\` of \`pos\` or \`fixedModifications\` causes the corresponding
  sequence to be returned unchanged.

## Value

A character vector with all fixed modifications applied, one element per
input sequence.

## Author

Guillaume Deflandre \<guillaume.deflandre@uclouvain.be\>

## Examples

``` r
addFixedModifications("ATCK")
#> [1] "ATC[+57.021464]K"
addFixedModifications(
    "ATCK",
    fixedModifications = c(Nterm = 304, C = 57.02146)
)
#> [1] "[+304]-ATC[+57.02146]K"
addFixedModifications(
    "ATSK",
    fixedModifications = c(A = -42, S = 79.966331)
)
#> [1] "A[-42]TS[+79.966331]K"
addFixedModifications(
    "ATK",
    fixedModifications = c(T = "Phospho", A = 42)
)
#> [1] "A[+42]T[Phospho]K"
addFixedModifications("A[+42]TK", fixedModifications = c(T = "Phospho"))
#> [1] "A[+42]T[Phospho]K"
addFixedModifications("AT[+79.966]AK", fixedModifications = c(A = 42))
#> [1] "A[+42]T[+79.966]A[+42]K"
addFixedModifications("ATA[+79.966]K", fixedModifications = c(A = 42))
#> [1] "A[+42]TA[+79.966][+42]K"

## Positional modification
addFixedModifications("ATK",
    fixedModifications = c(T = "Phospho"), pos = 2L
)
#> [1] "AT[Phospho]K"

## One positional modification per sequence
addFixedModifications(
    c("ATK", "PQTR"),
    fixedModifications = c("Phospho", 79.966),
    pos = c(2L, 3L)
)
#> [1] "AT[Phospho]K"  "PQT[+79.966]R"

## Apply two modifications to the same residue consecutively
"ATK" |>
    addFixedModifications(c(T = "+57.024")) |>
    addFixedModifications(c(T = "Phospho"))
#> [1] "AT[+57.024][Phospho]K"
```
