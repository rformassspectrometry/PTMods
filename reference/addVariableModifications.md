# Add variable modifications

Generates all possible combinations of variable modifications for
peptide sequences.

The annotation style of modifications already present in \`sequences\`
is preserved as-is. The annotation style of newly applied modifications
matches the style supplied in \`variableModifications\`: numeric values
produce deltaMass notation (e.g. \`\[+79.966331\]\`), while character
values are used verbatim (e.g. \`\[Phospho\]\` or \`\[UNIMOD:21\]\`).

## Usage

``` r
addVariableModifications(
  sequences,
  variableModifications = NULL,
  maxMods = Inf
)
```

## Arguments

- sequences:

  Character vector. Peptide sequences that may have modifications or
  not.

- variableModifications:

  Named \`numeric\` or \`character\`. If a \`character\` is given,
  values must be in UniMod name or UniMod ID format (e.g. \`"Phospho"\`,
  \`"UNIMOD:21"\`). The annotation style of the values is preserved in
  the output. Specifies which variable modifications are used on which
  amino acids. Supplying multiple modifications for the same amino acid
  is not supported; call the function consecutively instead:

        seq |>
          addVariableModifications(c(T = "+57.024")) |>
          addVariableModifications(c(T = "Phospho"))

- maxMods:

  Numeric. Indicates how many modifications can be applied at once.

## Value

A named list of character vectors with all possible combinations of
modifications, one element per input sequence.

## Author

Guillaume Deflandre \<guillaume.deflandre@uclouvain.be\>

## Examples

``` r
addVariableModifications(
    "APGHKA",
    variableModifications = c(A = 4, K = 5, S = 8),
    maxMods = 2
)
#> [1] "APGHKA"         "A[+4]PGHKA"     "APGHK[+5]A"     "APGHKA[+4]"    
#> [5] "A[+4]PGHK[+5]A" "A[+4]PGHKA[+4]" "APGHK[+5]A[+4]"
addVariableModifications(
    "APGHKA",
    variableModifications = c(A = 4, K = 5, S = 8),
    maxMods = 3
)
#> [1] "APGHKA"             "A[+4]PGHKA"         "APGHK[+5]A"        
#> [4] "APGHKA[+4]"         "A[+4]PGHK[+5]A"     "A[+4]PGHKA[+4]"    
#> [7] "APGHK[+5]A[+4]"     "A[+4]PGHK[+5]A[+4]"
addVariableModifications(
    c("ATK", "PAK"),
    variableModifications = c(T = "Phospho", A = "UNIMOD:1")
)
#> [1] "ATK"                    "A[UNIMOD:1]TK"          "AT[Phospho]K"          
#> [4] "A[UNIMOD:1]T[Phospho]K" "PAK"                    "PA[UNIMOD:1]K"         
addVariableModifications(
    "A[+10]KT[-5]",
    variableModifications = c(A = -20, T = "Phospho")
)
#> [1] "A[+10]KT[-5]"               "A[+10][-20]KT[-5]"         
#> [3] "A[+10]KT[-5][Phospho]"      "A[+10][-20]KT[-5][Phospho]"

## Apply two modifications to the same residue consecutively
"ATK" |>
    addVariableModifications(c(T = "+57.024")) |>
    addVariableModifications(c(T = "Phospho"))
#> [1] "ATK"                   "AT[Phospho]K"          "AT[+57.024]K"         
#> [4] "AT[+57.024][Phospho]K"
```
