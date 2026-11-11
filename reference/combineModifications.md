# Combines consecutive modifications on amino acids in peptide sequences.

Combines consecutive modifications on amino acids in peptide sequences.

## Usage

``` r
combineModifications(sequences)
```

## Arguments

- sequences:

  \`character\`. A character vector of peptide sequences with
  modifications in any format accepted by \`convertAnnotation\`.
  Multiple modifications per amino acid are allowed.

## Value

A character vector of the same length as \`sequences\` with consecutive
modifications on each amino acid summed into a single modification. This
summed modifications is written in delta mass format.

## Author

Guillaume Deflandre \<guillaume.deflandre@uclouvain.be\>

## Examples

``` r
combineModifications("AT[+10][-5]K")
#> [1] "AT[+5]K"
combineModifications("[+10]-[-5]-ATK")
#> [1] "[+5]-ATK"
combineModifications("AT[Phospho][UNIMOD:1]K")
#> [1] "AT[+121.976896]K"
combineModifications("EM[+15.9949][-10]EK[+42][-8]")
#> [1] "EM[+5.9949]EK[+34]"
combineModifications(c("AT[+10][-5]K", "EM[+15.9949][-10]EK"))
#> [1] "AT[+5]K"       "EM[+5.9949]EK"
```
