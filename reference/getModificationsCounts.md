# Extract modification counts from one or more peptide sequences

Given a character vector of peptide sequences containing modifications
in any annotation style (deltaMass, UniMod ID, or name), returns a named
integer vector with the count of each unique modification. Regardless of
the input annotation style, all modifications are converted to their
UniMod name in the output (e.g. \`\[+79.966\]\` and \`\[UNIMOD:21\]\`
both appear as \`Phospho\`).

## Usage

``` r
getModificationsCounts(sequences)
```

## Arguments

- sequences:

  A \`character()\` vector of peptide sequences in ProForma format.

## Value

A named integer vector where names are UniMod modification names and
values are their occurrence counts across all input sequences. Returns
an empty integer vector if no modifications are found.

## Author

Guillaume Deflandre \<guillaume.deflandre@uclouvain.be\>

## Examples

``` r
seqs <- c(
    "[+304]-AT[Phospho]K",
    "AC[Carbamidomethyl]T[+79.966]S[Phospho]K",
    "PEPTIDE"
)
getModificationsCounts(seqs)
#>            +304 Carbamidomethyl         Phospho 
#>               1               1               3 
```
