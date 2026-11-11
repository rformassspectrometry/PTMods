# Convert sequences from one annotation style to another

Converts modifications between different annotation formats for multiple
sequences at once. See the details and examples sections for more
information. The annotation styles are inferred from the
\`modifications\` dataframe (see \`?modifications\`).

## Usage

``` r
convertAnnotation(
  x,
  convertToStyle = c("deltaMass", "unimodId", "name"),
  massTolerance = 0.01,
  verbose = TRUE
)
```

## Arguments

- x:

  Character vector with peptide sequences in ProForma format

- convertToStyle:

  Character string specifying target format. Options: "deltaMass",
  "unimodId", "name"

- massTolerance:

  Numeric mass tolerance in Daltons for matching modifications (default:
  0.01). Used when converting from deltaMass.

- verbose:

  \`logical(1)\`. If \`FALSE\`, warnings about unrecognised
  modifications are silenced (default: \`TRUE\`).

## Value

Character vector with the sequences in the target annotation format

## Details

The function handles three main conversion scenarios:

- Name to deltaMass: `"M[Oxidation]PEPTIDE"` -\>
  `"M[+15.994915]PEPTIDE"`

- Name to unimodId: `"M[Oxidation]PEPTIDE"` -\> `"M[UNIMOD:35]PEPTIDE"`

- deltaMass to name: `"M[+15.995]PEPTIDE"` -\> `"M[Oxidation]PEPTIDE"`

- deltaMass to unimodId: `"M[+15.995]PEPTIDE"` -\>
  `"M[UNIMOD:35]PEPTIDE"`

- unimodId to name: `"M[UNIMOD:35]PEPTIDE"` -\> `"M[Oxidation]PEPTIDE"`

- unimodId to deltaMass: `"M[UNIMOD:35]PEPTIDE"` -\>
  `"M[+15.994915]PEPTIDE"`

## Author

Guillaume Deflandre \<guillaume.deflandre@uclouvain.be\>

## Examples

``` r
# Convert sequence from name to delta mass
convertAnnotation("M[Oxidation]PEPTIDE", convertToStyle = "deltaMass")
#> [1] "M[+15.994915]PEPTIDE"
# Result: "M[+15.994915]PEPTIDE"

# Name to Unimod ID
convertAnnotation("M[Oxidation]PEPTIDE", convertToStyle = "unimodId")
#> [1] "M[UNIMOD:35]PEPTIDE"
# Result: "M[UNIMOD:35]PEPTIDE"

# Delta mass to name
convertAnnotation("M[+15.995]PEPTIDE", convertToStyle = "name")
#> [1] "M[Oxidation]PEPTIDE"
# Result: "M[Oxidation]PEPTIDE"

# Multiple modifications
convertAnnotation("M[Oxidation]EVNES[Phospho]PEK", convertToStyle = "deltaMass")
#> [1] "M[+15.994915]EVNES[+79.966331]PEK"
# Result: "M[+15.994915]EVNES[+79.966331]PEK"

# Convert multiple sequences from name to delta mass
sequences <- c("M[Oxidation]PEPTIDE", "EVNES[Phospho]PEK", "PEPTIDE")
convertAnnotation(sequences, convertToStyle = "deltaMass")
#> [1] "M[+15.994915]PEPTIDE" "EVNES[+79.966331]PEK" "PEPTIDE"             
# Result: c("M[+15.994915]PEPTIDE", "EVNES[+79.966331]PEK", "PEPTIDE")

# Convert from delta mass to name
sequences <- c("M[+15.995]PEPTIDE", "S[+79.966]EQUENCE")
convertAnnotation(sequences, convertToStyle = "name")
#> [1] "M[Oxidation]PEPTIDE" "S[Phospho]EQUENCE"  
# Result: c("M[Oxidation]PEPTIDE", "S[Phospho]EQUENCE")

# Convert to Unimod IDs
sequences <- c("M[Oxidation]PEPTIDE", "C[Carbamidomethyl]PEPTIDE")
convertAnnotation(sequences, convertToStyle = "unimodId")
#> [1] "M[UNIMOD:35]PEPTIDE" "C[UNIMOD:4]PEPTIDE" 
# Result: c("M[UNIMOD:35]PEPTIDE", "C[UNIMOD:4]PEPTIDE")
```
