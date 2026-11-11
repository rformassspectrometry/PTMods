# Amminoacids data set.

\`data.frame\` of aminoacids from the unimod database.

## Usage

``` r
data("aminoacids")
```

## Format

A \`data.frame\` with 11 columns (OneLetter, ThreeLetter, FullName,
AvgMass, MonoMass, H, C, N, O S, Se) for the aminoacids. The
H/C/N/O/S/Se columns contain the number of elements that build the
aminoacid.

## Source

Taken from the unimod database: http://www.unimod.org/xml/unimod.xml.

## Details

It was created as follows:

“\` PTMods:::.createDataSets() “\`

## Examples

``` r
data(aminoacids)
head(aminoacids)
#>   OneLetter ThreeLetter      FullName  AvgMass  MonoMass  H C N O S Se
#> -         -                             0.0000   0.00000  0 0 0 0 0  0
#> A         A         Ala       Alanine  71.0779  71.03711  5 3 1 1 0  0
#> R         R         Arg      Arginine 156.1857 156.10111 12 6 4 1 0  0
#> N         N         Asn    Asparagine 114.1026 114.04293  6 4 2 2 0  0
#> D         D         Asp Aspartic acid 115.0874 115.02694  5 4 1 3 0  0
#> C         C         Cys      Cysteine 103.1429 103.00919  5 3 1 1 1  0
```
