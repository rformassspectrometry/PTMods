# Elements data set.

\`data.frame\` of chemical elements from the unimod database.

## Usage

``` r
data("elements")
```

## Format

A \`data.frame\` with 4 columns (Name, FullName, AvgMass, MonoMass) for
the chemical elements.

## Source

Taken from the unimod database: http://www.unimod.org/xml/unimod.xml.

## Details

It was created as follows:

“\` PTMods:::.createDataSets() “\`

## Examples

``` r
data(elements)
head(elements)
#>     Name  FullName   AvgMass  MonoMass
#> H      H  Hydrogen  1.007940  1.007825
#> 2H    2H Deuterium  2.014102  2.014102
#> Li    Li   Lithium  6.941000  7.016003
#> C      C    Carbon 12.010700 12.000000
#> 13C  13C  Carbon13 13.003355 13.003355
#> N      N  Nitrogen 14.006700 14.003074
```
