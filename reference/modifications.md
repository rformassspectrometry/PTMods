# Modifications data set.

\`data.frame\` of modifications from the unimod database.

## Usage

``` r
data("modifications")
```

## Format

A \`data.frame\` with 15 columns (Id (created by
Name:(Position-)Site(:NeutralLoss) because unimod id is not unique),
UnimodId, Name, Description, AvgMass, MonoMass, Site, Position,
Classification, SpecGroup, NeutralLoss, LastModification, Approved,
Hidden) for the modifications.

## Source

Taken from the unimod database: http://www.unimod.org/xml/unimod.xml.

## Details

It was created as follows:

“\` PTMods:::.createDataSets() “\`

## Examples

``` r
data(modifications)
head(modifications)
#>                              Id UnimodId   Name Description Composition AvgMass
#> Acetyl:K               Acetyl:K        1 Acetyl Acetylation H(2) C(2) O 42.0367
#> Acetyl:N-term     Acetyl:N-term        1 Acetyl Acetylation H(2) C(2) O 42.0367
#> Acetyl:C               Acetyl:C        1 Acetyl Acetylation H(2) C(2) O 42.0367
#> Acetyl:S               Acetyl:S        1 Acetyl Acetylation H(2) C(2) O 42.0367
#> Acetyl:P-N-term Acetyl:P-N-term        1 Acetyl Acetylation H(2) C(2) O 42.0367
#> Acetyl:T               Acetyl:T        1 Acetyl Acetylation H(2) C(2) O 42.0367
#>                 MonoMass   Site       Position     Classification SpecGroup
#> Acetyl:K        42.01056      K       Anywhere           Multiple         1
#> Acetyl:N-term   42.01056 N-term     Any N-term           Multiple         2
#> Acetyl:C        42.01056      C       Anywhere Post-translational         3
#> Acetyl:S        42.01056      S       Anywhere Post-translational         4
#> Acetyl:P-N-term 42.01056 N-term Protein N-term Post-translational         5
#> Acetyl:T        42.01056      T       Anywhere Post-translational         6
#>                 NeutralLoss        LastModified Approved Hidden
#> Acetyl:K              FALSE 2017-11-08 16:08:56     TRUE  FALSE
#> Acetyl:N-term         FALSE 2017-11-08 16:08:56     TRUE  FALSE
#> Acetyl:C              FALSE 2017-11-08 16:08:56     TRUE   TRUE
#> Acetyl:S              FALSE 2017-11-08 16:08:56     TRUE   TRUE
#> Acetyl:P-N-term       FALSE 2017-11-08 16:08:56     TRUE  FALSE
#> Acetyl:T              FALSE 2017-11-08 16:08:56     TRUE   TRUE
```
