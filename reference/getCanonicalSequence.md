# Get canonical sequence from any peptidoform

Remove any modifications annotation on a sequence. Annotation from any
type (so long as they are specified within brackets "\[..\]") are
removed.

## Usage

``` r
getCanonicalSequence(x)
```

## Arguments

- x:

  \`character\`, ProForma sequence.

## Value

\`character\`, a \`character\` vector with sequences cleaned of all
modifications.

## Examples

``` r
getCanonicalSequence(
    c(
        "EM[-15.9949]EVEES[+79.9663]PEK",
        "EK[Acetyl]EVEES[UNIMOD:21]PEK"
    )
)
#> [1] "EMEVEESPEK" "EKEVEESPEK"
```
