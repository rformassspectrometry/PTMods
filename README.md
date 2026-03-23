# PTMods: managing post-translational modifications in R

[![license](http://img.shields.io/badge/license-GPL%20%28%3E=%203%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)

A package to managing post-translational modifications in `R`.

The [reference manual](https://rformassspectrometry.github.io/PTMods/articles/PTMods.html) 
is a good way to get started with the package.

Do not hesitate to contact us through issues or pull requests !
We are open to suggestions on how to improve the package at all times. 

## Suggested packages

We highly suggest combining `PTMods` with 
[PSMatch](https://github.com/rformassspectrometry/PSMatch) and 
[Spectra](https://github.com/rformassspectrometry/Spectra) for a more
complete approach to analysing mass spectrometry data in R. 

## Installation instructions

To install the package from Bioconductor, make sure you have the
`BiocManager` package, available from CRAN, and then run

```r
BiocManager::install("PTMods")
```

## Credit

Code in this package is focused on the use of UniMod data and the 
[ProForma](https://github.com/HUPO-PSI/ProForma) standard Proteoform and 
Peptidoform Notation formats with the goal to integrate them within the 
[*R for Mass Spectrometry*](https://github.com/rformassspectrometry)
infrastructure.
