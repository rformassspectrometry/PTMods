# The PTMods package: a package to handle post-translational modifications

**Package**: *[PTMods](https://bioconductor.org/packages/3.24/PTMods)*\
**Authors**: Laurent Gatto \[aut\] (ORCID:
<https://orcid.org/0000-0002-1520-2268>), Sebastian Gibb \[aut\] (ORCID:
<https://orcid.org/0000-0001-7406-4443>), Guillaume Deflandre \[cre,
aut\] (ORCID: <https://orcid.org/0009-0008-1257-2416>)\
**Last modified:** 2026-11-12 13:14:52.900579\
**Compiled**: Thu Nov 12 13:16:00 2026

## Introduction

Post-translational modifications (PTMs) are covalent modifications that
happen after protein biosynthesis and that can have a major influence on
protein function. PTMs are biologically very important and can be
identified in a high-throughput way using mass spectrometry.

The `PTMods` package focuses on handling such PTMs and works in
conjuction with the
*[PSMatch](https://bioconductor.org/packages/3.24/PSMatch)* and
*[Spectra](https://bioconductor.org/packages/3.24/Spectra)* packages. It
distributes PTMs from the Unimod database and allows to convert PTM
annotations between multiple annotation formats. The package also allows
the seemless addition of fixed and/or variable PTMs using different
annotation formats.

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`PTMods`](https://github.com/RforMassSpectrometry/PTMods)`)`

### Installation

The package can be installed from Bioconductor with

\
`if`` ``(``!`[`require`](https://rdrr.io/r/base/library.html)`(`[`"BiocManager"`](https://bioconductor.github.io/BiocManager/)`, quietly ``=`` ``TRUE``)``)`` ``{`\
`    `[`install.packages`](https://rdrr.io/r/utils/install.packages.html)`(``"BiocManager"``)`\
`}`\
\
`BiocManager``::`[`install`](https://bioconductor.github.io/BiocManager/reference/install.html)`(``"PTMods"``)`

Alternatively, altough this isn’t recommended in general, it can also be
installed from Github.

\
`BiocManager``::`[`install`](https://bioconductor.github.io/BiocManager/reference/install.html)`(``"RforMassSpectrometry/PTMods"``)`

## The Unimod database

The `PTMods` package distributes the data from the Unimod database:

> Unimod is a public domain database, distributed under a [copyleft
> licence](http://www.gnu.org/licenses/license-list.md): “a copyright
> notice that permits unrestricted redistribution and modification,
> provided that all copies and derivatives retain the same permissions.”
>
> The aim is to create a community supported, comprehensive database of
> protein modifications for mass spectrometry applications. That is,
> accurate and verifiable values, derived from elemental compositions,
> for the mass differences introduced by all types of natural and
> artificial modifications. Other important information includes any
> mass change, (neutral loss), that occurs during MS/MS analysis, and
> site specificity, (which residues are susceptible to modification and
> any constraints on the position of the modification within the protein
> or peptide).

Source:
[http://www.unimod.org/unimod_help.html](http://www.unimod.org/unimod_help.md)
We try to keep the data up to date with the Unimod database. If you
believe it isn’t, do not hesitate to contact us through an [issue on the
GitHub page](https://github.com/rformassspectrometry/PTMods/issues/new).

The package provides the following objects.

- The `modifications` dataframe describing 3704 PTMs from the unimod
  database:

\
[`data`](https://rdrr.io/r/utils/data.html)`(``modifications``)`\
[`head`](https://rdrr.io/r/utils/head.html)`(``modifications``)`

    ##                              Id UnimodId   Name Description Composition AvgMass
    ## Acetyl:K               Acetyl:K        1 Acetyl Acetylation H(2) C(2) O 42.0367
    ## Acetyl:N-term     Acetyl:N-term        1 Acetyl Acetylation H(2) C(2) O 42.0367
    ## Acetyl:C               Acetyl:C        1 Acetyl Acetylation H(2) C(2) O 42.0367
    ## Acetyl:S               Acetyl:S        1 Acetyl Acetylation H(2) C(2) O 42.0367
    ## Acetyl:P-N-term Acetyl:P-N-term        1 Acetyl Acetylation H(2) C(2) O 42.0367
    ## Acetyl:T               Acetyl:T        1 Acetyl Acetylation H(2) C(2) O 42.0367
    ##                 MonoMass   Site       Position     Classification SpecGroup
    ## Acetyl:K        42.01056      K       Anywhere           Multiple         1
    ## Acetyl:N-term   42.01056 N-term     Any N-term           Multiple         2
    ## Acetyl:C        42.01056      C       Anywhere Post-translational         3
    ## Acetyl:S        42.01056      S       Anywhere Post-translational         4
    ## Acetyl:P-N-term 42.01056 N-term Protein N-term Post-translational         5
    ## Acetyl:T        42.01056      T       Anywhere Post-translational         6
    ##                 NeutralLoss        LastModified Approved Hidden
    ## Acetyl:K              FALSE 2017-11-08 16:08:56     TRUE  FALSE
    ## Acetyl:N-term         FALSE 2017-11-08 16:08:56     TRUE  FALSE
    ## Acetyl:C              FALSE 2017-11-08 16:08:56     TRUE   TRUE
    ## Acetyl:S              FALSE 2017-11-08 16:08:56     TRUE   TRUE
    ## Acetyl:P-N-term       FALSE 2017-11-08 16:08:56     TRUE  FALSE
    ## Acetyl:T              FALSE 2017-11-08 16:08:56     TRUE   TRUE

- The `aminoacids` dataframe describing 24 amino acids from the unimod
  database:

\
[`data`](https://rdrr.io/r/utils/data.html)`(``aminoacids``)`\
[`head`](https://rdrr.io/r/utils/head.html)`(``aminoacids``)`

    ##   OneLetter ThreeLetter      FullName  AvgMass  MonoMass  H C N O S Se
    ## -         -                             0.0000   0.00000  0 0 0 0 0  0
    ## A         A         Ala       Alanine  71.0779  71.03711  5 3 1 1 0  0
    ## R         R         Arg      Arginine 156.1857 156.10111 12 6 4 1 0  0
    ## N         N         Asn    Asparagine 114.1026 114.04293  6 4 2 2 0  0
    ## D         D         Asp Aspartic acid 115.0874 115.02694  5 4 1 3 0  0
    ## C         C         Cys      Cysteine 103.1429 103.00919  5 3 1 1 1  0

- The `elements` data.frame describing 40 chemical elements from the
  unimod database:

\
[`data`](https://rdrr.io/r/utils/data.html)`(``elements``)`\
[`head`](https://rdrr.io/r/utils/head.html)`(``elements``)`

    ##     Name  FullName   AvgMass  MonoMass
    ## H      H  Hydrogen  1.007940  1.007825
    ## 2H    2H Deuterium  2.014102  2.014102
    ## Li    Li   Lithium  6.941000  7.016003
    ## C      C    Carbon 12.010700 12.000000
    ## 13C  13C  Carbon13 13.003355 13.003355
    ## N      N  Nitrogen 14.006700 14.003074

Note that the goal of `PTMods` is not to place this information in a
biological context. There are other sources of biological information
about post-translational modifications:

- [RESID](https://proteininformationresource.org/resid/resid.shtml) for
  post-translational modifications
- [Uniprot/Swiss-Prot](https://www.uniprot.org/) protein sequence
  database
- [Prosite](https://prosite.expasy.org/) database of protein families
  and domains
- [Glycan Database](https://www.functionalglycomics.org/) from the
  Consortium for Functional Glycomics
- [PhosphoDB](https://phosphodb.hecklab.com/site/index) and
  [PhosphoSitePlus](https://www.phosphosite.org/homeAction.action) for
  phosphoryations
- …

Other compilations of modifications with a focus on mass spectrometry:

- [Delta Mass](https://abrf.org/resources/delta-mass/)
- [FindMod](https://web.expasy.org/findmod/)
- …

## Different types of PTM annotations

PTMs can be annotated using different syntaxes:

- The *deltaMass* syntax that describes the mass shift of an amino acid
  such as for example `"M[+15.994915]PEPTIDE"` for an oxidation.
- The *unimodId* that uses the unimod identifier such as for example
  `"M[UNIMOD:35]PEPTIDE"` for an oxidation.
- The modification’s name, such as for example `"M[Oxidation]PEPTIDE"`.

The `convertAnnotion()` function can be used to convert between these
different syntaxes:

\
[`convertAnnotation`](https://rformassspectrometry.github.io/PTMods/reference/convertAnnotation.md)`(``"M[Oxidation]PEPTIDE"``, convertToStyle ``=`` ``"deltaMass"``)`

    ## [1] "M[+15.994915]PEPTIDE"

\
[`convertAnnotation`](https://rformassspectrometry.github.io/PTMods/reference/convertAnnotation.md)`(``"M[Oxidation]PEPTIDE"``, convertToStyle ``=`` ``"unimodId"``)`

    ## [1] "M[UNIMOD:35]PEPTIDE"

\
[`convertAnnotation`](https://rformassspectrometry.github.io/PTMods/reference/convertAnnotation.md)`(``"M[+15.995]PEPTIDE"``, convertToStyle ``=`` ``"name"``)`

    ## [1] "M[Oxidation]PEPTIDE"

## Add modifications to a sequence

PTMs may be added (in any of the formats mentioned above) using the
[`addFixedModifications()`](https://rformassspectrometry.github.io/PTMods/reference/addFixedModifications.md)
and
[`addVariableModifications()`](https://rformassspectrometry.github.io/PTMods/reference/addVariableModifications.md)
functions.

The latter will provide all combinations of possible sequences with the
given modifications.

\
[`addFixedModifications`](https://rformassspectrometry.github.io/PTMods/reference/addFixedModifications.md)`(``"MPEPTIDE"``,`\
`    fixedModifications ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``M ``=`` ``"Oxidation"``)`\
`)`

    ## [1] "M[Oxidation]PEPTIDE"

\
[`addVariableModifications`](https://rformassspectrometry.github.io/PTMods/reference/addVariableModifications.md)`(``"MPEPTIDE"``,`\
`    variableModifications ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``M ``=`` ``-``15.995``, T ``=`` ``"Phospho"``)`\
`)`

    ## [1] "MPEPTIDE"                   "M[-15.995]PEPTIDE"         
    ## [3] "MPEPT[Phospho]IDE"          "M[-15.995]PEPT[Phospho]IDE"

For fixed modifications, N- and C-terminal modifications are also
accepted:

\
[`addFixedModifications`](https://rformassspectrometry.github.io/PTMods/reference/addFixedModifications.md)`(``"MPEPT[Phospho]IDE"``,`\
`    fixedModifications ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``Nterm ``=`` ``304``)`\
`)`

    ## [1] "[+304]-MPEPT[Phospho]IDE"

And so are positional modifications, where one modification is applied
on the residue at position `pos`:

\
[`addFixedModifications`](https://rformassspectrometry.github.io/PTMods/reference/addFixedModifications.md)`(``"APEP[Phospho]TIDE"``,`\
`    fixedModifications ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"Fluoro"``)``,`\
`    pos ``=`` ``1`\
`)`

    ## [1] "A[Fluoro]PEP[Phospho]TIDE"

Notice that for both functions, the input sequence is kept as is
(keeping the already present PTM “Phospho”) whilst the notation style
provided is also untouched. Combined with
[`convertAnnotation()`](https://rformassspectrometry.github.io/PTMods/reference/convertAnnotation.md),
it is possible to convert the sequence into the preferred notation style
if need be.

For variable modifications, it is possible to limit the number of
possible combinations using the `maxMods` parameter:

\
[`addVariableModifications`](https://rformassspectrometry.github.io/PTMods/reference/addVariableModifications.md)`(``"MPEPTIDE"``,`\
`    variableModifications ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``M ``=`` ``-``15.995``, T ``=`` ``"Phospho"``)``,`\
`    maxMods ``=`` ``1`\
`)`

    ## [1] "MPEPTIDE"          "M[-15.995]PEPTIDE" "MPEPT[Phospho]IDE"

Both fixed and variable modifications can be called at once using
[`addModifications()`](https://rformassspectrometry.github.io/PTMods/reference/addModifications.md).
Use the `convertToStyle` parameter to convert the final sequence into
your preferred notation style if need be:

\
[`addModifications`](https://rformassspectrometry.github.io/PTMods/reference/addModifications.md)`(``"MPEPTIDE"``,`\
`    fixedModifications ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``M ``=`` ``"Oxidation"``)``,`\
`    variableModifications ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``T ``=`` ``79.966``)``,`\
`    convertToStyle ``=`` ``"unimodId"`\
`)`

    ## [1] "M[UNIMOD:35]PEPTIDE"            "M[UNIMOD:35]PEPT[UNIMOD:21]IDE"

To remove any and all modifications on a sequence, use
[`getCanonicalSequence()`](https://rformassspectrometry.github.io/PTMods/reference/getCanonicalSequence.md):

\
[`getCanonicalSequence`](https://rformassspectrometry.github.io/PTMods/reference/getCanonicalSequence.md)`(``"[+304]-M[Oxidation]PEPT[UNIMOD:21]IDE"``)`

    ## [1] "MPEPTIDE"

## License

This vignette is distributed under a [CC BY-SA
license](https://creativecommons.org/licenses/by-sa/4.0/).

## Session information

\
[`sessionInfo`](https://rdrr.io/r/utils/sessionInfo.html)`(``)`

    ## R version 4.6.1 (2026-06-24)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.5 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] PTMods_1.1.2     BiocStyle_2.41.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] cli_3.6.6           knitr_1.52          rlang_1.3.0        
    ##  [4] xfun_0.62           otel_0.2.0          textshaping_1.0.5  
    ##  [7] jsonlite_2.0.0      htmltools_0.5.9     ragg_1.5.2         
    ## [10] sass_0.4.10         rmarkdown_2.32      evaluate_1.0.5     
    ## [13] jquerylib_0.1.4     fastmap_1.2.0       yaml_2.3.12        
    ## [16] lifecycle_1.0.5     bookdown_0.48       BiocManager_1.30.27
    ## [19] compiler_4.6.1      fs_2.1.0            htmlwidgets_1.6.4  
    ## [22] systemfonts_1.3.2   digest_0.6.39       R6_2.6.1           
    ## [25] bslib_0.12.0        tools_4.6.1         pkgdown_2.2.1.9000 
    ## [28] cachem_1.1.0        desc_1.4.3
