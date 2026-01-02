test_that(".aminoacids", {
    skip_if_not_installed("xml2")
    x <- xml2::read_xml(paste0(
        '<umod:unimod xmlns:umod="http://www.unimod.org/xmlns/schema/unimod_2">',
        '<umod:amino_acids>',
        '<umod:aa title="A" three_letter="Ala" full_name="Alanine" mono_mass="71.037114" avge_mass="71.0779">',
        '<umod:element symbol="H" number="5"/>',
        '<umod:element symbol="C" number="3"/>',
        '<umod:element symbol="N" number="1"/>',
        '<umod:element symbol="O" number="1"/>',
        '</umod:aa>',
        '<umod:aa title="C" three_letter="Cys" full_name="Cysteine" mono_mass="103.009185" avge_mass="103.1429">',
        '<umod:element symbol="H" number="5"/>',
        '<umod:element symbol="C" number="3"/>',
        '<umod:element symbol="N" number="1"/>',
        '<umod:element symbol="O" number="1"/>',
        '<umod:element symbol="S" number="1"/>',
        '</umod:aa>',
        '</umod:amino_acids>',
        '</umod:unimod>'
    ))
    d <- data.frame(
        OneLetter=c("A", "C"),
        ThreeLetter=c("Ala", "Cys"),
        FullName=c("Alanine", "Cysteine"),
        AvgMass=c(71.0779, 103.1429),
        MonoMass=c(71.037114, 103.009185),
        H=c(5L, 5L),
        C=c(3L, 3L),
        N=c(1L, 1L),
        O=c(1L, 1L),
        S=c(0L, 1L),
        row.names=c("A", "C"),
        stringsAsFactors=FALSE
    )
    expect_equal(PTMods:::.aminoacids(x), d)
})
