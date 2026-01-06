test_that(".elements", {
    skip_if_not_installed("xml2")
    x <- xml2::read_xml(paste0(
        '<umod:unimod xmlns:umod="http://www.unimod.org/xmlns/schema/unimod_2">',
        '<umod:elements>',
        '<umod:elem title="H" full_name="Hydrogen" avge_mass="1.00794" mono_mass="1.007825035"/>',
        '<umod:elem title="2H" full_name="Deuterium" avge_mass="2.014101779" mono_mass="2.014101779"/>',
        '</umod:elements>',
        '</umod:unimod>'))
    d <- data.frame(
        Name=c("H", "2H"),
        FullName=c("Hydrogen", "Deuterium"),
        AvgMass=c(1.00794, 2.014101779),
        MonoMass=c(1.007825035, 2.014101779),
        row.names=c("H", "2H"),
        stringsAsFactors=FALSE
    )
    expect_equal(PTMods:::.elements(x), d)
})
