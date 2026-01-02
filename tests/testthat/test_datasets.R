context("datasets")

# may move to longtests/
# see: https://stat.ethz.ch/pipermail/bioc-devel/2017-November/012327.html

skip_if_not_installed("xml2")

xml <- PTMods:::.unimodDb()

test_that("aminoacids", {
    data(aminoacids, package="PTMods")
    expect_equal(PTMods:::.aminoacids(xml), aminoacids)
})

test_that("elements", {
    data(elements, package="PTMods")
    expect_equal(PTMods:::.elements(xml), elements)
})

test_that("modifications", {
    data(modifications, package="PTMods")
    expect_equal(PTMods:::.modifications(xml), modifications)
})
