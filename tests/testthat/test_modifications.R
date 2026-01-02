context("modifications")

skip_if_not_installed("xml2")
x <- xml2::read_xml(paste0(
    '<umod:unimod xmlns:umod="http://www.unimod.org/xmlns/schema/unimod_2">',
    '<umod:modifications>',
     '<umod:mod title="Acetyl" full_name="Acetylation" username_of_poster="unimod" group_of_poster="admin" date_time_posted="2002-08-19 19:17:11" date_time_modified="2008-02-15 05:20:02" approved="1" record_id="1">',
    '<umod:specificity hidden="1" site="S" position="Anywhere" classification="Post-translational" spec_group="4"/>',
     '<umod:specificity hidden="1" site="C" position="Anywhere" classification="Post-translational" spec_group="3"/>',
    '<umod:specificity hidden="0" site="N-term" position="Any N-term" classification="Multiple" spec_group="2">',
    '<umod:misc_notes>GIST acetyl light</umod:misc_notes>',
    '</umod:specificity>',
    '<umod:delta mono_mass="42.010565" avge_mass="42.0367" composition="H(2) C(2) O"> <umod:element symbol="H" number="2"/> <umod:element symbol="C" number="2"/> <umod:element symbol="O" number="1"/>',
    '</umod:delta>',
    '<umod:xref>',
    '<umod:text>14730666</umod:text>',
    '<umod:source>PubMed PMID</umod:source>',
    '<umod:url/>',
    '</umod:xref>',
    '<umod:xref>',
    '<umod:text>15350136</umod:text>',
    '<umod:source>PubMed PMID</umod:source>',
    '<umod:url/>',
    '</umod:xref>',
    '</umod:mod>',
    '</umod:modifications>',
    '</umod:unimod>'))
node <- xml2::xml_find_first(x, "//umod:mod[@record_id=\"1\"]")

sp <- cbind(
    hidden=c("0", "1", "1"),
    site=c("N-term", "C", "S"),
    position=c("Any N-term", rep("Anywhere", 2)),
    classification=c("Multiple", rep("Post-translational", 2)),
    spec_group=as.character(2:4)
)

test_that(".characterId", {
    expect_error(PTMods:::.characterId(1:10, "foo", "bar", "0"), "name")
    expect_error(PTMods:::.characterId("foo", 1:10, "bar", "0"), "site")
    expect_error(PTMods:::.characterId("foo", "bar", 1:10, "0"), "position")
    expect_error(PTMods:::.characterId("foo", "bar", "0", 1:10), "neutralLoss")
    expect_equal(
        PTMods:::.characterId(
            rep("Acetyl", 4),
            c("K", "S", "N-term", "N-term"),
            c("Anywhere", "Anywhere", "Any N-term", "Protein N-term"),
            c("0", "1", "0", "0")
        ),
        c("Acetyl:K", "Acetyl:S:NL", "Acetyl:N-term", "Acetyl:P-N-term")
    )
    expect_equal(
        PTMods:::.characterId(
            rep("Acetyl", 4),
            c("K", "S", "N-term", "N-term"),
            c("Anywhere", "Anywhere", "Any N-term", "Any N-term"),
            c("0", "1", "0", "0")
        ),
        c("Acetyl:K", "Acetyl:S:NL", "Acetyl:N-term:1", "Acetyl:N-term:2")
    )
})

test_that(".delta", {
    delta <- c(avge_mass="42.0367", mono_mass="42.010565",
               composition="H(2) C(2) O")
    expect_equal(PTMods:::.delta(node), delta)
})

test_that(".modifications", {
    m <- data.frame(Id=c("Acetyl:N-term", "Acetyl:C", "Acetyl:S"),
                    UnimodId=1, Name="Acetyl", Description="Acetylation",
                    Composition="H(2) C(2) O",
                    AvgMass=42.0367, MonoMass=42.010565,
                    Site=c("N-term", "C", "S"),
                    Position=factor(c("Any N-term", "Anywhere", "Anywhere")),
                    Classification=factor(c("Multiple", "Post-translational",
                                            "Post-translational")),
                    SpecGroup=2:4,
                    NeutralLoss=FALSE,
                    LastModified="2008-02-15 05:20:02",
                    Approved=TRUE, Hidden=c(FALSE, TRUE, TRUE),
                    stringsAsFactors=FALSE)
    rownames(m) <- m$Id
    expect_equal(PTMods:::.modifications(x), m)
})

test_that(".neutralLoss", {
    node <- xml2::read_xml(paste0(
'<umod:unimod xmlns:umod="http://www.unimod.org/xmlns/schema/unimod_2">',
'<umod:modifications>',
'<umod:mod title="Phospho" full_name="Phosphorylation" username_of_poster="unimod"',
'   group_of_poster="admin"',
'   date_time_posted="2002-08-19 19:17:11"',
'   date_time_modified="2011-11-25 10:55:54"',
'   approved="1"',
'   record_id="21">',
'<umod:specificity hidden="0" site="T" position="Anywhere" classification="Post-translational" spec_group="1">',
'<umod:NeutralLoss mono_mass="97.976896" avge_mass="97.9952" flag="false" composition="H(3) O(4) P">',
'<umod:element symbol="H" number="3"/>',
'<umod:element symbol="O" number="4"/>',
'<umod:element symbol="P" number="1"/>',
'</umod:NeutralLoss>',
'<umod:NeutralLoss mono_mass="0.000000" avge_mass="0.0000" flag="false" composition="0"/>',
'</umod:specificity>',
'<umod:delta mono_mass="79.966331" avge_mass="79.9799" composition="H O(3) P">',
'<umod:element symbol="H" number="1"/>',
'<umod:element symbol="O" number="3"/>',
'<umod:element symbol="P" number="1"/>',
'</umod:delta>',
'</umod:mod>',
'</umod:modifications>',
'</umod:unimod>'))
    nl <- matrix(c("T", "1", "97.9952", "97.976896", "H(3) O(4) P"), nrow=1,
                 dimnames=list(c(), c("site", "spec_group",
                                      "avge_mass", "mono_mass", "composition")))
    expect_equal(PTMods:::.neutralLoss(node), nl)
})

test_that(".specificity", {
    expect_equal(PTMods:::.specificity(node), sp)

    sp2 <- cbind(
        hidden="1",
        site="S",
        position="Anywhere",
        classification="Post-translational",
        spec_group="4"
    )
    node2 <- xml2::xml_find_first(xml2::read_xml(paste0(
'<umod:unimod xmlns:umod="http://www.unimod.org/xmlns/schema/unimod_2">',
'<umod:modifications>',
'<umod:mod title="Acetyl" full_name="Acetylation" username_of_poster="unimod" group_of_poster="admin" date_time_posted="2002-08-19 19:17:11" date_time_modified="2008-02-15 05:20:02" approved="1" record_id="1">',
'<umod:specificity hidden="1" site="S" position="Anywhere" classification="Post-translational" spec_group="4"/>',
'</umod:mod>',
'</umod:modifications>',
'</umod:unimod>')), "//umod:mod")
    expect_equal(PTMods:::.specificity(node2), sp2)
})

test_that(".title", {
    title <- c(record_id="1", title="Acetyl", full_name="Acetylation",
               date_time_modified="2008-02-15 05:20:02", approved="1")
    expect_equal(PTMods:::.title(node), title)
})
