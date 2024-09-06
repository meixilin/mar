test_that("MARcalc reproduces arabidopsis results", {
    zM = 0.324
    zMci = c(0.238, 0.41)

    # load old mares results
    load('testdata/tmpobjects/mares-arabidopsis.rda')
    # load new marsampling output
    load('testdata/mardf-arabidopsis.rda')
    colnames(mares)[6:7] = c('Asq','A')
    oldmar = MARcalc(mares, Atype = 'Asq')
    oldemar = MARcalc(mares, Mtype = 'E', Atype = 'Asq')
    newmar = MARcalc(mardf, Atype = "Asq")
    newmarz = newmar$sigConf[2,1]
    # print(newmar)
    # print(oldmar)
    expect_equal(oldmar$sigConf[2,1], zM, tolerance = 0.0005)
    expect_equal(unname(oldmar$sigConf[2,5:6]), zMci, tolerance = 0.0005)
    expect_true(newmarz > zMci[1] & newmarz < zMci[2])
})

test_that("MARcalc reproduces joshua results", {
    zM = 0.128
    zMci = c(0.109, 0.147)

    # load old mares results
    load('testdata/tmpobjects/mares-joshua.rda')
    # load new marsampling output
    load('testdata/mardf-joshua.rda')
    colnames(mares)[6:7] = c('Asq','A')
    oldmar = MARcalc(mares, Atype = 'Asq')
    newmar = MARcalc(mardf, Atype = "Asq")
    newmarz = newmar$sigConf[2,1]
    expect_equal(oldmar$sigConf[2,1], zM, tolerance = 0.0005)
    expect_equal(unname(oldmar$sigConf[2,5:6]), zMci, tolerance = 0.0005)
    expect_true(newmarz > zMci[1] & newmarz < zMci[2])
})




