## buildGrowthCurves() picks the PSP plot-years each growth curve is fitted to, and buildModels() fits
## a curve only if there are at least `minimumPlots` of them. Until 3.0.2.9005 two things were wrong:
## - the check counted rows of species x plot-year. Under "focal", a plot-year where another species
##   was also co-dominant (> 20% of plot biomass) counted twice: once for the focal species and once
##   for the other species, pooled as "Other", whose rows are not fitted. Under "pairwise", every
##   plot-year counted once per species of the pair.
## - species codes were matched as substrings, so a code that is part of another (Pice_eng in
##   Pice_eng_gla) drew that species' plots.

latinOf <- c(Abie_las = "Abies lasiocarpa", Pice_eng = "Picea engelmannii",
             Pice_eng_gla = "Picea engelmannii x glauca", Pice_mar = "Picea mariana",
             Popu_tre = "Populus tremuloides")

## sppEquiv for the LandR codes `spp`
equivFor <- function(spp) data.table::data.table(Latin_full = unname(latinOf[spp]), LandR = spp)

## PSP records as prepPSPaNPP() hands them to buildGrowthCurves(), with the columns it reads or drops:
## here one record per species x plot x measurement year, with that species' biomass on a 1-ha plot
pspRecords <- function(species, plot, year, age, biomass) {
  data.table::data.table(
    MeasureID = paste(plot, year, sep = "_"), OrigPlotID1 = plot, MeasureYear = year,
    Latin_full = unname(latinOf[species]), SpBiomassEq = "x", source = "x", PSP = "x",
    Elevation = 0, PlotSize = 1, standAge = age, biomass = biomass)
}

## biomass along a growth curve the fits pin down at once
growth <- function(age) 20000 * (1 - exp(-0.02 * age))^2 * exp(stats::rnorm(length(age), 0, 0.1))

## Pice_mar on 17 plots, 8 of them measured twice: 25 plot-years. Popu_tre is co-dominant (a third of
## plot biomass) in 17 of them, those on the first 9 plots. The old count gave Pice_mar 25 + 17 = 42
## rows (as for black spruce in LandWeb's WesternAlbertaUpland: 42 rows, 25 plot-years on 17 plots).
focalPSP <- function() {
  withr::local_seed(1)
  plots <- sprintf("p%02d", 1:17)
  plot <- c(plots, plots[1:8])
  year <- rep(c(2000L, 2010L), c(17L, 8L))
  age <- seq(15, 175, by = 10)[match(plot, plots)] + (year - 2000L)
  biomass <- growth(age)
  coDom <- plot %in% plots[1:9]
  rbind(pspRecords("Pice_mar", plot, year, age, biomass),
        pspRecords("Popu_tre", plot[coDom], year[coDom], age[coDom], biomass[coDom] / 2))
}

fitFocal <- function(psp, minimumPlots) {
  suppressMessages(buildGrowthCurves(psp, speciesCol = "LandR",
                                     sppEquiv = equivFor(c("Pice_mar", "Popu_tre")),
                                     quantileAgeSubset = 99, minimumSampleSize = minimumPlots,
                                     speciesFittingApproach = "focal"))
}

test_that("minimumPlots counts the focal species' plot-years, not the rows of co-dominant 'Other'", {
  ## regression: the old count gave Pice_mar 42 >= 40 rows, so its curve was fitted
  withr::local_package("data.table")
  w <- capture_warnings(gcs <- fitFocal(focalPSP(), minimumPlots = 40))
  expect_identical(unname(gcs$Pice_mar), "insufficient data")
  expect_identical(unname(gcs$Popu_tre), "insufficient data")
  expect_match(w, "Pice_mar: 25 plot-years", fixed = TRUE, all = FALSE)
  expect_match(w, "Popu_tre: 17 plot-years", fixed = TRUE, all = FALSE)
})

test_that("a species with exactly minimumPlots plot-years is fitted, with one fewer it is not", {
  withr::local_package("data.table")
  psp <- focalPSP()
  gcs <- suppressWarnings(fitFocal(psp, minimumPlots = 25))
  expect_named(gcs$Pice_mar, c("originalData", "simData", "NonLinearModel"))
  expect_s3_class(gcs$Pice_mar$NonLinearModel$Pice_mar, "nlrob")
  ## regression: 17 plot-years, but 34 rows with the co-dominant Pice_mar as "Other"
  expect_identical(unname(gcs$Popu_tre), "insufficient data")
  gcs <- suppressWarnings(fitFocal(psp, minimumPlots = 26))
  expect_identical(unname(gcs$Pice_mar), "insufficient data")
})

test_that("a species of sppEquiv without plot-years is not fitted", {
  withr::local_package("data.table")
  w <- capture_warnings(gcs <- suppressMessages(buildGrowthCurves(
    focalPSP(), speciesCol = "LandR", sppEquiv = equivFor(c("Pice_mar", "Popu_tre", "Abie_las")),
    quantileAgeSubset = 99, minimumSampleSize = 40, speciesFittingApproach = "focal")))
  expect_identical(unname(gcs$Abie_las), "insufficient data")
  expect_match(w, "Abie_las: 0 plot-years", fixed = TRUE, all = FALSE)
})

test_that("under pairwise, minimumPlots counts the pair's plot-years once, not once per species", {
  ## regression: the old count gave the pair 2 x 17 = 34 rows, so its curves were fitted
  withr::local_package("data.table")
  withr::local_seed(1)
  age <- seq(15, 175, by = 10)
  biomass <- growth(age)
  plots <- sprintf("p%02d", seq_along(age))
  psp <- rbind(pspRecords("Abie_las", plots, 2000L, age, biomass / 2),
               pspRecords("Pice_eng", plots, 2000L, age, biomass / 2))
  w <- capture_warnings(gcs <- suppressMessages(buildGrowthCurves(
    psp, speciesCol = "LandR", sppEquiv = equivFor(c("Abie_las", "Pice_eng")),
    quantileAgeSubset = 99, minimumSampleSize = 18, speciesFittingApproach = "pairwise")))
  expect_identical(unname(gcs[["Abie_las__Pice_eng"]]), "insufficient data")
  expect_match(w, "Abie_las__Pice_eng: 17 plot-years", fixed = TRUE, all = FALSE)
})

test_that("focal: a species code that is part of another does not draw that species' plots", {
  ## regression: matched as a substring, Pice_eng also drew the Pice_eng_gla plots, as "Other"
  withr::local_package("data.table")
  ## the PSP records each species' fit is given
  local_mocked_bindings(buildModels = function(species, psp, ...) psp)
  age <- seq(20, 100, by = 20)
  psp <- rbind(pspRecords("Pice_eng", paste0("eng", 1:5), 2000L, age, 100 * age),
               pspRecords("Pice_eng_gla", paste0("gla", 1:5), 2000L, age, 100 * age),
               pspRecords(c("Pice_eng", "Pice_eng_gla"), "mixed", 2000L, 60, c(4000, 3000)))
  out <- suppressMessages(buildGrowthCurves(
    psp, speciesCol = "LandR", sppEquiv = equivFor(c("Pice_eng", "Pice_eng_gla")),
    quantileAgeSubset = 99, minimumSampleSize = 1, speciesFittingApproach = "focal"))
  expect_setequal(unique(out$Pice_eng$OrigPlotID1), c(paste0("eng", 1:5), "mixed"))
  expect_setequal(unique(out$Pice_eng_gla$OrigPlotID1), c(paste0("gla", 1:5), "mixed"))
  ## the other species is "Other" only on the plot it shares with the focal species
  expect_identical(unique(out$Pice_eng$OrigPlotID1[out$Pice_eng$speciesTemp == "Other"]), "mixed")
})

test_that("pairwise: a pair whose name is part of another pair's does not draw that pair's plots", {
  ## regression: matched as a substring, "Abie_las__Pice_eng" also drew the plots of
  ## "Abie_las__Pice_eng_gla"
  withr::local_package("data.table")
  local_mocked_bindings(buildModels = function(species, psp, ...) psp)
  ## Abie_las is recorded first on every plot, so the pairs are named "Abie_las__<other>"
  stands <- function(other, plots) {
    rbind(pspRecords("Abie_las", plots, 2000L, 50, 5000), pspRecords(other, plots, 2000L, 50, 5000))
  }
  psp <- rbind(stands("Pice_eng", paste0("eng", 1:3)), stands("Pice_eng_gla", paste0("gla", 1:3)))
  out <- suppressMessages(buildGrowthCurves(
    psp, speciesCol = "LandR", sppEquiv = equivFor(c("Abie_las", "Pice_eng", "Pice_eng_gla")),
    quantileAgeSubset = 99, minimumSampleSize = 1, speciesFittingApproach = "pairwise"))
  expect_setequal(names(out), c("Abie_las__Pice_eng", "Abie_las__Pice_eng_gla"))
  expect_setequal(unique(out[["Abie_las__Pice_eng"]]$OrigPlotID1), paste0("eng", 1:3))
  expect_setequal(as.character(unique(out[["Abie_las__Pice_eng"]]$speciesTemp)),
                  c("Abie_las", "Pice_eng"))
  expect_setequal(unique(out[["Abie_las__Pice_eng_gla"]]$OrigPlotID1), paste0("gla", 1:3))
})
