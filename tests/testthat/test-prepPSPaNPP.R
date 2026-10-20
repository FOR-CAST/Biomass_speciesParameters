## prepPSPaNPP() keeps the trees of at least `minDBH` cm DBH. Until 3.1.0.9002 it filtered after the join
## that adds PSPplot's own `minDBH` column (each plot's recorded minimum DBH, as PSPclean returns it),
## and inside `[` that column masked the argument: each tree was compared with its plot's minimum, so
## the parameter had no effect (BC plots, which record trees from 2 cm, kept them at `minDBH = 9`).

## one plot measured once, with 10 trees each of 3, 8.9, 9 and 15 cm DBH (the forest filter needs
## >= 30 trees at the first measurement); `plotMinDBH` is the plot's own recorded minimum DBH
smallTreesPSP <- function(plotMinDBH) {
  dbh <- rep(c(3, 8.9, 9, 15), each = 10)
  list(
    measure = data.table::data.table(
      MeasureID = 1L, OrigPlotID1 = "p1", MeasureYear = 2000L, source = "BC", TreeNumber = seq_along(dbh),
      DBH = dbh, Height = 10, Species = "Picea glauca", status = "A", first_tree_year = 2000L,
      last_tree_year = 2010L, diff_dbh = 0),
    plot = data.table::data.table(
      MeasureID = 1L, OrigPlotID1 = "p1", MeasureYear = 2000L, source = "BC", minDBH = plotMinDBH,
      baseYear = 2000L, baseSA = 50, PlotSize = 0.1)
  )
}

## the DBH classes of the trees prepPSPaNPP() keeps
keptDBH <- function(minDBH, plotMinDBH) {
  p <- smallTreesPSP(plotMinDBH)
  out <- suppressMessages(prepPSPaNPP(
    studyAreaANPP = NULL, PSPgis = NULL, PSPmeasure = p$measure, PSPplot = p$plot,
    sppEquivLong = data.table::data.table(Latin_full = "Picea glauca", SpBiomassEq = "white spruce"),
    useHeight = TRUE, biomassModel = "Lambert2005", PSPperiod = c(1990, 2020), minDBH = minDBH))
  sort(unique(out$DBH))
}

test_that("regression: minDBH, not the plot's recorded minimum DBH, decides which trees are kept", {
  ## prepPSPaNPP calls biomassCalculation(), setkey() ... unqualified, as it does attached in a simList
  withr::local_package("pemisc")
  withr::local_package("data.table")
  ## a plot that records trees from 2 cm, as BC's do: the trees under 9 cm go at 9 and stay at 0
  expect_identical(keptDBH(minDBH = 0, plotMinDBH = 2), c(3, 8.9, 9, 15))
  expect_identical(keptDBH(minDBH = 9, plotMinDBH = 2), c(9, 15))
  ## a plot whose recorded minimum is 12 cm keeps the 9 cm trees it holds anyway
  expect_identical(keptDBH(minDBH = 9, plotMinDBH = 12), c(9, 15))
})
