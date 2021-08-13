#devtools::load_all("/project/wig/persia/CRUP_package/CRUP/R/")
context("Find condition-specific enhancers")
library(GenomicRanges)
library(rtracklayer)

##################################################################
# create input data
##################################################################

files <- c(system.file("extdata", "Condition1.H3K4me1.bam", package="crupR"),
           system.file("extdata", "Condition1.H3K4me3.bam", package="crupR"),
           system.file("extdata", "Condition1.H3K27ac.bam", package="crupR"),
           system.file("extdata", "Condition2.H3K4me1.bam", package="crupR"),
           system.file("extdata", "Condition2.H3K4me3.bam", package="crupR"),
           system.file("extdata", "Condition2.H3K27ac.bam", package="crupR"))

inputs <- rep(system.file("extdata", "Condition1.Input.bam", package="crupR"), 3)

inputs2 <- rep(system.file("extdata", "Condition2.Input.bam", package="crupR"), 3)

metaData <- data.frame(HM = rep(c("H3K4me1","H3K4me3","H3K27ac"),2),
                       condition = c(1,1,1,2,2,2), replicate = c(1,1,1,1,1,1),
                       bamFile = files, inputFile = c(inputs, inputs2))

metaData1 <- subset(metaData, condition == 1)
metaData2 <- subset(metaData, condition == 2)

#pred1 <-get(load("condition1.predictions.Rdata"))
pred1 <- readRDS(system.file("extdata", "condition1_predictions.rds", package = "crupR"))
#pred1 <- c(list("metaData" = metaData1), pred1)

#c1.half1 <- rtracklayer::import(system.file("extdata", "condition1_prediction-part1.bw", package = "crupR"),
#                                format = "BigWig")
#c1.half2 <- rtracklayer::import(system.file("extdata", "condition1_prediction-part2.bw", package = "crupR"),
#                                 format = "BigWig")
#c1.pred <- suppressWarnings(c(c1.half1, c1.half2))
pred1 <- list("metaData" = metaData1, "data_matrix" = pred1)
#pred2 <-get(load("condition2.predictions.Rdata"))
pred2 <- readRDS(system.file("extdata", "condition2_predictions.rds", package = "crupR"))
#pred2 <- c(list("metaData" = metaData2), pred2)
#c2.half1 <- rtracklayer::import(system.file("extdata", "condition2_prediction-part1.bw", package = "crupR"),
#                                format = "BigWig")
#c2.half2 <- rtracklayer::import(system.file("extdata", "condition2_prediction-part2.bw", package = "crupR"),
#                                format = "BigWig")
#c2.pred <- suppressWarnings(c(c2.half1, c2.half2))
pred2 <- list("metaData" = metaData2, "data_matrix" = pred2)

preds <- list(pred1, pred2)

##################################################################
# test enhancerDynamics()
##################################################################

testthat::test_that("the error messages of getDynamics() work",{

  testthat::expect_error(crupR::getDynamics(data = preds, w_0 = 1.2, cores = 10),
                         "1.2 is not in range [0,1].", fixed = TRUE)

  testthat::expect_error(crupR::getDynamics(data = preds, threshold = -0.5, cores = 10),
                         "-0.5 is not in range [0, 1].", fixed = TRUE)
})

testthat::test_that("getDynamics runs as expected",{
  #dynamics.expected <- get(load("differential_enhancers.Rdata"))
  dynamics.expected <- readRDS(system.file("extdata", "differential_enhancers.rds", package = "crupR"))
  dynamics <- crupR::getDynamics(data = preds, cores = 2)
  testthat::expect_equal(length(dynamics), 2)
  testthat::expect_equal(length(dynamics$sumFile), 1)
  testthat::expect_equal(dynamics$sumFile$cond1_1, dynamics.expected$sumFile$cond1_1, tolerance = 1e-9)
})



