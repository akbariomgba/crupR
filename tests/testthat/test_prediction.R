#devtools::load_all("/project/wig/persia/CRUP_package/CRUP/R/")
context("Enhancer Prediction")
library(GenomicRanges)
library(rtracklayer)
##################################################################
# create input data
##################################################################

files.short <- c(system.file("extdata", "Condition2.H3K4me1.bam", package="crupR"),
           system.file("extdata", "Condition2.H3K4me3.bam", package="crupR"),
           system.file("extdata", "Condition2.H3K27ac.bam", package="crupR"))

inputs <- rep(system.file("extdata", "Condition2.Input.bam", package="crupR"), 3)


#create the metaData frame
metaData.short <- data.frame(HM = c("H3K4me1","H3K4me3","H3K27ac"),
                       condition = c(2,2,2), replicate = c(1,1,1),
                       bamFile = files.short, inputFile = inputs)


#the expected values for the prediction
#pred.half1 <- rtracklayer::import(system.file("extdata", "condition2_prediction-part1.bw", package = "crupR"),
 #                                 format = "BigWig")
#pred.half2 <- rtracklayer::import(system.file("extdata", "condition2_prediction-part2.bw", package = "crupR"),
#                                  format = "BigWig")
#pred.expected <- suppressWarnings(c(pred.half1, pred.half2))

pred.expected <- readRDS(file = system.file("extdata", "condition2_predictions.rds", package="crupR"))
#the actual normalized counts

#data_matrix <- get(load("condition1_normalized.Rdata"))
#data_matrix <- readRDS(system.file("extdata", "condition2_normalized.rds", package = "crupR"))
#mcols_counts <- readRDS(system.file("extdata", "condition2_normalized_mcols.rds", package = "crupR"))
data_matrix <- readRDS(file = system.file("extdata", "condition2_normalized.rds", package="crupR"))
#GenomicRanges::mcols(data_matrix) <- mcols_counts
#create a list
norm = list(metaData = metaData.short, normalized = data_matrix)

##################################################################
# test enhancerPrediction()
##################################################################
#pred.expected <- readRDS(system.file("extdata", "condition2_predictions.rds", package = "crupR"))

testthat::test_that("the error messages of getEnhancers() work", {

  testthat::expect_error(crupR::getEnhancers(data = norm, classifier = "/project/wig/persia/classifier/",cores = 4),
                         "The directory of the classifier /project/wig/persia/classifier/ is not a valid directory")

  testthat::expect_error(crupR::getEnhancers(data = norm, cutoff = 1.3, cores = 4),
                         "1.3 is not in range [0,1].", fixed = TRUE)
})

testthat::test_that("getEnhancers() runs as expected",{
  pred <- crupR::getEnhancers(data = norm, cores = 2)
  pred_short <- pred$data_matrix[which(as.character(GenomicRanges::seqnames(pred$data_matrix)) == "chr8")]
  testthat::expect_equal(length(pred), 2)
  testthat::expect_equal(pred_short$score, pred.expected$score, tolerance = 1e-5)
})

