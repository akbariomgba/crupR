#' Predicts the occurence of enhancers using the normalized HMs counts.
#'
#' @param data List containing the metaData of the experiments that were normalized and also the normalized counts as GRanges object (basically the output of the normalize() step)
#' @param classifier The path of the classifier to use for the prediction. When set to default, the default clasifier is used.
#' @param cutoff Cutoff for the probabilities. Default is 0.5.
#' @param distance Maximum distance (bp) for peak clustering. Default is 12500
#' @param cores Number of cores to use
#' @return A list containing the meta data of the experiments whose enhancers were predicted and
#' the enhancer probabilities for each 100 bp bin in the genome as a GRanges object
#' @examples
#' #first recreate the output of crupR::normalize (so skip this)
#' files <- c(system.file("extdata", "Condition2.H3K4me1.bam", package="crupR"),
#'           system.file("extdata", "Condition2.H3K4me3.bam", package="crupR"),
#'           system.file("extdata", "Condition2.H3K27ac.bam", package="crupR"))
#' inputs <- rep(system.file("extdata", "Condition2.Input.bam", package="crupR")) 
#' metaData <- data.frame(HM = c("H3K4me1","H3K4me3","H3K27ac"),
#'                     condition = c(2,2,2), replicate = c(1,1,1),
#'                     bamFile = files, inputFile = inputs)
#' data_matrix <- readRDS(system.file("extdata", "condition2_normalized.rds", package="crupR"))
#' norm <- list(metaData = metaData, normalized = data_matrix)
#' #let's run the actual function
#' getEnhancers(data = norm, cores = 2)
#'
#' @export
#' @importFrom preprocessCore normalize.quantiles.use.target
#' @import randomForest
#' @importFrom stats predict
#' @importFrom GenomicRanges mcols elementMetadata end reduce
#' @import parallel

getEnhancers <- function(data, classifier = "default", cutoff = 0.5, distance = 12500, cores = 1){
  ##################################################################
  # start run time
  ##################################################################

  start_time <- Sys.time()

  ##################################################################
  # check input parameter
  ##################################################################

  # check classifier directory
  if ((classifier != "default") & (!dir.exists(classifier))){
    message <- paste0("The directory of the classifier ",classifier, " is not a valid directory")
    stop(message);
  }

  # check cutoff
   if (!is.null(cutoff) & (cutoff < 0 | cutoff > 1) ) {
    message <- paste0(cutoff, " is not in range [0,1].")
    stop(message);
  }

  ##################################################################
  # define input parameter
  ##################################################################

  startPart("List input parameter")

  cat(skip(), "classifier: ",classifier, "\n")
  cat(skip(), "distance: ",distance, "\n")
  cat(skip(), "cutoff: ",cutoff, "\n")
  cat(skip(), "cores: ",cores, "\n")

  endPart()

  ##################################################################
  # Read classifier files and ecdf
  ##################################################################

  startPart("Get classifier and empirical distribution function")

  if(classifier != "default"){
    classifier_file1 <- file.path(classifier, "active_vs_inactive.rds")
  }else{
    classifier_file1 <- system.file("extdata", "active_vs_inactive.rds", package = "crupR")
  }
  check_file(classifier_file1)
  classifier1 <- readRDS(classifier_file1)
  features1 <- unique(gsub('_.*','',names(classifier1$forest$xlevels)))


  if(classifier != "default"){
    classifier_file2 <- file.path(classifier, "enhancer_vs_active_promoter.rds")
  }else{
    classifier_file2 <- system.file("extdata", "enhancer_vs_active_promoter.rds", package = "crupR")
  }
  check_file(classifier_file2)
  classifier2 <- readRDS(classifier_file2)
  features2 <- unique(gsub('_.*','',names(classifier2$forest$xlevels)))

  features_all <- unique(c(features1, features2))

  # ecdf: will be used for quantile normalization
  ecdf_file <- system.file("extdata", "ecdf.rds", package = "crupR")
  check_file(ecdf_file)
  ecdf <- readRDS(ecdf_file)

  endPart()

  ##################################################################
  # Quantile normalization
  ##################################################################

  startPart("Quantile normalize summarized data matrix")

  # load data file
  data_matrix <-data$normalized

  # normalize for each feature in data matrix
  data_matrix_norm <- data_matrix
  for (feature in features_all) {
    cat(paste0(skip(), "quantile normalize counts for feature ", feature))

    feature_norm <- get_targetQuantileNorm(ecdf[[feature]])
    GenomicRanges::mcols(data_matrix_norm)[,feature] <- preprocessCore::normalize.quantiles.use.target(matrix(GenomicRanges::mcols(data_matrix_norm)[,feature]),
                                                                        feature_norm)
    done()
  }

  endPart()


  ##################################################################
  # Extend data matrix
  ##################################################################

  startPart("Extend data matrix (+/- 5 bins)")

  # create extended data matrix
  cat(paste0(skip(), "create extended data matrix"))

  # original matrix
  data_matrix_ext  <- extend_dataMatrix(N = 5, df = data.frame(data_matrix), f = features1)#?
  zero.idx <- which(rowSums(data_matrix_ext[,-c(seq_len(3))]) == 0)
  rm(data_matrix_ext)

  # normalized matrix
  data_matrix_norm_ext <- extend_dataMatrix(N = 5, df = data.frame(data_matrix_norm), f = features_all)
  data_matrix_norm_ext <- data_matrix_norm_ext[-zero.idx,]
  done()

  endPart()

  ##################################################################
  # make predictions
  ##################################################################

  startPart("Get enhancer probabilities for each bin")

  mid=round(nrow(data_matrix_norm_ext)/2,0)
  prediction1 <- predict(classifier1, data_matrix_norm_ext[seq_len(mid),], type = "prob")[,2]
  prediction1 <- c(prediction1, predict(classifier1, data_matrix_norm_ext[(mid+1):nrow(data_matrix_norm_ext),], type = "prob")[,2])

  prediction2 <- predict(classifier2, data_matrix_norm_ext[seq_len(mid),], type = "prob")[,2]
  prediction2 <- c(prediction2, predict(classifier2, data_matrix_norm_ext[(mid+1):nrow(data_matrix_norm_ext),], type = "prob")[,2])

  rm(data_matrix_norm_ext)

  GenomicRanges::elementMetadata(data_matrix) <- NULL
  GenomicRanges::mcols(data_matrix)["prob"] <- 0
  GenomicRanges::mcols(data_matrix[-zero.idx])["prob"] <- prediction1 * prediction2

  endPart()

  ##################################################################
  # save predicitons
  ##################################################################
  #create here a new granges object  ==> how it should be in bw format
   score.matrix <- data_matrix
   colnames(GenomicRanges::elementMetadata(data_matrix)) <- "score"
   #print("old seqlenghts")
   #print(GenomeInfoDb::seqlengths(data_matrix))
   #print("supplied")
   #print(GenomicRanges::end(GenomicRanges::reduce(data_matrix)))
   #GenomeInfoDb::seqlengths(data_matrix) <- GenomicRanges::end(GenomicRanges::reduce(data_matrix))
   #print("new seqlengths")
   #print(GenomeInfoDb::seqlengths(data_matrix))
  ##################################################################
  # print run time
  ##################################################################

  run_time <- Sys.time() - start_time

  startPart("Run time")
  cat(paste0(skip(), format(run_time), "\n"))
  endPart()
  outlist <- list("metaData" = data$metaData ,"data_matrix" = data_matrix)
  return(outlist)
}
