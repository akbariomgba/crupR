#' Clusters the active enhancer into condition-specific groups
#'
#' @param data List containing the metadata of the experiments and
#' the bin-wise enhancer prediction values as a GRanges object
#' (calculated in the prior step)
#' @param w_0 The minimum difference between the normalized prediction means
#' that two enhancers need in order to be included in the clustering.
#' Default is 0.5.
#' @param threshold Threshold for the p-values calculated during the clustering.
#'  Default is 0.5.
#' @param window Number of bins +/- the current bin that should be included
#' when calculating the p-values. Default is 10.
#' @param cores Number of cores to use
#' @return A list containing the meta data of the experiments and
#' the condition-specific clusters in a GRanges object
#' @examples
#'
#' files <- c(system.file("extdata", "Condition1.H3K4me1.bam", package="crupR"),
#'            system.file("extdata", "Condition1.H3K4me3.bam", package="crupR"),
#'            system.file("extdata", "Condition1.H3K27ac.bam", package="crupR"),
#'            system.file("extdata", "Condition2.H3K4me1.bam", package="crupR"),
#'            system.file("extdata", "Condition2.H3K4me3.bam", package="crupR"),
#'            system.file("extdata", "Condition2.H3K27ac.bam", package="crupR"))
#' in <- rep(system.file("extdata", "Condition1.Input.bam", package="crupR"),3)
#' in2 <- rep(system.file("extdata","Condition2.Input.bam", package="crupR"),3)
#' metaData <- data.frame(HM = rep(c("H3K4me1", "H3K4me3", "H3K27ac"),2),
#'                 condition = c(1,1,1,2,2,2), replicate = c(1,1,1,1,1,1),
#'                 bamFile = files, inputFile = c(in, in2))
#' metaData1 <- subset(metaData, condition == 1)
#' metaData2 <- subset(metaData, condition == 2)
#' pred1 <- readRDS(system.file("extdata", "condition1_predictions.rds",
#'                                package = "crupR"))
#' pred1 <- list("metaData" = metaData1, "data_matrix" = pred1)
#' pred2 <- readRDS(system.file("extdata", "condition2_predictions.rds",
#'                                package = "crupR"))
#' pred2 <- list("metaData" = metaData2, "data_matrix" = pred2)
#' predictions <- list(pred1, pred2)
#' getDynamics(data = predictions, cores = 2)
#'
#' @export
#' @importFrom GenomicRanges GRanges elementMetadata merge mcols granges start end reduce findOverlaps
#' @importFrom dplyr mutate lag filter setdiff
#' @importFrom matrixStats rowVars colMedians colMaxs
#' @import parallel
#' @importFrom IRanges values
#' @importFrom magrittr extract %>%
#' @importFrom stats ks.test
#' @importFrom utils combn getFromNamespace
#' @importFrom S4Vectors queryHits subjectHits

getDynamics <- function(data, w_0 = 0.5, threshold = 0.05, window = 10, cores = 1){

  ##################################################################
  # start run time
  ##################################################################

  start_time <- Sys.time()

  ##################################################################
  # define fixed parameters
  ##################################################################

  GR_header <- c("seqnames", "start","end","width")

  # axis and legend labels:
  legend_label <- "Probabilities"
  y_axis_label <- "Differential Enhancer"
  x_axis_label <- ""

  # prefixes:
  condition_prefix  <- "Condition "
  ID_prefix         <- "cond"
  sample_prefix     <- "Rep"

  ##################################################################
  # check input parameter
  ##################################################################

  # get meta data

  metaData <- do.call("rbind", lapply(data, function(x) x$metaData))
  conds <- unique(metaData$condition)
  num_cond <- length(conds)
  IDs <- list()
  for(i in seq_len(num_cond)){
    cond_subset = subset(metaData, condition == conds[i])
    IDs[[i]] <- paste0(ID_prefix,conds[i], "_", unique(cond_subset$replicate))
  }
  comb <- (combn(seq(length(unique(metaData$condition))),2))
  comb.list <- lapply(seq_len(ncol(comb)), function(i) comb[,i])


  if(threshold < 0 | threshold > 1){
    message <- paste0(threshold, " is not in range [0, 1].")
    stop(message);
  }


  if(w_0 < 0 | w_0 > 1){
    message <- paste0(w_0, " is not in range [0,1].")
    stop(message);
  }

  ##################################################################
  # list input parameter
  ##################################################################

  startPart("List input parameter")

  cat(skip(), "w_0: ",w_0, "\n")
  cat(skip(), "threshold: ",threshold, "\n")
  cat(skip(), "cores: ",cores, "\n")
  cat(skip(), "window: ",window, "\n")

  endPart()


  ##################################################################
  # read and combine probabilities:
  ##################################################################

  startPart("Read enhancer probabilites for all samples and conditions")
  probs <- GenomicRanges::GRanges()
  ids <- unlist(IDs)
  #print(ids)

  for(i in seq_len(length(data))){
    this <- data[[i]]
    matrix <- this$data_matrix
    colnames(GenomicRanges::elementMetadata(matrix)) <- ids[i]
    probs <- GenomicRanges::merge(probs, matrix, all = TRUE)
  }


  endPart()

  ##################################################################
  # get empricial p-values
  # (list with empirical p-values per pair comparison and
  # direction of group mean differences)
  ##################################################################

  startPart("Calculate (pairwise) empirical p-values")
  p <- get_pairwisePvalues(p = probs, I = IDs, w_0 = w_0, W = window,
                           p.thres = threshold, C = cores)
  endPart()

  ##################################################################
  # Call dynamically changing enhancer peaks:
  ##################################################################

  startPart("Get condition-specific enhancer peaks")
  cat(paste0(skip(), "build significance peak pattern"))
  probs <- get_cluster(p = probs, pvalues = p, I = IDs)
  done()

  if ((length(unique(probs$pattern)) == 1) && (unique(probs$pattern) == paste(rep("0", num_cond),collapse=""))) {
    message <- paste0(skip(),
                    "no significant peaks found between any two conditions.\n")
    stop(message);
  }

  cat(paste0(skip(), "combine peaks by significance pattern"))
  peaks <- get_ranges(p = probs, I = IDs, W = window, C = cores)
  rm(probs)
  done()

  cat(paste0(skip(), "assign clusters"))
  peaks <- name_patterns(gr = peaks, C = cores)
  done()

  ##################################################################
  # print run time
  ##################################################################

  out.list = list()
  out.list[['metaData']] <- metaData
  out.list[['sumFile']] <- peaks
  run_time <- Sys.time() - start_time
  startPart("Run time")
  cat(paste0(skip(), format(run_time), "\n"))
  endPart()

  return(out.list)

}
