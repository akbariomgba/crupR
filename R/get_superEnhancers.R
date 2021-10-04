#' Uses the enhancer predictions to find superenhancers.
#'
#' @param data List containing the meta data of the experiments and the
#' enhancer prediction values in a GRanges object calculated
#' by enhancerPrediction()
#' @param cutoff Cutoff for the enhancer probabilities. Default is 0.5.
#' @param distance Maximum distance (bp) for peak clustering. Default is 12500
#' @param cores Number of cores to use
#' @return A list containing the meta data of the experiments,
#' the enhancer prediction values as a GRanges object (like the input data),
#' the enhancer peak calls as a GRanges object (can be exported as a bedGraph)
#' and the cluster of peaks in a GRanges object (can be exported as bed file)
#' @examples
#'
#' # first recreate the output of crupR:getPredictions
#' files <- c(system.file('extdata', 'Condition2.H3K4me1.bam', package='crupR'),
#'         system.file('extdata', 'Condition2.H3K4me3.bam', package='crupR'),
#'         system.file('extdata', 'Condition2.H3K27ac.bam', package='crupR'))
#' inp <- rep(system.file('extdata', 'Condition2.Input.bam',package='crupR'), 3)
#' #create the metaData frame
#' metaData <- data.frame(HM = c('H3K4me1','H3K4me3','H3K27ac'),
#'                         condition = c(2,2,2), replicate = c(1,1,1),
#'                         bamFile = files, inputFile = inp)
#' pred <- readRDS(system.file('extdata', 'condition2_predictions.rds',
#'                         package='crupR'))
#' prediction <- list(metaData = metaData, data_matrix = pred)
#' #run the function
#' getSE(data = prediction, cores = 2)
#' @export
#' @importFrom GenomicRanges start reduce mcols findOverlaps width seqnames
#' @import parallel
#' @importFrom S4Vectors queryHits subjectHits

getSE <- function(data, cutoff = 0.5, distance = 12500, cores = 1) {

    data_matrix <- data$data_matrix

    #check if the parameters are valid
    if (cutoff < 0 | cutoff > 1) {
        message <- paste0(cutoff, " is not a valid cutoff. Please choose a cutoff between 0 and 1.")
        stop(message)
    }

    if (distance < 0) {
        message <- paste0(distance, " is not a valid distance. Please choose a distance greater than 0.")
        stop(message)
    }

    startPart("Call enhancer peaks and cluster")

    # single enhancer peaks
    cat(paste0(skip(), "define single enhancers"))
    peaks <- get_enhancerPeaks(data_matrix, cutoff, cores)

    cat(paste0(skip(), length(peaks), " single enhancer peak(s) found"))
    done()

    if (length(peaks) > 0) {
        # enhancer cluster
        cat(paste0(skip(), "define cluster of enhancer peaks"))
        cluster <- get_enhancerCluster(peaks, distance, cores)
        cat(paste0(skip(), length(cluster), " enhancer cluster found"))
        #GenomicRanges::start(peaks) <- GenomicRanges::start(peaks) - 1
        done()
    } else {
        cluster <- NULL
    }

    if ((!is.null(cluster)) & (length(cluster) > 0)) {
        #GenomicRanges::start(cluster) <- GenomicRanges::start(cluster) - 1
    }

    endPart()

    output <- list(metaData = data$metaData, data_matrix = data$data_matrix, peaks = peaks,
        cluster = cluster)
    return(output)
}
