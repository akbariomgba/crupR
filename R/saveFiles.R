#' Function to save the different output objects of each step
#'
#' @param data An output list of one of the steps 2 to 4 or
#' the get_superEnhancer step
#' @param modes Formats in which the GRanges object should be saved.
#' Following modes are available:
#' for the output of getEnhancers(): 'bigWig' for a bigWig file, 'rds'
#' for an .rds file
#' for the output of getDynamics(): 'beds' for saving each cluster in
#' a seperate bed file
#' for the output of getTargets(): 'UCSC' for a UCSC interaction file
#' for the output of getSE(): 'bedGraph' for saving the peak calls in a
#' bedGraph file, 'bed' for saving the clusters of peaks in a bed file
#' @param outdir Output directory in which the files should be saved
#' @param nearest Only relevant, if you want to save the output of
#' enhancerTargets. Specifies if the output was produced by the nearest
#' gene mode of the function or not. Default is false.
#' @return Nothing
#' @examples
#' \dontrun{
#' saveFiles(data = prediction_1_1, modes = c('bigWig', 'rds'),
#'             outdir = '/example/dir/')
#' saveFiles(data = dynamics, modes = 'beds',
#'             outdir ='/example/dir/')
#' saveFiles(data = targets, modes = 'UCSC',
#'             outdir = '/example/dir/', nearest = TRUE)}
#' @export
#' @importFrom rtracklayer export
#' @importFrom GenomicRanges mcols seqnames
#' @importFrom utils write.table


saveFiles <- function(data, modes, outdir, nearest = FALSE) {

    if (!dir.exists(outdir)) {
        #check the directory
        message <- paste0("Directory ", outdir, " doesn't exist!")
        stop(message)
    }

    #after the prediction step
    modes.pred <- c("bigWig", "bedGraph", "rds", "bed")
    #usethis::use_package('rtracklayer')

    if ("bedGraph" %in% modes) {
        out.bedgraph <- paste0(outdir, "singleEnh.bedGraph")
        rtracklayer::export(data$peaks, out.bedgraph)
    }
    if ("bed" %in% modes) {
        out.bed <- paste0(outdir, "clusterEnh.bed")
        rtracklayer::export(data$cluster, out.bed)
    }
    if ("bigWig" %in% modes) {
        out.bw <- paste0(outdir, "prediction.bw")
        #colnames(GenomicRanges::elementMetadata(object)) <- 'score'
        #GenomeInfoDb::seqlengths(object) <-
        #GenomicRanges::end(GenomicRanges::reduce(object))
        rtracklayer::export(data$data_matrix, out.bw)
    }
    if ("rds" %in% modes) {
        out.rds <- paste0(outdir, "prediction.rds")
        saveRDS(data$data_matrix, out.rds)
    }


    #after the dynamics step

    if ("beds" %in% modes) {
        peaks <- data$sumFile
        clusters <- unique(peaks$cluster)
        GR_header_short <- c("seqnames", "start", "end")

        for (c in clusters) {
            if (!grepl("r", c)) {
                out.bed <- paste0(outdir, paste0("dynamicEnh__cluster_", c, ".bed"))
                write.table(data.frame(peaks)[which(peaks$cluster == c), GR_header_short],
                    file = out.bed, quote = FALSE, row.names = FALSE, col.names = FALSE,
                    sep = "\t")
            }
        }
    }

    #after the targets step
    if ("UCSC" %in% modes) {
        #usethis::use_package('GenomicRanges')

        out.interaction <- ""
        units <- data$RegulatoryUnits
        header <- "track type=interact name=\"Dynamic Promoter-Enhancer Pairs\" description=\"Dynamic Promoter-Enhancer Pairs\" interactDirectional=true visibility=full"
        enhancer <- as.matrix(data.frame(units)[, c("start", "end")])

        promoter <- c()
        score <- 0
        if (nearest == FALSE) {
            out.interaction <- paste0(outdir, paste0("RegulatoryUnits.interaction"))
            promoter <- as.matrix(data.frame(GenomicRanges::mcols(units)$CORRELATED_GENE_PROMOTER_START,
                GenomicRanges::mcols(units)$CORRELATED_GENE_PROMOTER_END))
            score <- GenomicRanges::mcols(units)$CORRELATION
        } else {
            out.interaction <- paste0(outdir, paste0("RegulatoryUnitsNearestGene.interaction"))
            promoter <- as.matrix(data.frame(GenomicRanges::mcols(units)$NEAREST_GENE_PROMOTER_START,
                GenomicRanges::mcols(units)$NEAREST_GENE_PROMOTER_END))
            score <- GenomicRanges::mcols(units)$DISTANCE_TO_NEAREST
        }

        interaction <- cbind(as.character(GenomicRanges::seqnames(units)), apply(cbind(enhancer,
            promoter), 1, min), apply(cbind(enhancer, promoter), 1, max), rep(".",
            length(units)), rep(0, length(units)), score, rep(".", length(units)),
            rep("#7A67EE", length(units)), as.character(GenomicRanges::seqnames(units)),
            enhancer, rep(".", length(units)), rep(".", length(units)), as.character(GenomicRanges::seqnames(units)),
            promoter, rep(".", length(units)), rep(".", length(units)))

        write.table(header, file = out.interaction, quote = FALSE, col.names = FALSE,
            row.names = FALSE, sep = "\t")
        write.table(interaction, file = out.interaction, quote = FALSE, row.names = FALSE,
            col.names = FALSE, sep = "\t", append = TRUE)
    }

}
