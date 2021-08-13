#' Find the target genes of the regulatory, condition-specific clusters
#'
#' @param data Outputfile of the previous step (getDynamics())
#' @param expression GRanges object containing the gene expression counts of RNAseq experiments for each condition and its replicates
#' @param genome Genome used in the .bam files of the RNAseq experiments. Possible options are 'mm9', 'mm10', 'hg19' and 'hg38'.
#' @param TAD.file Path to the TAD file to use for finding the target genes. If set to 'default', the default file is used (only if the 'mm10' genome was used)
#' @param threshold Threshold for correlation between cluster and gene. Default is 0.9.
#' @param nearest Logical: if set, the nearest gene is taken to build the regulatory regions.
#' @param cores Number of cores to use
#' @return A list containing the meta data of the experiments
#' and a GRanges object containing the dynamic regulatory units
#' @examples
#' 
#' #first get the output of crupR::getDynamics so skip this
#' dynamics <- readRDS(system.file("extdata", "differential_enhancers.rds", package="crupR"))
#' #load your GRanges object containing the gene expressions counts and run the function
#' expr <- readRDS(system.file("extdata", "expressions.rds", package="crupR"))
#' getTargets(data <- dynamics, expression <- expr, genome <- "mm10", cores = 2)
#'
#' @export
#' @importFrom GenomicRanges mcols GRanges makeGRangesFromDataFrame start end promoters strand seqnames nearest distance
#' @importFrom GenomicFeatures genes
#' @import parallel
#' @importFrom GenomeInfoDb seqlevels seqlengths
#' @importFrom IRanges IRanges subsetByOverlaps %within%
#' @import BSgenome
#' @import TxDb.Mmusculus.UCSC.mm9.knownGene
#' @import TxDb.Mmusculus.UCSC.mm10.knownGene
#' @import TxDb.Hsapiens.UCSC.hg19.knownGene
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @importFrom utils read.table

getTargets <- function(data, expression, genome, TAD.file = 'default', threshold = 0.9, nearest = FALSE, cores = 1){

  ##################################################################
  # start run time
  ##################################################################

  start_time <- Sys.time()

  ##################################################################
  # define fixed parameters;
  ##################################################################

  ID_prefix <- "cond"
  GR_header <- c("seqnames", "start","end","width")
  GR_header_short <- c("seqnames", "start","end")
  DF_header <- c("chr", "start","end")

  ##################################################################
  # check input parameter
  ##################################################################

  genome_values <- c('hg19', 'mm10', 'mm9', 'hg38')
  # check genome
  if (!(genome %in% genome_values)) {
    message <- paste0("Genome " , genome, " currently not provided. Choose one of: hg19 , hg38 , mm10 , mm9.");
    stop(message);
    }

  # check values
  if(threshold < 0.5 | threshold > 1) {
    message <- paste0(threshold, " is not in range [0.5,1].")
    stop(message);
  }


  if (nearest == FALSE){
      if (TAD.file == 'default' & genome == "mm10") {
          TAD   <- system.file("extdata", "mESC_mapq30_KR_all_TADs.bed", package = "crupR")
      }
      else if (TAD.file == 'default')  {
          message <- paste0("You have to provide your own file with TAD domains (fitting to the genome choice).")
          stop(message);
      } else { #TAD.file != 'defualt'
          TAD <- TAD.file
      }
  }

  ##################################################################
  # define input parameter
  ##################################################################

  startPart("List input parameter")

  if (nearest == TRUE){
      cat(paste0("Will choose the nearest gene to a each differential region to build a regulatory unit.\n"))
  }

  cat(skip(), "cores: ",cores, "\n")

  endPart()

  ##################################################################
  # libraries
  ##################################################################

  startPart("Load packages")

  txdb <- ""

  if (genome == "mm10") txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
  if (genome == "mm9") txdb  <- TxDb.Mmusculus.UCSC.mm9.knownGene::TxDb.Mmusculus.UCSC.mm9.knownGene
  if (genome == "hg19") txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  if (genome == "hg38") txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene


  endPart()

  ##################################################################
  # get normalized gene expression values
  ##################################################################

    startPart("Get normalized gene expression values from saved rds file")
    expr.gr <- expression
    IDs <- colnames(GenomicRanges::mcols(expr.gr)[,-1])
    endPart()

  ##################################################################
  # get correlation:
  ##################################################################

  startPart("Build Regulatory Units")

  # dynamic enhancer regions:
  regions.gr <- data$sumFile
  GenomeInfoDb::seqlevels(regions.gr) <- paste0("chr", gsub("chr|Chr","", GenomeInfoDb::seqlevels(regions.gr)))
  #cluster.U <- which(GenomicRanges::mcols(regions.gr)$cluster == 'U')
  #if(length(cluster.U) > 0) {regions.gr <- regions.gr[-cluster.U]}

  # TAD/domain regions:
  TAD.gr <- GenomicRanges::GRanges()

  if (nearest == FALSE){
      TAD.df <- read.table(TAD, col.names = GR_header_short)
      TAD.df <- TAD.df[which((TAD.df$end-TAD.df$start) > 0),]
      TAD.gr <- GenomicRanges::makeGRangesFromDataFrame(TAD.df)
      GenomeInfoDb::seqlevels(TAD.gr) = paste0("chr", gsub("chr|Chr","", GenomeInfoDb::seqlevels(TAD.gr)))
      units <- get_units(regions.gr,
                          expr.gr,
                          TAD.gr,
                          IDs,
                          cores,
                          threshold,
                          txdb)
  }else{
      units <- get_units(regions.gr,
                          expr.gr,
                          TAD.gr,
                          IDs,
                          cores,
                          threshold,
                          txdb,
                          nearest = TRUE)
  }


  out.txt <- ""

  cat(paste0(skip(), "Prepare output file"))

  output <- list("metaData" = data$metaData,  "RegulatoryUnits" = units)
  done()
  endPart()

  ##################################################################
  # print run time
  ##################################################################

  run_time <- Sys.time() - start_time

  startPart("Run time")
  cat(paste0(skip(), format(run_time), "\n"))
  endPart()

  return(output)
}
