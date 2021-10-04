#' Normalize the ChIPseq counts.
#'
#' @param metaData A data frame containing all important information about
#' the ChIPseq experiments, i,e. for which histone modification they were
#' conducted, path to the file, which condition or replicate, etc.
#' @param condition The number of the condition to normalize.
#' @param replicate The number of the replicate to normalize.
#' @param genome A string, the name of the genome that was used while mapping.
#' Possible options: mm9, mm10, hg19, hg38.
#' @param mapq A number, Minimum mapping quality of the reads that should be
#' included. Default is 10.
#' @param sequencing A string, what type of sequencing was used.
#' Possible options: paired and single
#' @param input.free A boolean, Whether or not to run the function
#' in the inputfree mode, in case there is no input experiment available.
#' Default is False.
#' @param chroms A vector of strings. Specify the relevant chromosomes
#' if needed. Then only the reads on these chromosomes are considered for
#' normalization and used in further analysis. Default is 'all'.
#' @param cores Number of cores to use. Default is 1.
#' @return A list containing the meta data of the experiments that were normalized and the normalized counts in a GRanges object
#' @examples
#'
#' files <- c(system.file('extdata', 'Condition1.H3K4me1.bam', package='crupR'),
#'         system.file('extdata', 'Condition1.H3K4me3.bam', package='crupR'),
#'         system.file('extdata', 'Condition1.H3K27ac.bam', package='crupR'),
#'         system.file('extdata', 'Condition2.H3K4me1.bam', package='crupR'),
#'         system.file('extdata', 'Condition2.H3K4me3.bam', package='crupR'),
#'         system.file('extdata', 'Condition2.H3K27ac.bam', package='crupR'))
#'
#' inp <- c(rep(system.file('extdata','Condition1.Input.bam',package='crupR'),3),
#'         rep(system.file('extdata','Condition2.Input.bam',package='crupR'),3))
#'
#' metaData <- data.frame(HM = rep(c('H3K4me1','H3K4me3','H3K27ac'),2),
#'         condition = c(1,1,1,2,2,2), replicate = c(1,1,1,1,1,1),
#'         bamFile = files, inputFile = inp)
#'
#' normalize(metaData = metaData, condition = 1, replicate = 1,
#'         genome = 'mm10', sequencing = 'paired', cores = 2)
#'
#' @export
#' @importFrom bamsignals bamProfile
#' @importFrom Rsamtools scanBamHeader BamFile
#' @importFrom stats knots
#' @importFrom GenomicRanges seqinfo mcols tileGenome width
#' @importFrom GenomeInfoDb seqlevelsStyle seqlevels
#' @importFrom GenomicAlignments summarizeOverlaps
#' @importFrom BSgenome.Mmusculus.UCSC.mm10 BSgenome.Mmusculus.UCSC.mm10
#' @importFrom BSgenome.Mmusculus.UCSC.mm9 BSgenome.Mmusculus.UCSC.mm9
#' @importFrom BSgenome.Hsapiens.UCSC.hg19 BSgenome.Hsapiens.UCSC.hg19
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 BSgenome.Hsapiens.UCSC.hg38

normalize <-
    function(metaData,
            condition,
            replicate,
            genome,
            mapq = 10,
            sequencing,
            input.free = FALSE,
            chroms = "all",
            cores = 1) {
    ##################################################################
    ##start  run time
    ##################################################################


    start_time <- Sys.time()

    ##################################################################
    ##check input parameter
    ##################################################################

    # check sequencing parameter
    if (!(sequencing %in% c("paired", "single"))) {
        mes <-
        "Sequencing parameter is not valid.\n Choose one of: paired , single."
        stop(mes)
    }

    # check summary file
    m <-
        metaData[which(metaData$condition == condition &
                        metaData$replicate == replicate),]  #get relevant rows

    if (nrow(m) != 3) {
        if (nrow(m) == 0) {
            message <-
            paste0(
            "The chosen combination of condition and replicate is not valid.\n There are no files for condition ",
            condition,
            " replicate ",
            replicate)
        } else if (nrow(m) > 3) {
            message <-
            paste0("There are too many bamfiles for condition ",
                    condition,
                    " replicate ",
                    replicate)
        } else {
            # there are only 1 or 2 files
            message <-
                paste0("There are not enough bamfiles for condition ",
                    condition,
                    " replicate ",
                    replicate)
        }
        stop(message)
    }

    for (path in m$bamFile) {
        check_file(path)
    }

    if (!(input.free)) {
        for (path in m$inputFile) {
        check_file(path)
        }
    }

    # check genome
    if (!(genome %in% genome_values)) {
        message <-
            paste0("Genome ",
            genome,
            " currently not provided. Choose one of: hg19 , hg38 , mm10 , mm9.")
        stop(message)
    }

    ##################################################################
    ##List input parameter
    ##################################################################

    startPart("List input parameter")

    # input <- c(as.vector(m$bamFile)) #,as.vector(m$inputFile))
    input <-
        c(as.vector(m[which(m$HM == "H3K4me1"),]$bamFile)[1],
            as.vector(m[which(m$HM ==
                        "H3K4me3"),]$bamFile)[1],
            as.vector(m[which(m$HM == "H3K27ac"),]$bamFile)[1])
    if (input.free == FALSE) {
        input <-
            c(
                input,
                as.vector(m[which(m$HM == "H3K4me1"),]$inputFile)[1],
                as.vector(m[which(m$HM == "H3K4me3"),]$inputFile)[1],
                as.vector(m[which(m$HM ==
                    "H3K27ac"),]$inputFile)[1]
        )
    }

    for (i in seq_along(input)) {
        input[i] <- normalizePath(input[i])
    }
    cat(paste0(skip(), "all ChIPseq files exist:\n"))
    cat(skip(), "H3K4me1 experiment: ", input[1], "\n")
    cat(skip(), "H3K4me3 experiment: ", input[2], "\n")
    cat(skip(), "H3K27ac experiment: ", input[3], "\n")
    if (input.free == FALSE) {
        cat(skip(), "Input experiment(s): ", input[4:length(input)], "\n")
    } else {
        cat(skip(), "Running the input free mode", "\n")
    }
    cat(skip(), "mapq: ", mapq, "\n")
    cat(skip(), "genome: ", genome, "\n")
    cat(skip(), "sequencing: ", sequencing, "\n")
    cat(skip(), "cores: ", cores, "\n")
    endPart()

    ##################################################################
    ## get the txdb file
    ##################################################################

    if (genome == "mm10") {
        txdb <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
    }
    if (genome == "mm9") {
        txdb <- BSgenome.Mmusculus.UCSC.mm9::BSgenome.Mmusculus.UCSC.mm9
    }
    if (genome == "hg19") {
        txdb <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    }
    if (genome == "hg38") {
        txdb <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    }

    ##################################################################
    ##get summarized informationof ChiP-seq experiments
    ##################################################################

    startPart("Get information of ChIP-seq experiments")

    # get and check feature set
    features.valid <- c("H3K4me1", "H3K4me3", "H3K27ac")
    bam_files_HM <- input[seq_len(3)]
    if (input.free == TRUE) {
        bam_files_input <- NULL
    } else {
        bam_files_input <- unique(input[4:length(input)])
    }

    endPart()

    ##################################################################
    ##get and adjust binned genome
    ##################################################################

    startPart("Prepare the binned genome")

    # bin the genome
    cat(paste0(skip(), "get binned genome for ", genome))
    binned_genome <- get_binned_genome(txdb, tilewidth <- 100)
    done()

    if (chroms != "all") {
        p <-
            which(as.character(GenomicRanges::seqnames(binned_genome)) %in% chroms)
        binned_genome <- binned_genome[p]
    }


    # get prefix of chromosome names
    cat(paste0(skip(), "get prefix of chromosomes from  file ", bam_files_HM[1]))
    bam_header <-
        Rsamtools::scanBamHeader(bam_files_HM[1])[[1]]$targets
    done()


    # gives information about path of bam and bai file (Object)
    bf <- Rsamtools::BamFile(bam_files_HM[1])

    # gives length of different chromosomes (SeqInfo class)
    si <- GenomicRanges::seqinfo(bf)

    # check if prefix of chromosome names ('chr1' or '1')
    cat(paste0(skip(), "adjust prefix of chromosome names in binned genome"))
    if (!("chr1" %in% names(bam_header))) {
        GenomeInfoDb::seqlevelsStyle(binned_genome) <-
            GenomeInfoDb::seqlevelsStyle(si)[1]
    }


    done()
    endPart()

    ##################################################################
    ##get counts from ChIP-seq experiments
    ##################################################################

    startPart("Get summarized counts from ChIP-seq experiments")

    # prepare feature names:
    names <- c(features.valid, paste0("Input_", features.valid))

    if (input.free == TRUE)
        names <- features.valid

    if (length(unique(bam_files_input)) == 1) {
        bam_files_input <- bam_files_input[1]
        names <- c(features.valid, "Input_All")
    }

    # get counts for ChIP-seq HM experiments
    bams <- c(bam_files_HM, bam_files_input)


    if (sequencing == "paired") {
        counts <-
            parallel::mcmapply(
            bamsignals::bamProfile,
            bampath = bams,
            MoreArgs = list(
                gr = binned_genome,
                binsize = 100,
                mapqual = mapq,
                ss = FALSE,
                paired.end = "midpoint",
                filteredFlag = 1024,
                verbose = FALSE
            ),
            mc.cores = cores,
            SIMPLIFY = FALSE
        )

    } else if (sequencing == "single") {
        counts <-
            parallel::mcmapply(
            bamsignals::bamProfile,
            bampath = bams,
            MoreArgs = list(
                gr = binned_genome,
                binsize = 100,
                mapqual = mapq,
                ss = FALSE,
                shift = 100,
                filteredFlag = 1024,
                verbose = FALSE
            ),
            mc.cores = cores,
            SIMPLIFY = FALSE
        )
    }

    if (length(counts) == 4) {
        counts <-
            list(unlist(as.list(counts[[1]])),
                unlist(as.list(counts[[2]])),
                unlist(as.list(counts[[3]])),
                unlist(as.list(counts[[4]])))
    } else if (input.free == TRUE) {
        counts <-
            list(unlist(as.list(counts[[1]])), unlist(as.list(counts[[2]])),
                unlist(as.list(counts[[3]])))
    } else {
        # then it should have length 6
        counts <-
            list(
                unlist(as.list(counts[[1]])),
                unlist(as.list(counts[[2]])),
                unlist(as.list(counts[[3]])),
                unlist(as.list(counts[[4]])),
                unlist(as.list(counts[[5]])),
                unlist(as.list(counts[[6]]))
            )
    }

    names(counts) <- names
    endPart()

    ##################################################################
    ##get counts from ChIP-seq experiments
    ##################################################################

    if (input.free == FALSE) {
        startPart("Normalize histone modifications by Input")

        if ("Input_All" %in% names(counts)) {
            normalized_counts <-
                lapply(features.valid, function(x)
                log2((counts[[x]] +
                    1) /
                (counts[[paste0("Input_All")]] + 1)))
        } else {
            normalized_counts <-
                lapply(features.valid, function(x)
                log2((counts[[x]] +
                    1) /
                (counts[[paste0("Input_", x)]] + 1)))
        }
        endPart()
    } else {
        normalized_counts <- counts
    }


    ##################################################################
    ##create data matrix for all normalized ChIP-seq experiments
    ##################################################################

    startPart("Create summarized data matrix")

    cat(paste0(skip(), "create matrix"))
    GenomicRanges::mcols(binned_genome) <-
        matrix(
            unlist(normalized_counts),
            ncol = 3,
            byrow = FALSE,
            dimnames = list(NULL, features.valid))
    GenomeInfoDb::seqlevels(binned_genome) <-
        paste0("chr", gsub("chr|Chr", "", GenomeInfoDb::seqlevels(binned_genome)))
    done()

    # include log2 H3K4me1/H3K4me3 ratio
    cat(paste0(skip(), "include log2 H3K4me1/H3K4me3 ratio"))
    nominator <-
        GenomicRanges::mcols(binned_genome)[, "H3K4me1"] + abs(min(GenomicRanges::mcols(binned_genome)[,
                                                                                                    "H3K4me1"])) + 1
    denominator <-
        GenomicRanges::mcols(binned_genome)[, "H3K4me3"] + abs(min(GenomicRanges::mcols(binned_genome)[,
                                                                                                    "H3K4me3"])) + 1
    GenomicRanges::mcols(binned_genome)[, "ratio"] <-
        log2(nominator / denominator)
    done()

    done()

    endPart()

    ##################################################################
    ##print run time
    ##################################################################

    run_time <- Sys.time() - start_time

    startPart("Run time")
    cat(paste0(skip(), format(run_time), "\n"))
    endPart()
    output <- list(metaData = m, normalized = binned_genome)
    return(output)
}
