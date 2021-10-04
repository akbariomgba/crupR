
##################################################################
# definition: valid parameter values:
##################################################################

sequencing_values <- c("paired", "single")
genome_values <- c("hg19", "mm10", "mm9", "hg38")

##################################################################
#definition: output messages
##################################################################

done <- function() {
    cat(".. done.\n")
}
skip <- function() {
    cat("\t ..")
}
startPart <- function(m) {
    cat(paste0("\n--- ", m, " ---\n\n"))
}
endPart <- function() {
    cat("\n\t>>> All done!\n")
}

##################################################################
#definition: standard GRanges/DataFrame header:
##################################################################

GR_header <- c("seqnames", "start", "end", "width")
GR_header_short <- c("seqnames", "start", "end")
DF_header <- c("chr", "start", "end")

##################################################################
# function: re-header
##################################################################

reheader_DF <- function(DF, header) {
    colnames(DF)[seq_along(header)] <- header
    return(DF)
}

##################################################################
# function: check if file exists
##################################################################

check_file <- function(f) {
    if (!(file.exists(f))) {
        message <- paste0("File ", f, " does not exist.\n")
        stop(message)
    }
}

##################################################################
# function: check if outdir exists
#
#(d = outdir, alt_d = alternative dir)
##################################################################

check_outdir <- function(d, alt_d) {
    if (is.null(d))
        d <- dirname(alt_d)
    if (!dir.exists(d)) {
        cat(paste0("Output directory '", d, "'is not a valid directory. \n
            Output directory is set to ",
            dirname(alt_d)))
        d <- paste0(dirname(alt_d), "/")
    }
    return(d)
}

##################################################################
# function: partition genome into bins
##################################################################

get_binned_genome <- function(txdb, width) {

    # get binned genome from txdb object
    binned <- GenomicRanges::tileGenome(GenomeInfoDb::seqinfo(txdb), tilewidth = width,
        cut.last.tile.in.chrom = TRUE)

    # only take autosomes and X chromosome
    GenomeInfoDb::seqlevels(binned, pruning.mode = "coarse") <- GenomeInfoDb::seqlevels(txdb)[grep("^chr[0-9]{,2}$|chrX$",
        GenomeInfoDb::seqlevels(txdb))]

    # delete last region in each chromosome that is smaller than 100
    return(binned[-which(GenomicRanges::width(binned) != 100)])
}

##################################################################
# function: get summarized counts in defined bins --> kann eigentlich weg
##################################################################

#get_bamProfile <- function(bampath, gr, mapq, sequencing){

# gives information about path of bam and bai file (Object)
#  bf <- Rsamtools::BamFile(bampath)

# gives length of different chromosomes (SeqInfo class)
#  si <- GenomicRanges::seqinfo(bf)

# fix chromosome prefix
# seqlevels(gr) <- gsub('chr','',seqlevels(gr))
#  if (grepl('chr', seqlevels(si)[1])) {
#    seqlevels(gr) <- paste0('chr', seqlevels(gr))
#  }

# filter chromosomes
#filter <- IRanges::'%in%'(seqnames(gr), seqnames(si))
# filter seqnames(gr) in seqnames(si)
#  gr <- gr[filter]

# count and summarize
#  if (sequencing == 'paired') {
#    sapply(bamsignals::bamProfile( bampath = bampath,
#                        gr = gr,
#                        binsize = 100,
#                        mapqual = mapq,,
#                        paired.end = 'midpoint',
#                        filteredFlag = 1024,
#                        verbose = FALSE),
#            function(x) x)
#  } else if (sequencing == 'single') {
#    sapply( bamsignals::bamProfile( bampath = bampath,
#                        gr = gr,
#                        binsize = 100,
#                        mapqual = mapq,
#                        ss = FALSE,
#                        shift = 100,
#                        filteredFlag = 1024,
#                        verbose = FALSE),
#            function(x) x)
#  }
#}

##################################################################
# function: quantile normalize with target
##################################################################

get_targetQuantileNorm <- function(ecdf) {

    ### recreate HM count vector from ecdf of mESC
    x.ref <- knots(ecdf)
    y.ref <- ecdf(x.ref)

    # 26337756 = nb of 100 bp bins in mm10
    temp <- c(0, round(y.ref * 26337756, digits = 0))
    ret <- unlist(lapply(seq_along(x.ref), function(x) rep(x.ref[x], temp[x +
        1] - temp[x])))

    return(ret)
}


##################################################################
# function: create extended data matrix (plus/minus) bins
##################################################################

extend_dataMatrix <- function(N, df, f) {

    N_regions <- nrow(df)

    # make extended feature vector
    f_ext <- NULL
    for (i in seq_len(N)) {
        f_ext <- c(f_ext, paste0(f, "_left", i))
        f_ext <- c(f_ext, paste0(f, "_right", i))
    }

    # prepare new matrix
    df_ext <- cbind(df[, c(GR_header_short, f)], matrix(0, nrow = N_regions,
        ncol = length(f_ext)))
    colnames(df_ext) <- c(DF_header, f, f_ext)

    # make extended data matrix
    region_ext <- (N + 1):(N_regions - N)

    for (i in seq_len(N)) {
        df_ext[region_ext, f_ext[((2 * i - 2) * length(f) + 1):((2 * i -
            1) * length(f))]] <- df[region_ext - i, f]
        df_ext[region_ext, f_ext[((2 * i - 1) * length(f) + 1):(2 * i *
            length(f))]] <- df[region_ext + i, f]
    }

    return(df_ext)
}


#############################################################################
# function: sort peak candidates and remove overlapping peaks
# used in peak_calling function
#############################################################################

sort_peaks <- function(peaks) {
    # sort peaks according to probability score
    peaks <- peaks[sort(GenomicRanges::mcols(peaks)$score, decreasing = TRUE, index.return = TRUE)$ix]
    count <- 0
    while (length(peaks) > (count + 1)) {

        count <- count + 1
        overlap.to <- S4Vectors::subjectHits(GenomicRanges::findOverlaps(query = peaks[count],
            subject = peaks))

        if (length(overlap.to) == 1)
            next

        delete.index <- sort(overlap.to, decreasing = FALSE)[-1]
        peaks <- peaks[-delete.index]
    }
    return(peaks)
}


#############################################################################
# function: call peaks (from probabilities)
#############################################################################

get_enhancerPeaks <- function(gr, cutoff, cores) {

    # define and merge all regions under background probability
    p <- which(GenomicRanges::mcols(gr)$score > cutoff)
    candidates <- GenomicRanges::reduce(gr[p])

    # possible 100 bp peaks
    peaks <- gr[S4Vectors::queryHits(GenomicRanges::findOverlaps(gr,
        candidates))]

    # create 1100 bp peak regions
    GenomicRanges::start(peaks) <- GenomicRanges::start(peaks) - 500
    GenomicRanges::width(peaks) <- 1100

    # call peaks for each chromosome separately
    out <- parallel::mclapply(split(peaks, GenomicRanges::seqnames(peaks)),
        sort_peaks, mc.cores = cores)

    return(do.call("c", unname(out)))
}


##################################################################
# function: call super enhancer candidates from called peaks
##################################################################

get_enhancerCluster <- function(peaks, peak.gap, cores) {

    # cluster candidates (neighboring enhancers within 12.5 kb)
    red.peaks <- GenomicRanges::reduce(peaks, min.gapwidth = peak.gap)
    cluster <- red.peaks[which(GenomicRanges::width(red.peaks) > 1100)]

    # sort cluster accroding to maximum probability
    if (length(cluster) > 0) {
        sort.max <- unlist(parallel::mclapply(seq_along(cluster), function(x) max(peaks$score[S4Vectors::subjectHits(GenomicRanges::findOverlaps(cluster[x],
            peaks))]), mc.cores = cores))
        cluster <- cluster[sort(sort.max, index.return = TRUE, decreasing = TRUE)$ix]
    }
    return(cluster)
}


##################################################################
# function: compute Kolmogorow-Smirnow test statistic
#
#
##################################################################

get_KS.STATISTIC <- function(d, W) {
    #parameters
    # for KS TEST:
    N <- W * 2 + 1
    Nplus <- rep(1/N, N * 2)
    Nminus <- -1/N

    d.sort <- sort.int(d, index.return = TRUE, method = "quick")
    z <- cumsum(replace(Nplus, d.sort$ix <= N, Nminus))
    return(max(abs(z[c(diff(d.sort$x) != 0, TRUE)])))
}


##################################################################
# function: get indices of regions
#whose mean difference is higher than w_0
##################################################################

get_idx <- function(c1m, c2m, w_0, W) {

    z <- rle(abs(c1m - c2m) < w_0) %>%
        unclass() %>%
        as.data.frame() %>%
        dplyr::mutate(end = cumsum(lengths), start = c(1, dplyr::lag(end)[-1] +
            1)) %>%
        magrittr::extract(c(1, 2, 4, 3)) %>%
        dplyr::filter(values == TRUE) %>%
        dplyr::filter(lengths >= (W * 2 + 1))

    return(dplyr::setdiff(seq(length(c1m)), unlist(apply(z[, c("start", "end")],
        1, function(x) seq(x[1], (x[2] - W + 1))))))
}

##################################################################
# function: get pairwise p values
##################################################################
# get pairwise pvalues:
get_pairwisePvalues <- function(p, I, w_0, W, p.thres, C) {
    #parameters
    N <- W * 2 + 1
    KS.FACTOR <- sqrt(N * 0.5)

    comb <- combn(seq(length(I)), 2)
    ret <- list()

    for (i in seq(ncol(comb))) {
        p.value <- rep(NA, length(p))

        # define IDs:
        i1 <- unlist(I[comb[1, i]])
        i2 <- unlist(I[comb[2, i]])

        # mean weighted by sd:
        if (length(i2) > 1) {
            c1m <- rowMeans(as.matrix(GenomicRanges::mcols(p)[, i1]))
            c1s <- (1 - sqrt(matrixStats::rowVars(as.matrix(GenomicRanges::mcols(p)[,
                i1]))) + 1e-07)
        } else {
            c1m <- GenomicRanges::mcols(p)[, i1]
            c1s <- rep(1, length(c1m))
        }
        if (length(i2) > 1) {
            c2m <- rowMeans(as.matrix(GenomicRanges::mcols(p)[, i2]))
            c2s <- (1 - sqrt(matrixStats::rowVars(as.matrix(GenomicRanges::mcols(p)[,
                i2]))) + 1e-07)
        } else {
            c2m <- GenomicRanges::mcols(p)[, i2]
            c2s <- rep(1, length(c2m))
        }

        # adjust mean for standard deviations
        z <- cbind(c1m/c1s, c2m/c2s)
        # get indices with minimal difference of w_0:
        idx <- get_idx(c1m, c2m, w_0, W)
        idx <- idx[(W + 1):(length(idx) - W)]

        if (length(idx) == 0)
            next

        # get means over replicates:
        D <- unlist(parallel::mclapply(idx, FUN = function(x) get_KS.STATISTIC(z[(x -
            W):(x + W), ], W), mc.cores = C, mc.allow.recursive = FALSE))
        C_pKS2 <- utils::getFromNamespace("C_pKS2", "stats")
        #p.value[idx] <- 1 - .Call( stats:::C_pKS2, p = KS.FACTOR * D, tol = 0.000001)
        p.value[idx] <- 1 - .Call(C_pKS2, p = KS.FACTOR * D, tol = 1e-06)
        p.adjust <- p.adjust(p.value, method = "bonferroni")

        idx.significant <- which(p.adjust <= p.thres)
        if (length(idx.significant) == 0) {
            message <- "No regions with significant p-values were found!"
            stop(message)
        }

        this.result <- GenomicRanges::granges(p)[idx.significant]

        # I think, we don't need all of this:
        if (length(idx != 0)) {

            GenomicRanges::mcols(this.result) <- data.frame(p.value = p.value[idx.significant],
                p.adj = p.adjust[idx.significant], p.direction = unlist(lapply(idx.significant,
                    function(x) (mean(z[(x - W):(x + W), 1]) - mean(z[(x - W):(x +
                    W), 2]) <= 0))), idx = idx.significant, comparison = paste0(unique(gsub("_.*",
                    "", unlist(I[comb[, i]]))), collapse = ","))
        }
        # 'direction' of probabilties:
        ret[[paste0(comb[, i], collapse = ",")]] <- this.result[which(this.result$p.direction ==
            TRUE)]
        ret[[paste0(rev(comb[, i]), collapse = ",")]] <- this.result[which(this.result$p.direction ==
            FALSE)]
        #the same
        ret[[paste0(comb[, i], collapse = ",")]] <- this.result[which(this.result$p.direction ==
            TRUE)]
        ret[[paste0(rev(comb[, i]), collapse = ",")]] <- this.result[which(this.result$p.direction ==
            FALSE)]
    }

    return(ret)
}

##################################################################
# function: get the pairwise pattern
##################################################################

get_cluster <- function(p, pvalues, I) {

    # specify all possible pairwise combinations:
    comb <- combn(seq(seq_along(I)), 2)
    comb.f <- cbind(comb, comb[c(2, 1), ])

    # create emtpy pattern matrix:
    idx <- unique(unlist(lapply(pvalues, function(x) x$idx)))
    m <- matrix(0, nrow = length(idx), ncol = ncol(comb.f))
    colnames(m) <- apply(comb.f, 2, function(x) paste(x, collapse = ","))
    rownames(m) <- idx

    # fill matrix:
    for (i in colnames(m)) m[as.character(pvalues[[i]]$idx), i] <- 1

    # add to probabilities:
    p$sign <- FALSE
    p$sign[idx] <- TRUE

    p$pattern <- paste0(rep(0, ncol(comb.f)), collapse = "")
    p$pattern[idx] <- apply(m, 1, function(x) paste0(x, collapse = ""))

    return(p)
}

##################################################################
# function: get summarized ranges
##################################################################

get_ranges <- function(p, I, W, C) {

    # get significant positions:
    p.sign <- p[which(p$sign == TRUE)]

    #extend to actual size:
    GenomicRanges::start(p.sign) <- GenomicRanges::start(p.sign) - W * 100
    GenomicRanges::end(p.sign) <- GenomicRanges::end(p.sign) + W * 100

    # split by cluster:
    p.split <- split(p.sign, p.sign$pattern)
    p.red <- lapply(p.split, GenomicRanges::reduce)
    pattern <- rep(names(p.red), unlist(lapply(p.red, length)))
    p.red <- do.call("c", unname(p.red))
    p.red$pattern <- pattern

    # get maximum probability value within each range for each condition:
    overlap <- GenomicRanges::findOverlaps(p, p.red)

    p$range <- NA
    p$range[S4Vectors::queryHits(overlap)] <- S4Vectors::subjectHits(overlap)
    p.split2 <- split(p, p$range)

    #FUN = function(x) matrixStats::colMaxs(as.matrix(GenomicRanges::mcols(x)[,unlist(I)])),
    #FUN = function(x) colMeans(as.matrix(mcols(x)[,unlist(I)])),

    GenomicRanges::mcols(p.red)[, unlist(I)] <- do.call("rbind", parallel::mclapply(p.split2,
        FUN = function(x) matrixStats::colMedians(as.matrix(GenomicRanges::mcols(x)[,
            unlist(I)])), mc.cores = C, mc.allow.recursive = FALSE))
    return(p.red)
}

##################################################################
# function: rename the clusters
##################################################################

name_patterns <- function(gr, C) {
    patterns <- unique(gr$pattern)
    sizes <- vapply(seq_along(patterns), function(i) length(which(gr$pattern ==
        patterns[i])), FUN.VALUE = numeric(1))
    frame <- data.frame(pattern = patterns, cluster_size = sizes)
    sorted.frame <- frame[order(-frame$cluster_size), ]
    hashtab <- list()
    patterns <- as.vector(sorted.frame$pattern)
    for (i in seq_along(patterns)) {
        hashtab[patterns[i]] <- paste0("c_", i)
    }  #vapply!
    clusters <- parallel::mclapply(gr$pattern, function(x) hashtab[[x]],
        mc.cores = C)
    gr$cluster <- clusters
    return(gr)
}

##################################################################
# function: adjust y positioning for heatmap (y axis)
##################################################################

adjust <- function(v) {

    mid <- v[1]
    res <- c(mid)
    for (i in 2:length(v)) {
        #seq?
        mid <- (mid + v[i - 1]/2 + v[i]/2)
        res <- c(res, mid)
    }
    return(res)
}

##################################################################
# function: plot summary of the K first differential enhancer regions
##################################################################
#' plots boxplots of the median enhancer prediction values of the enhancers in the condition-specific clusters
#'
#' @param dynamic_enhancers The getDynamics() output file of containing
#' the GRanges object with the differential enhaners
#' @param num_plots Maximal number of cluster whose plots should be
#' displayed (clusters are sorted by their sizes). This parameter can be
#' set in case the number of clusters is very high (20+). Default is 20.
#' @return a figure containing a boxplot for each cluster
#' @examples
#'
#' #get the output of crupR::getDynamics()
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
#' clusters <- readRDS(system.file('extdata', 'differential_enhancers.rds', package='crupR'))
#' dynamics <- list(metaData = metaData, sumFile = clusters)
#' plotSummary(dynamics)#plot the clusters
#'
#' @export
#' @importFrom ggplot2 ggplot aes facet_wrap geom_boxplot theme_bw guides xlab ylab
#' @importFrom reshape2 melt
#' @importFrom GenomicRanges mcols


plotSummary <- function(dynamic_enhancers, num_plots = 20) {
    p <- dynamic_enhancers$sumFile

    if (length(unique(p$cluster)) > num_plots) {
        valid_clusters <- paste0("c_", seq_len(num_plots))
        p <- p[which(p$cluster %in% valid_clusters)]
    }
    clusters <- p$cluster
    p$cluster <- NULL
    d <- suppressMessages(reshape2::melt(data.frame(GenomicRanges::mcols(p))))
    d$condition <- gsub("_.*", "", d$variable)
    d$replicate <- gsub(".*_", "", d$variable)
    tab <- table(p$pattern)
    d$label <- paste0(clusters, " (", d$pattern, "): ",
        tab[d$pattern], " regions")



    ggplot2::ggplot(d, ggplot2::aes(x = condition, y = value,
        col = condition, shape = replicate)) + ggplot2::facet_wrap(label ~
        ., ncol = round(length(unique(d$label))/4 + 0.3)) +
        ggplot2::geom_boxplot() + ggplot2::theme_bw() +
        ggplot2::guides(shape = FALSE, col = FALSE) + ggplot2::xlab("") +
        ggplot2::ylab("max. enh. prob. in range")
}

##################################################################
# function: get summarized counts in defined bins
##################################################################


#get_bamOverlaps <- function(files, IDs, txdb, singleEnd){
#usethis::use_package('GenomicAlignments')  # for summarizeOverlaps()
#usethis::use_package('Rsamtools')          # for seqinfo()
#usethis::use_package('DESeq2')
#usethis::use_package('GenomeInfoDb')
#usethis::use_package('GenomicFeatures')
# gives information about path of bam and bai file (Object)
#  bf = Rsamtools::BamFileList(unlist(files))

# gives length of different chromosomes (SeqInfo class)
#  si = GenomeInfoDb::seqinfo(bf)

# get exons and genes
#  expr.gr <- GenomicFeatures::genes(txdb)
#  GenomeInfoDb::seqlevels(expr.gr) = paste0('chr', gsub('chr','', GenomeInfoDb::seqlevels(expr.gr)))

#  exons <- GenomicFeatures::exonsBy(txdb, by = 'gene')
#  exons <- exons[which((names(exons) %in% mcols(expr.gr)$gene_id) == T)]

#  # fix chromosome prefix:
#  GenomeInfoDb::seqlevels(exons) <- gsub('chr','', GenomeInfoDb::seqlevels(exons))
#  if (grepl('chr', GenomeInfoDb::seqlevels(si)[1])) {
#    GenomeInfoDb::seqlevels(exons) <- paste0('chr', GenomeInfoDb::seqlevels(exons))
#  }

#  se <- GenomicAlignments::summarizeOverlaps(exons,
#                          bf,
#                            mode='Union',
#                         singleEnd = singleEnd,
#                          fragments = setdiff(c(FALSE,TRUE), singleEnd))

# get counts and summarize per gene:
#  counts.per.exon <- data.frame(assays(se)$counts,
#                                gene = rownames(assays(se)$counts))
#  counts.split <- split(counts.per.exon,
#                        counts.per.exon$gene)
#  counts.per.gene <- (do.call(rbind,lapply(counts.split, function(x) colSums(x[,-dim(x)[2]]))))

# create new symmarized experiment object:
#  se0 <- SummarizedExperiment( assays = SimpleList(counts = counts.per.gene),
#                               colData = names(bf)
#  )

# stabilize the variance across the mean:
# (The transformed data should be approximated variance stabilized
# and also includes correction for size factors or normalization factors)

#  dds <- suppressMessages(DESeqDataSet(se0, ~ 1))
#  vsd <- varianceStabilizingTransformation(dds,
#                                            blind = TRUE)

#  counts_vst <- assay(vsd)
#  counts_raw <- assay(dds)
#  expr_raw.gr <- expr.gr

#  for (i in 1:dim(counts_vst)[2]) {
#        mcols(expr.gr)[,unlist(IDs)[i]] <- counts_vst[,i]
#        mcols(expr_raw.gr)[,unlist(IDs)[i]] <- counts_raw[,i]
#  }

#  expr.gr <- expr.gr[which(rowVars(counts_vst) > 0)]
#  expr_raw.gr <- expr_raw.gr[which(rowVars(counts_vst) > 0)]

#  return(list(raw = expr_raw.gr, vst = expr.gr))
#}


##################################################################
# function: create TAD if not defined
##################################################################

check_TAD <- function(t, TAD, this.region, regions) {

    if (length(t) == 0) {

        precedeT <- GenomicRanges::precede(TAD, this.region)
        followT <- GenomicRanges::follow(TAD, this.region)

        if (length(which(!is.na(precedeT))) == 0) {
            start <- 1
        } else {
            start <- GenomicRanges::end(TAD[rev(which(!is.na(precedeT)))[1]])
            if (length(start) == 0)
                start <- 1
        }

        if (length(which(!is.na(followT))) == 0) {
            end <- max(GenomicRanges::end(regions[which(GenomicRanges::seqnames(regions) ==
                GenomicRanges::seqnames(this.region))]))
        } else {
            end <- GenomicRanges::start(TAD[(which(!is.na(followT)))[1]])
            if (length(end) == 0)
                end <- max(GenomicRanges::end(regions[which(GenomicRanges::seqnames(regions) ==
                    GenomicRanges::seqnames(this.region))]))
        }

        t <- GenomicRanges::GRanges(GenomicRanges::seqnames(this.region), IRanges::IRanges(start,
            width = (end - start + 1)))
    }
    return(t)
}

##################################################################
# function: correlate probabilities and gene exression counts
##################################################################

get_correlation <- function(i, threshold, regions.gr, expr.gr, TAD.gr, IDs) {
    interactions <- data.frame(stringsAsFactors = FALSE)
    # get region, associated TAD and genes within TAD
    this.region <- regions.gr[i]
    this.TAD <- IRanges::subsetByOverlaps(TAD.gr, this.region)
    this.TAD <- check_TAD(this.TAD, TAD.gr, this.region, regions.gr)
    this.genes.idx <- IRanges::"%within%"(expr.gr, this.TAD)
    if (sum(this.genes.idx) > 0) {

        this.genes <- expr.gr[this.genes.idx, ]
        cor <- apply(as.matrix(GenomicRanges::mcols(this.genes)[, IDs]), 1, function(x) cor(x,
            as.numeric(unlist(GenomicRanges::mcols(this.region)[, IDs]))))

        for (c in seq_along(cor)) {
            if (!is.na(cor[c]) && cor[c] >= threshold) {

                promoter_start <- GenomicRanges::start(GenomicRanges::promoters(this.genes[c]))
                promoter_end <- GenomicRanges::end(GenomicRanges::promoters(this.genes[c]))

                if (as.character(GenomicRanges::strand(this.genes[c])) == "-") {
                    promoter_start <- GenomicRanges::end(GenomicRanges::promoters(this.genes[c]))
                    promoter_end <- GenomicRanges::start(GenomicRanges::promoters(this.genes[c]))
                }

                if (length(this.TAD) != length(this.genes[c])) {
                    gene.TAD <- IRanges::subsetByOverlaps(TAD.gr, this.genes[c])
                    gene.TAD <- check_TAD(gene.TAD, TAD.gr, this.genes[c], this.genes)
                } else {
                    gene.TAD <- this.TAD
                }
                #TAD_COORDINATES = paste0(this.TAD),
                interactions <- rbind(interactions, data.frame(data.frame(this.region)[,
                    c(GR_header_short, "cluster", IDs)], TAD_COORDINATES = paste0(gene.TAD),
                    CORRELATED_GENE = paste(GenomicRanges::mcols(this.genes)[c, "gene_id"]),
                    CORRELATED_GENE_CHR = GenomicRanges::seqnames(this.genes[c]), CORRELATED_GENE_PROMOTER_START = promoter_start,
                    CORRELATED_GENE_PROMOTER_END = promoter_end, CORRELATION = cor[c]))
            }
        }
        if (length(interactions) > 0)
            return(GenomicRanges::makeGRangesFromDataFrame(interactions, keep.extra.columns = TRUE))
    }
}


##################################################################
# function: get nearest gene for a region
##################################################################

get_nearest_gene <- function(i, regions.gr, genes,
    IDs) {

    this.region <- regions.gr[i]

    nearest <- genes[GenomicRanges::nearest(this.region,
        genes)]
    distance.to.nearest <- GenomicRanges::distance(this.region,
        nearest)

    promoter_start <- GenomicRanges::start(GenomicRanges::promoters(nearest))
    promoter_end <- GenomicRanges::end(GenomicRanges::promoters(nearest))

    if (as.character(GenomicRanges::strand(nearest)) ==
        "-") {
        promoter_start <- GenomicRanges::end(GenomicRanges::promoters(nearest))
        promoter_end <- GenomicRanges::start(GenomicRanges::promoters(nearest))
    }

    interactions <- data.frame(data.frame(this.region)[,
        c(GR_header_short, "cluster", IDs)],
        NEAREST_GENE = GenomicRanges::mcols(nearest)[,
            "gene_id"], NEAREST_GENE_CHR = GenomicRanges::seqnames(nearest),
        NEAREST_GENE_PROMOTER_START = promoter_start,
        NEAREST_GENE_PROMOTER_END = promoter_end,
        DISTANCE_TO_NEAREST = distance.to.nearest)

    return(GenomicRanges::makeGRangesFromDataFrame(interactions,
        keep.extra.columns = TRUE))
}


##################################################################
# function: correlate probabilities with gene expression values
#(in same TAD; per cluster)
##################################################################

get_units <- function(regions.gr, expr.gr,
    TAD.gr, IDs, cores, threshold, txdb,
    nearest = FALSE) {

    if (nearest == FALSE) {
        # get correlation for each differential region:
        list <- parallel::mclapply(seq(length(regions.gr)),
            function(x) get_correlation(x,
                threshold, regions.gr,
                expr.gr, TAD.gr, IDs),
            mc.cores = cores)
        #print(list)
    } else {

        # get nearest gene:
        genes <- GenomicFeatures::genes(txdb)
        #GenomeInfoDb::seqlevels(genes) = paste0('chr', gsub('chr|Chr','', GenomeInfoDb::seqlevels(genes)))#what is this doing?
        GenomeInfoDb::seqlevels(genes,
            pruning.mode = "coarse") <- GenomeInfoDb::seqlevels(regions.gr)
        GenomeInfoDb::seqlengths(regions.gr) <- GenomeInfoDb::seqlengths(genes)
        list <- parallel::mclapply(seq(length(regions.gr)),
            function(x) get_nearest_gene(x,
                regions.gr, genes, IDs),
            mc.cores = cores)
    }
    #if(is.empty(list) == TRUE){
    #  message <- 'No target genes were found for the clusters'
    #  stop(message);
    #}
    units <- do.call("c", unname(unlist(list)))
    return(units)
}
