% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{plotSummary}
\alias{plotSummary}
\title{plots boxplots of the median enhancer prediction values of the enhancers in the condition-specific clusters}
\usage{
plotSummary(D, num_plots = 9)
}
\arguments{
\item{D}{The getDynamics() output file of containing the GRanges object with the differential enhaners}

\item{num_plots}{Maximal number of cluster whose plots should be displayed (clusters are sorted by their sizes). 
This parameter can be set in case the number of clusters is very high (20+). Default is 20.}
}
\value{
a ggplot2 object containing a boxplot for each cluster
}
\description{
plots boxplots of the median enhancer prediction values of the enhancers in the condition-specific clusters
}
\examples{
#get the output of getDynamics()
files <- c(system.file("extdata", "Condition1.H3K4me1.bam", package="crupR"),
           system.file("extdata", "Condition1.H3K4me3.bam", package="crupR"),
           system.file("extdata", "Condition1.H3K27ac.bam", package="crupR"),
           system.file("extdata", "Condition2.H3K4me1.bam", package="crupR"),
           system.file("extdata", "Condition2.H3K4me3.bam", package="crupR"),
           system.file("extdata", "Condition2.H3K27ac.bam", package="crupR"))
inputs <- rep(system.file("extdata", "Condition1.Input.bam", package="crupR"), 3)
inputs2 <- rep(system.file("extdata", "Condition2.Input.bam", package="crupR"), 3)  
metaData <- data.frame(HM = rep(c("H3K4me1", "H3K4me3", "H3K27ac"),2),
                       condition = c(1,1,1,2,2,2), replicate = c(1,1,1,1,1,1),
                       bamFile = files, inputFile = c(inputs, inputs2))
clusters <- readRDS(system.file("extdata", "differential_enhancers.rds", package="crupR"))
dynamics <- list("metaData" = metaData, "sumFile" = clusters)
plotSummary(dynamics)

}
