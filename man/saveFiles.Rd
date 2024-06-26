% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/saveFiles.R
\name{saveFiles}
\alias{saveFiles}
\title{Function to save the different output objects of each step}
\usage{
saveFiles(data, modes, outdir, nearest = FALSE)
}
\arguments{
\item{data}{An output list of one of the steps 2 to 4 or the get_superEnhancer step}

\item{modes}{Formats in which the GRanges object should be saved. Following modes are available:
for the output of getEnhancers(): "bigWig" for a bigWig file, "rds" for an .rds file
for the output of getDynamics(): "beds" for saving each cluster in a seperate bed file
for the output of getTargets(): "UCSC" for a UCSC interaction file
for the output of getSE(): "bedGraph" for saving the peak calls in a bedGraph file, "bed" for saving the clusters of peaks in a bed file}

\item{outdir}{Output directory in which the files should be saved}

\item{nearest}{Only relevant, if you want to save the output of enhancerTargets. Specifies if the output was produced by the nearest gene mode of the function or not. Default is false.}
}
\value{
Nothing
}
\description{
Function to save the different output objects of each step
}
\examples{
\donttest{
saveFiles(data = prediction_1_1, modes = c("bigWig", "rds"), outdir = "/example/dir/")
saveFiles(data = dynamics, modes = "beds", outdir ="/example/dir/")
saveFiles(data = targets, modes = "UCSC", outdir = "/example/dir/", nearest = TRUE)}
}
