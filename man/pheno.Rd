\name{pheno}
\alias{pheno}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
	data
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
	FID, IID, and igrowth of the HapMap samples used in the intrinsic growth rate example. Column 3 header should be the same as \code{\link{phenoname}}.
}
\usage{data(pheno)}
\format{
  A data frame with 534 observations on the following 3 variables.
  \describe{
    \item{\code{FID}}{a character vector}
    \item{\code{IID}}{a character vector}
    \item{\code{igrowth}}{a numeric vector}
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
	Im, HK. et al. (2012) Mixed effects modeling of proliferation rates in cell-based models: consequence for pharmacogenomics and cancer. PLoS Genetics 8(2):e1002525.

	http://www.plosgenetics.org/article/info\%3Adoi\%2F10.1371\%2Fjournal.pgen.1002525

}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(pheno)
str(pheno)
hist(pheno[,3])
}
\keyword{datasets}
