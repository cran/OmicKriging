\name{Xcovamat}
\alias{Xcovamat}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
	data
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
	matrix of covariates (10 principal components from the genetic relationship matrix and 10 principal components from the gene expression correlation matrix) from the intrinsic growth rate example with the IID list as the rownames
}
\usage{data(Xcovamat)}
\format{
  The format is:
 num [1:159, 1:20] 0.0828 0.079 0.0771 0.0807 0.0787 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:159] "NA10846" "NA10847" "NA12144" "NA12145" ...
  ..$ : chr [1:20] "PC1.SNPs" "PC2.SNPs" "PC3.SNPs" "PC4.SNPs" ...
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
        Im, HK. et al. (2012) Mixed effects modeling of proliferation rates in cell-based models: consequence for pharmacogenomics and cancer. PLoS Genetics 8(2):e1002525.

        http://www.plosgenetics.org/article/info\%3Adoi\%2F10.1371\%2Fjournal.pgen.1002525

}
\examples{
data(Xcovamat)
str(Xcovamat)
dim(Xcovamat)
Xcovamat[1:10,]
}
\keyword{datasets}
