\name{corlist}
\alias{corlist}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
	data
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
List of 2 correlation matrices, \code{corlist[1]} is the genetic relationship matrix of the intrinsic growth rate samples calculated by GCTA from 43,182 regulome SNPs. \code{corlist[2]} is the gene expression correlation matrix from the intrinsic growth rate samples. The function \code{\link{readcorlist}} generated this list.
}
\usage{data(corlist)}
\format{
  The format is:
List of 2
 $ : num [1:159, 1:159] 1.062 0.0909 0.0697 0.0884 0.085 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:159] "NA12750" "NA11831" "NA12146" "NA11882" ...
  .. ..$ : chr [1:159] "NA12750" "NA11831" "NA12146" "NA11882" ...
 $ : num [1:159, 1:159] 1 0.969 0.959 0.964 0.961 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:159] "NA10846" "NA10847" "NA12144" "NA12145" ...
  .. ..$ : chr [1:159] "NA10846" "NA10847" "NA12144" "NA12145" ...
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
data(corlist)
str(corlist)
}
\keyword{datasets}
