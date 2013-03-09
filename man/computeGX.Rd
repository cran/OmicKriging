\name{computeGX}
\alias{computeGX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
	Compute the correlation matrix for an "omics" dataset
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
	This function computes the correlation matrix of all individuals in an "omics" dataset (e.g. transcriptome, metabolome, methylome, etc.).

}
\usage{
computeGX(genefile, outheader, idfile = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genefile}{
%%     ~~Describe \code{genefile} here~~
	path to "omics" dataset (i.e. column labels in a transcriptome dataset: FID, IID, gene1, gene2, ...)
}
  \item{outheader}{
%%     ~~Describe \code{outheader} here~~
	output path and file name
}
  \item{idfile}{
%%     ~~Describe \code{idfile} here~~
	path to a file listing the subset of individuals to include in the correlation matrix calculation (2 columns labeled: FID, IID)
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
	a matrix of individual correlations and output files filename.grm.gz and filename.grm.id
}
\references{
%% ~put references to the literature/web site here ~
        Wheeler HE, et al. (2013) Poly-Omic Prediction of Complex Traits: OmicKriging. arXiv:1303.1788 
        http://arxiv.org/abs/1303.1788
}
\author{
%%  ~~who you are~~
	Hae Kyung Im 

	Heather E. Wheeler
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%}
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.

%## The function is currently defined as
%function (genefile, outheader, idfile = NULL) 
%{
%    genedata = read.delim(genefile, sep = "", as.is = T, header = T)
%    if (!is.null(idfile)) {
%        iddata = read.table(idfile, header = T, as.is = T)
%        genedata = merge(iddata, genedata, by.x = c("FID", "IID"), 
%            by.y = c("FID", "IID"))
%    }
%    cordata.id = genedata[c("FID", "IID")]
%    genemat = as.matrix(genedata[, !(names(genedata) \%in\% c("FID", 
%        "IID"))])
%    if (sum(is.na(genemat)) > 0) 
%        print("Impute first")
%    cormat = cor(t(genemat))
%    ng = ncol(genemat)
%    cordata = cbind(c(col(cormat)), c(row(cormat)), ng, c(cormat))
%    cordata = subset(cordata, cordata[, 1] >= cordata[, 2])
%    write.table(cordata, file = outheader \%&\% ".grm", col = F, 
%        row = F, quote = F)
%    system("gzip " \%&\% outheader \%&\% ".grm ")
%    write.table(cordata.id, file = outheader \%&\% ".grm.id", col = F, 
%        row = F, quote = F)
%    return(cormat)
%  }
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
