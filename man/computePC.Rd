\name{computePC}
\alias{computePC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
	Perform principal components analysis on a correlation matrix	
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
	This function calls the software GCTA to calculate the first 10 principal components from a correlation matrix (filename.grm.gz and filename.grm.id files are needed).
}
\usage{
computePC(fullheader, idfile=idfile, gctaname = "gcta")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fullheader}{
%%     ~~Describe \code{fullheader} here~~
	path to .grm.gz and .grm.id files 
}
  \item{idfile}{
%%     ~~Describe \code{idfile} here~~
	 path to a file listing the subset of individuals to include in the principal components analysis (2 columns labeled: FID, IID)	
}
  \item{gctaname}{
%%     ~~Describe \code{gctaname} here~~
	path to gcta executable
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
	This function creates output files filename.eigenval and filename.eigenvec
}
\references{
%% ~put references to the literature/web site here ~
        Wheeler HE, et al. (2013) Poly-Omic Prediction of Complex Traits: OmicKriging. arXiv:1303.1788 
        http://arxiv.org/abs/1303.1788

        Yang J, Lee S, Goddard M, Visscher P (2011) GCTA: a tool for genome-wide complex trait analysis. American Journal of Human Genetics 88: 76-82.
	http://www.complextraitgenomics.com/software/gcta/
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
%function (fullheader, idfile=idfile, gctaname = "gcta") 
%{
%    runPC = gctaname \%&\% " --grm " \%&\% fullheader \%&\% " --keep " \%&\% 
%        idfile \%&\% " --pca 10 --out " \%&\% fullheader
%    print(runPC)
%    system(runPC)
%  }
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
