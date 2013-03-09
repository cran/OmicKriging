\name{computeGRM}
\alias{computeGRM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Compute the genetic relatedness matrix
}
\description{
  	This function calls the software GCTA to calculate the genetic relationship matrix from all autosomal SNP genotype data in PLINK binary PED file format.
}
\usage{
computeGRM(bedfullheader, grmfullheader = "genomic", 
gctaname = "gcta", idfile = NULL, snpfile = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bedfullheader}{
%%     ~~Describe \code{bedfullheader} here~~
	path to PLINK binary .bed/.bim/.fam files
}
  \item{grmfullheader}{
%%     ~~Describe \code{grmfullheader} here~~
	output path and file name
}
  \item{gctaname}{
%%     ~~Describe \code{gctaname} here~~
	path to gcta executable
}
  \item{idfile}{
%%     ~~Describe \code{idfile} here~~
	path to a file listing the subset of individuals to include in GRM calculation (2 columns labeled: FID, IID)
}
  \item{snpfile}{
%%     ~~Describe \code{snpfile} here~~
	path to a file listing a subset of SNPs, one per line and no column header, to included in the GRM calculation, default is to include all autosomal SNPs in the .bim file 
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
	This function creates the GCTA GRM output files filename.grm.gz and filename.grm.id	
}
\references{
%% ~put references to the literature/web site here ~
	Wheeler HE, et al. (2013) Poly-Omic Prediction of Complex Traits: OmicKriging. arXiv:1303.1788 
	http://arxiv.org/abs/1303.1788

	Yang J, Lee S, Goddard M, Visscher P (2011) GCTA: a tool for genome-wide complex trait analysis. American Journal of Human Genetics 88: 76-82.
	http://www.complextraitgenomics.com/software/gcta/

	Purcell, S. et al. (2007) PLINK: a tool set for whole-genome association and population-based linkage analyses. American Journal of Human Genetics 81:559-75.
	http://pngu.mgh.harvard.edu/~purcell/plink/
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
%function (bedfullheader, grmfullheader = "genomic", gctaname = "gcta", 
%    idfile = NULL) 
%{
%    if (is.null(idfile)) 
%        idfile = bedfullheader \%&\% ".fam"
%    runGRM = gctaname \%&\% " --bfile " \%&\% bedfullheader \%&\% " --autosome --make-grm --keep " \%&\% 
%        idfile \%&\% " --out " \%&\% grmfullheader
%    system(runGRM)
%  }
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
