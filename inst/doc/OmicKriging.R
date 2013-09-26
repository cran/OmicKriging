### R code from vignette source 'OmicKriging.Rnw'

###################################################
### code chunk number 1: OmicKriging.Rnw:46-54
###################################################
library(OmicKriging)

"%&%" <- function(a, b) paste(a, b, sep="")
gdsFile <-"gdsTemp.gds"
ok.dir <-  path.package('OmicKriging') %&% "/doc/vignette_data/"
bFile <- ok.dir %&% "ig_genotypes"
expFile <- ok.dir %&% "ig_gene_subset.txt.gz"
phenoFile <- ok.dir %&% "ig_pheno.txt"


###################################################
### code chunk number 2: OmicKriging.Rnw:61-62
###################################################
pheno <- read.table(phenoFile, header = T)


###################################################
### code chunk number 3: OmicKriging.Rnw:69-70
###################################################
grmMat <- read_GRMBin(bFile)


###################################################
### code chunk number 4: OmicKriging.Rnw:77-78
###################################################
convert_genotype_data(bFile = bFile, gdsFile = gdsFile)


###################################################
### code chunk number 5: OmicKriging.Rnw:85-86
###################################################
grmMat <- make_GRM(gdsFile = gdsFile)


###################################################
### code chunk number 6: OmicKriging.Rnw:96-97
###################################################
gxmMat <- make_GXM(expFile = expFile)


###################################################
### code chunk number 7: OmicKriging.Rnw:108-113
###################################################
pcMatXM <- make_PCs_irlba(gxmMat, n.top = 10)

pcMatGM <- make_PCs_irlba(grmMat, n.top = 10)

pcMat <- cbind(pcMatGM, pcMatXM[match(rownames(pcMatGM), rownames(pcMatXM)),])


###################################################
### code chunk number 8: OmicKriging.Rnw:127-133
###################################################
result <- krigr_cross_validation(pheno.df = pheno,
	cor.list = list(grmMat, gxmMat),
	h2.vec = c(0.5, 0.5),
	covar.mat = pcMat,
	ncore = 2,
	nfold = "LOOCV")


