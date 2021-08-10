\name{GOtest}
\docType{package}
\alias{GOtest}
\title{Gene Ontology Functional Enrichment Test}
\description{
  This function performs functional enrichment analyses.
}
\usage{
GOtest(x, go, query.population=NULL, background='common', 
	name.x='Input', name.go='Category',
	method=c('hypergeometric','GSEA','logitreg'), adj="BH",
	ReportOverlap=TRUE, ncores=1, gsea.alpha=1, 
	permutations=ifelse(method=='GSEA',1000,0), iseed=12345)
}
\arguments{
 \item{x}{a data.frame of query gene sets, with the first column the gene symbols and the second column the set names. A third column is needed for the \code{GSEA} and the \code{logitreg} method.}
 \item{go}{a data.frame of GO gene sets, with the first column the gene symbols and the second column being GO categories. The GO category can be written in a format of System:Term which will be decoded in the output.}
 \item{query.population}{a population of the genes from which the query gene sets are collected.}
 \item{background}{an integer or a keyword specifying the background for \code{hypergeometric} test method. See \code{Details}.}
 \item{name.x}{a character string specifing the name of the query.}
 \item{name.go}{a character string specifing the name of the function gene sets.}
 \item{method}{algorithm to performing the enrichment test. See \code{Details}.}
 \item{adj}{approach to correct for multiple tests for methods \code{hypergeometric} and \code{logitreg}. See function \code{\link[stats]{p.adjust}}.}
 \item{ReportOverlap}{whether to output overlapping elements.}
 \item{ncores}{number of CPU cores to be used for parallel computing. Using multiple cores might not lead to a significant speed gain for \code{hypergeometric} or \code{logitreg}.}
 \item{gsea.alpha}{power to scale the weights in \code{GSEA}: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 or larger (over-weighted)}
 \item{permutations}{number of permutations for computing significance and controlling FDR for \code{GSEA} and \code{logitreg} method.}
 \item{iseed}{seed for random number generation in permutations.}
}
\details{
If \code{method} is "\code{hypergeometric}", \code{x} must have at least two columns, with the first column the gene IDs and the second column the group name(s) for each gene in the row. Multiple tests will be controled by Benjamini-Hochberg's FDR.

If \code{method} is "\code{GSEA}" (see reference 1), \code{x} must have three columns, the first column contains gene IDs, the second column specifies phenotype name or group name, and the third column specifies gene-phenotype correlation strengths (eg, log fold change of differential expression, t-test statistics, or minus log P values), whose sign and absolute values can be used to rank the genes within a phenotype/group.

If \code{method} is "\code{logitreg}" (logistic regression; see reference 2), \code{x} should have at least three columns. The first three columns are the same as those for "\code{GSEA}". Any additional columns may contain covariates.

For method "\code{GSEA}", P value and multiple tests correction will be based on permutation analyses. Permutations can also be used for the "\code{logitreg}" method (not recommended as logistic regression is very slow).

The argument \code{background} specifies the size of the gene set background for \code{hypergeometic} test. If \code{background} is an integer, it must be larger than the size of the union of query gene set and annotated gene set. Alternatively, the size of the gene set background will be determined by one of the following keywords:
\itemize{
\item \code{query}: the genes supplied by \code{query.population};
\item \code{annotation}: the genes present in database \code{go};
\item \code{intersection}: the intersection between \code{query.population} and \code{go};
\item \code{common}: the same as \code{intersection};
\item \code{union}: the union of \code{query.population} and \code{go};
\item \code{query.population} should not be NULL if \code{background} is \code{query}, \code{common} or \code{union}.
}
}
\value{
A data.frame with components:
 \item{System}{GO gene set system classifications}
 \item{Category}{GO gene set names}
 \item{Input}{Query gene set names}
 \item{Overlap.Size}{Size of the overlap between GO gene set and query gene set}
 \item{Input.Size}{Size of the query gene set}
 \item{Category.Size}{GO gene set size}
 \item{Background.Size}{Size of the background population for \code{hypergeometric} test}
 \item{FE}{Fold enrichment of the overlap for \code{hypergeometric} test}
 \item{ES}{\code{GSEA} enrichment score}
 \item{log.odds}{\code{logitreg} log odds}
 \item{SE}{\code{logitreg} standard error of log odds}
 \item{Z}{\code{logitreg} Z-statistic}
 \item{Pvalue}{P value significance of the overlap}
 \item{P.adj}{Adjusted P value}
 \item{Overlap.Items}{Overlapping items}
The column names starting with \code{Input} and \code{Category} will be specified by options \code{name.x} and \code{name.go}.
}
\references{
1. Subramanian et al (2005). Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. \emph{Proceedings of the National Academy of Sciences of the United States of America} 102(43): 15545-15550.
2. Sartor et al (2009). LRpath: a logistic regression approach for identifying enriched biological groups in gene expression data. \emph{Bioinformatics} 25(2): 211-217.
}
\seealso{\code{\link{msigdb.gsea}}, \code{\link{plotOverlap}}, \code{\link{p.adjust}}, \code{\link{curated.genesets}}, \code{\link{GSEA}}}
\examples{
###Example usage 1: the hypergeometric test with MSigDB GO/pathway annotations

#In this example, we will load pre-installed 23 functional gene sets 
#curated by the MacArthur's Lab (https://github.com/macarthur-lab/gene_lists) 
# and then apply the hypergeometric test to evaluate the overlap in MSigDB GO/pathways.

library(GOtest)
MAGenes=curated.genesets(c('MacArthur'))
head(MAGenes)

#Some enrichment analysis methods, eg the hypergeometric test, 
#require a gene universe or a population of gene background. 
#Here we will use a pre-installed set of approved symbols for protein-coding genes by HGNC.

universe=curated.genesets(c('HGNC_universe'))$Gene
str(universe)

#For details about function curated.genesets, check help ?curated.genesets.

#Now let us run enrichment of MacArthur gene sets against the MSigDB canonical pathways.

result=msigdb.gsea(x=MAGenes, query.population=universe, genesets=c('c2.cp'),
	background='query', name.x='MacArthur')
head(result)

###Example usage 2: weighted enrichment tests

#We will again make use of the MacArthur gene set and the gene universe
#of HGNC approved symbols, so make sure they have been loaded as in Example 1.

library(GOtest)
MAGenes=curated.genesets(c('MacArthur'))
universe=curated.genesets(c('HGNC_universe'))$Gene
(n=length(universe))

#In this example, we will try both the hypergeometric test and weighted enrichment tests, 
#including GSEA and logistic regression, by geneating a toy dataset through simulation 
#of random gene-phenotype associations.

set.seed(123)
toy=data.frame(Gene=universe, Phenotype='Simulated', Z=rnorm(n,0,1), stringsAsFactors=FALSE)

#Select genes with absolute Z value larger than 3 
#and separate them into up and down groups based on the sign of Z value, 
#then run the hypergeometric test on both groups against the MacAuther gene sets:

toy.3=toy[abs(toy$Z)>3,]
toy.3$Direction=ifelse(toy.3$Z > 0, 'Up','Down')
fit1=GOtest(x=toy.3[,c('Gene','Direction')], go=MAGenes, query.population=universe,
	background='query', name.x='Toy', name.go='MacArthur', method='hypergeometric')

#As expected, no significant enrichment identified:

head(fit1)

\dontrun{
#Next, we are going to run weighted enrichment tests on the full test dataset 
#by using GSEA or logistic regression. First, run GSEA:

fit2=GOtest(x=toy, go=MAGenes, name.x='Toy', name.go='MacArthur', method='GSEA')
head(fit2)
#Again there is no significant enrichment. Let us check the GSEA running 
#enrichment score plot for the top 10 MacArthur terms:
plotGseaEnrTable(GseaTable=fit2[1:10,], x=toy, go=MAGenes)

#Run logistic regression:

fit3=GOtest(x=toy, go=MAGenes, name.x='Toy', name.go='MacArthur', method='logitreg')
head(fit3)
}
}
