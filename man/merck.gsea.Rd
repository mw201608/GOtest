\name{merck.gsea}
\alias{merck.gsea}
\title{Gene Set Enrichment Analysis Using Merck Collections}
\description{
  These functions perform gene set enrichment analyses based on Merck collections.
}
\usage{
merck.gsea(x, query.population=NULL, genesets=c('GO CC', 'GO MF', 'GO BP',
 'KEGG PATHWAYS', 'INGENUITY PATHWAYS', 'GENEGO PATHWAYS'),
 background='common', name.x='Input', name.go='Merck',
 method=c('hypergeometric', 'GSEA', 'logitreg'), adj="BH", species=c('human', 'mouse', 'rat'),
 ReportOverlap=TRUE, ncores=1, gsea.alpha=1,
 permutations=ifelse(method=='GSEA', 1000, 0), iseed=12345)
}
\arguments{
 \item{x}{a data.frame of query gene signatures. See \code{Details}.}
 \item{query.population}{a population of the genes from which the query gene signatures are collected.}
 \item{genesets}{Merck gene set collections to be tested.}
 \item{background}{an integer or a keyword. See explanation in function \code{GOtest}.}
 \item{name.x}{a character string specifing the name of the query.}
 \item{name.go}{a character string specifing the name of the function gene sets.}
 \item{method}{algorithm to performing the enrichment test. See \code{Details}.}
 \item{adj}{approach to correct for multiple tests. See function \code{p.adjust}}
 \item{species}{Specifying the species for the input gene signatures.}
 \item{ReportOverlap}{whether to output overlapping elements.}
 \item{ncores}{number of CPU cores to be used for parallel computing.}
 \item{gsea.alpha}{power to scale the weights: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 or larger (over-weighted)}
 \item{permutations}{number of permutations for computing P value significance and controlling FDR for \code{GSEA} and \code{logitreg} methods. See \code{\link{GOtest}} for details.}
 \item{iseed}{seed for random number generation in permutations.}
}
\details{
Method details can be found in function \code{\link{GOtest}}.
}
\value{
See output of function \code{GOtest}.
}
\seealso{\code{\link{GOtest}}, \code{\link{plotOverlap}}, \code{\link{p.adjust}}, \code{\link{GSEA}}}
\examples{
###Example usage 2: weighted enrichment tests

#We will again make use of the MacArthur gene set and the gene universe of HGNC approved symbols,
# so make sure they have been loaded as in Example 1.
\dontrun{
library(GOtest)
MAGenes=curated.genesets(c('MacArthur'))
universe=curated.genesets(c('HGNC_universe'))$Gene
(n=length(universe))

#In this example, we will try both the hypergeometric test and weighted enrichment tests, 
#including GSEA and logistic regression, by geneating a toy dataset
#through simulation of random gene-phenotype associations.

set.seed(123)
toy = data.frame(Gene=universe, Phenotype='Simulated', Z=rnorm(n, 0, 1), stringsAsFactors=FALSE)

#Next, we are going to run GSEA weighted enrichment tests on the full test dataset by using 4 CPU cores on a multicore system with at least 4 cores.

fit = merck.gsea(x = toy, query.population = universe, backgroup = 'query', name.x = 'Toy', method = 'GSEA', permutations = 1000, ncores = 4)
head(fit)

}
}
