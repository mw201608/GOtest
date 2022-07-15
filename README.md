# GOtest: Gene Ontology and Set Enrichment Test

## Dependency
R package `msigdbi` is required to run this tool. You can install it from github by
```
devtools::install_github("mw201608/msigdb")
```

## Installation
```
devtools::install_github("mw201608/GOtest")
```

## Example usage 1: hypergeometric test with MSigDB GO/pathway annotations

In this example, we will load pre-installed 23 functional gene sets curated by the MacArthur's Lab (https://github.com/macarthur-lab/gene_lists) and then apply the hypergeometric test to evaluate the overlap in the MSigDB GO/pathway annotations. 
```
library(GOtest)
MAGenes = curated.genesets(c('MacArthur'))
head(MAGenes)
```
Some enrichment analysis methods, eg the hypergeometric test, require a gene universe or a population of gene background. Here we will use a pre-installed set of approved symbols for protein-coding genes by HGNC.
```
universe=curated.genesets(c('HGNC_universe'))$Gene
str(universe)
```
For details about function `curated.genesets`, check help `?curated.genesets`.

Now let us run enrichment of MacArthur gene sets against the canonical pathways in MSigDB annotation.

```
result = msigdb.gsea(x=MAGenes, query.population=universe, genesets=c('c2.cp'), background='query', name.x='MacArthur')
head(result)
```

## Example usage 2: weighted enrichment tests

We will again make use of the MacArthur gene set and the gene universe of HGNC approved symbols, so make sure they have been loaded as in Example 1.
```
library(GOtest)
MAGenes = curated.genesets(c('MacArthur'))
universe = curated.genesets(c('HGNC_universe'))$Gene
(n = length(universe))
```
### Create a toy dataset
In this example, we will try both the hypergeometric test and weighted enrichment tests, including GSEA and logistic regression, by geneating a toy dataset through simulation of random gene-phenotype associations. 
```
set.seed(123)
toy = data.frame(Gene = universe, Phenotype = 'Simulated', Z = rnorm(n, 0, 1), stringsAsFactors = FALSE)
```

Select genes with absolute Z value larger than 3 and separate them into up and down groups based on the sign of Z value.
```
toy.3 = toy[abs(toy$Z) > 3, ]
toy.3$Direction = ifelse(toy.3$Z > 0, 'Up', 'Down')
```

### Hypergeometric test
Run the hypergeometric test on both groups against the MacAuther gene sets:
```
fit1 = GOtest(x = toy.3[, c('Gene', 'Direction')], go = MAGenes, query.population = universe, background = 'query', name.x = 'Toy', name.go = 'MacArthur', method = 'hypergeometric')
```
As expected, no significant enrichment identified:
```
head(fit1)
```

To test the enrichment in the MSigDB gene sets using the hypergeometric test, run
```
fit2 = msigdb.gsea(x = toy.3[, c('Gene', 'Direction')], genesets = c('c2.cp'), query.population = universe, background = 'query')
head(fit2)
```

### GSEA
Next, we are going to run weighted enrichment tests on the full test dataset by using GSEA or logistic regression. First, run GSEA:
```
fit3 = GOtest(x = toy, go = MAGenes, name.x = 'Toy', name.go = 'MacArthur', query.population = universe, background = 'query', method = 'GSEA')
head(fit3)
```
Again there is no significant enrichment. Let us check the GSEA running enrichment score plot for the top 10 MacArthur terms:
```
plotGseaEnrTable(GseaTable = fit2[1:10, ], x = toy, go = MAGenes)
```
![GSEA running enrichment score example](gsea.res.png)

To test the enrichment in the MSigDB gene sets using GSEA, run
```
fit4 = msigdb.gsea(x = toy, genesets = c('c2.cp'), query.population = universe, background = 'query', name.x = 'Toy', method = 'GSEA', permutations = 1000)
head(fit4)
```

### Logistic regression
Run logistic regression:
```
fit5 = GOtest(x = toy, go = MAGenes, name.x = 'Toy', name.go = 'MacArthur', query.population = universe, background = 'query', method = 'logitreg')
head(fit5)
```
