msigdb.gsea=function(x, query.population=NULL, genesets=c('c2.cp', 'c5.go.bp', 'c5.go.cc', 'c5.go.mf'),
 background='common', name.x='Input', name.go='MSigDB', method=c('hypergeometric', 'GSEA', 'logitreg'),
 adj="BH", species=c('human','mouse'), ReportOverlap=TRUE, ncores=1, gsea.alpha=1, permutations=ifelse(method=='GSEA', 1000, 0), iseed=12345){
	go=msigdb.genesets(sets=genesets, type='symbols', species=species, return.data.frame=TRUE)
	GOtest(x=x,go=go,query.population=query.population,background=background,name.x=name.x,name.go=name.go,method=method,adj=adj,ReportOverlap=ReportOverlap,ncores=ncores,gsea.alpha=gsea.alpha,permutations=permutations, iseed=iseed)
}
