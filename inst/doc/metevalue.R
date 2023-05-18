## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  metevalue.[DMR](
#    methyrate,                # Output file name of [DMR]
#    [DMR].output,             # Output file name of [DMR] with e-value of each region
#    adjust.methods = "BH",    # Adjust methods of e-value
#    sep = "\t",               # seperator, default is the TAB key.
#    bheader = FALSE           # A logical value indicating whether the [DMR].output file
#                              # contains the names of the variables as its first line
#  )

## ----eval=FALSE---------------------------------------------------------------
#  # Here  `[DMR]` coudle be one of `methylKit`, `biseq`, `DMRfinder` or `metilene`.
#  method_in_use = "[DMR]"
#  result = evalue_buildin_var_fmt_nm(
#            methyrate,              # Data frame of the methylation rate
#            DMR_evalue_output,      # Data frame of output data corresponding to the
#                                    # "method" option
#            method = method_in_use) # DMR: "metilene", "biseq", "DMRfinder" or "methylKit"
#  result = list(a = result$a,
#                b = result$b,
#                a_b = evalue_buildin_sql(result$a, result$b, method = method_in_use))
#  result = varevalue.metilene(result$a, result$b, result$a_b)

## ----eval=FALSE---------------------------------------------------------------
#  library(metevalue)
#  
#  ####Simulation Data ####
#  set.seed(1234)
#  
#  # methylKit is a Bioconductor package
#  library(methylKit)
#  file.list=list( system.file("extdata",
#                              "test1.myCpG.txt", package = "methylKit"),
#                  system.file("extdata",
#                              "test2.myCpG.txt", package = "methylKit"),
#                  system.file("extdata",
#                              "control1.myCpG.txt", package = "methylKit"),
#                  system.file("extdata",
#                              "control2.myCpG.txt", package = "methylKit") )
#  
#  
#  # read the files to a methylRawList object: myobj
#  myobj=methRead(file.list,
#                 sample.id=list("test1","test2","ctrl1","ctrl2"),
#                 assembly="hg18",
#                 treatment=c(1,1,0,0),
#                 context="CpG"
#  )
#  
#  meth=unite(myobj, destrand=FALSE)
#  meth.C <- getData(meth)[,seq(6,ncol(meth),3)]
#  meth.T <- getData(meth)[,seq(7,ncol(meth),3)]
#  mr <- meth.C/(meth.C + meth.T)
#  chr_pos = getData(meth)[,1:2]
#  methyrate = data.frame(chr_pos,mr)
#  names(methyrate) = c('chr', 'pos', rep('g1',2), rep('g2',2))
#  region<-tileMethylCounts(myobj)
#  meth<-unite(region,destrand=F)
#  myDiff<-calculateDiffMeth(meth)
#  met_all<-getMethylDiff(myDiff,type="all")
#  
#  example_tempfiles = tempfile(c("rate_combine", "methylKit_DMR_raw"))
#  tempdir()
#  write.table(methyrate, file=example_tempfiles[1], row.names=F, col.names=T, quote=F, sep='\t')
#  write.table (met_all, file=example_tempfiles[2], sep ="\t", row.names =F, col.names =T, quote =F)

## ----eval=FALSE---------------------------------------------------------------
#  result = metevalue.methylKit(example_tempfiles[1], example_tempfiles[2], bheader = T)
#  str(result)

## ----eval=FALSE---------------------------------------------------------------
#  result = evalue_buildin_var_fmt_nm(methyrate, met_all, method="methylKit")
#  result = list(a = result$a,
#                b = result$b,
#                a_b = evalue_buildin_sql(result$a, result$b, method="methylKit"))
#  result = varevalue.metilene(result$a, result$b, result$a_b)
#  str(result)

## ----eval=FALSE---------------------------------------------------------------
#  library(BiSeq)
#  library(dplyr)
#  data(rrbs)
#  rrbs.rel <- rawToRel(rrbs)
#  methyrate <- methLevel(rrbs.rel)
#  methyrate <- data.frame(methyrate)
#  methyrateq = cbind(rows = as.numeric(row.names(methyrate)), methyrate)
#  methypos = data.frame(rows = as.numeric(row.names(methyrate)), rowRanges(rrbs))
#  methyrate = left_join(methypos, methyrateq)
#  methyrate = methyrate[,c(2,3,7:16)]
#  names(methyrate) <- c('chr','pos',rep('g1',5),rep('g2',5))
#  
#  rrbs.clust.unlim <- clusterSites(object = rrbs,perc.samples = 3/4,min.sites = 20,max.dist = 100)
#  
#  clusterSitesToGR(rrbs.clust.unlim)
#  ind.cov <- totalReads(rrbs.clust.unlim) > 0
#  
#  quant <- quantile(totalReads(rrbs.clust.unlim)[ind.cov])
#  rrbs.clust.lim <- limitCov(rrbs.clust.unlim, maxCov = quant)
#  predictedMeth <- predictMeth(object = rrbs.clust.lim)
#  
#  test<- predictedMeth[, colData(predictedMeth)$group == "test"]
#  control <- predictedMeth[, colData(predictedMeth)$group == "control"]
#  mean.test <- rowMeans(methLevel(test))
#  mean.control <- rowMeans(methLevel(control))
#  
#  betaResults <- betaRegression(formula = ~group,link = "probit",object = predictedMeth,type = "BR")
#  vario <- makeVariogram(betaResults)
#  vario.sm <- smoothVariogram(vario, sill = 0.9)
#  
#  locCor <- estLocCor(vario.sm)
#  clusters.rej <- testClusters(locCor)
#  clusters.trimmed <- trimClusters(clusters.rej)
#  DMRs <- findDMRs(clusters.trimmed,max.dist = 100,diff.dir = TRUE)
#  
#  
#  example_tempfiles = tempfile(c('rate_combine', 'BiSeq_DMR'))
#  write.table(methyrate, example_tempfiles[1], row.names=F, col.names=T, quote=F, sep='\t')
#  write.table(DMRs, example_tempfiles[2], quote=F, row.names = F,col.names = F, sep = '\t')

## ----eval=FALSE---------------------------------------------------------------
#  hh <- data.frame(DMRs)
#  result = evalue_buildin_var_fmt_nm(rate_combine, hh, method="biseq")
#  result = list(a = result$a,  b = result$b,
#                a_b = evalue_buildin_sql(result$a, result$b,method="biseq"))
#  result = varevalue.metilene(result$a, result$b, result$a_b)
#  str(result)

## ----eval=FALSE---------------------------------------------------------------
#  rate_combine <- read.table("rate_combine_DMRfinder", header = T)
#  head(rate_combine)
#  
#  DMRs <- read.table("DMRfinder_DMR", header = T)
#  head(DMRs)
#  

## ----eval=FALSE---------------------------------------------------------------
#  result <- evalue.DMRfinder('rate_combine_DMRfinder', 'DMRfinder_DMR')
#  head(result)

## ----eval=FALSE---------------------------------------------------------------
#  result = evalue_buildin_var_fmt_nm(rate_combine, DMRs, method="DMRfinder")
#  result = list(a = result$a,
#                b = result$b,
#                a_b = evalue_buildin_sql(result$a, result$b, method="DMRfinder"))
#  result = varevalue.metilene(result$a, result$b, result$a_b)
#  head(result)

## ----eval=FALSE---------------------------------------------------------------
#  input <- read.table("metilene.input", header = T)
#  head(input)
#  
#  out <- read.table("metilene.out", header = F)
#  head(out)
#  

## ----eval=FALSE---------------------------------------------------------------
#  result <- evalue.metilene('metilene.input', 'metilene.out')
#  head(result)

## ----eval=FALSE---------------------------------------------------------------
#  result = evalue_buildin_var_fmt_nm(input, out, method="metilene")
#  result = list(a = result$a,
#                b = result$b,
#                a_b = evalue_buildin_sql(result$a, result$b, method="metilene"))
#  result = varevalue.metilene(result$a, result$b, result$a_b)
#  head(result)

## ----eval=FALSE---------------------------------------------------------------
#  input <- read.table("metilene.input", header = T)
#  e_value <- varevalue.signle_general(methyrate=input, group1_name='g1', group2_name='g2', chr='chr21', start=9439679, end=9439679)
#  head(e_value)

## ----eval=FALSE---------------------------------------------------------------
#  input <- read.table("desq_out", header = T)
#  data_e <- metevalue.RNA_general(input, group1_name='treated', group2_name='untreated')
#  head(data_e)

