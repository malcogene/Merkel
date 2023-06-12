
dir.create("data")

load(file.path(getwd(),
               "../data/MCC4cohorts.gs.temp.rda"))
load(file.path(getwd(),
               "../data/merkel.compendium.rda"))
load(file.path(getwd(),
               "../data/p.KEGG.PID.BioCarta.RData"))

gs.primaryMeta <- gs.temp

tmp1 <- gs.temp[[1]]$pheno$`metastasis, primary, or cell line:ch1`
inx1 = rownames(gs.temp[[1]]$pheno)[which(tmp1 == "primary tumor")]
inx2 = rownames(gs.temp[[1]]$pheno)[which(tmp1 == "metastatic tumor")]
gs.primaryMeta[[1]]$x <- gs.temp[[1]]$x[, c(inx1, inx2)]
gs.primaryMeta[[1]]$y <-
  data.frame(y = c(rep(0, length(inx1)), rep(1, length(inx2))))

tmp1 <- gs.temp[[2]]$pheno$`primary_tumor:ch1`
inx1 = rownames(gs.temp[[2]]$pheno)[which(tmp1 == "yes")]
inx2 = rownames(gs.temp[[2]]$pheno)[which(tmp1 == "no")]
gs.primaryMeta[[2]]$x <- gs.temp[[2]]$x[, c(inx1, inx2)]
gs.primaryMeta[[2]]$y <-
  data.frame(y = c(rep(0, length(inx1)), rep(1, length(inx2))))

tmp1 <- gs.temp[[3]]$pheno$characteristics_ch1.1
inx1 = rownames(gs.temp[[3]]$pheno)[which(tmp1 == unique(tmp1)[4])]
inx2 = rownames(gs.temp[[3]]$pheno)[which(tmp1 %in% unique(tmp1)[c(2, 3, 6, 7)])]
gs.primaryMeta[[3]]$x <- gs.temp[[3]]$x[, c(inx1, inx2)]
gs.primaryMeta[[3]]$y <-
  data.frame(y = c(rep(0, length(inx1)), rep(1, length(inx2))))

tmp1 <- gs.temp[[4]]$pheno$source_name_ch1
inx1 = rownames(gs.temp[[4]]$pheno)[which(tmp1 == "MCC_Good prognosis_tumor cells")]
inx2 = rownames(gs.temp[[4]]$pheno)[which(tmp1 == "MCC_Poor prognosis_tumor cells")]
gs.primaryMeta[[4]]$x <- gs.temp[[4]]$x[, c(inx1, inx2)]
gs.primaryMeta[[4]]$y <-
  data.frame(y = c(rep(0, length(inx1)), rep(1, length(inx2))))

d.merged <- summ.gene(gs.primaryMeta)
s <- d.merged

seed <- getOption("mySeedV")
metares <-
  .metaGene(
    s,
    .05,
    nperm = 5,
    bedArgs = list(label.sig = c(0.01, 2)),
    seed = seed
  )
table.outRes

set.seed(seed)
key.l <- rep(list(c(0, 1)), 4)
label.l <- rep(list(c("C", "E")), 4)

s <- .tarExtract(d.merged, key.l, label.l)
names(s)


inx = which(p.KEGG.PID.BioCarta$source %in% c("KEGG", "PID"))
pw.in <-
  p.KEGG.PID.BioCarta$entrez_gene_ids[inx]
names(pw.in) <- do.call("c", p.KEGG.PID.BioCarta$pathway[inx])

alphas.seq = seq(0.01, 1, length = 10)
devCohorts = rep(1, 4)

res <-
  biPDSmlv2(
    s,
    CVal = "LOOCV",
    GlobalOp = "epsgo",
    devCohorts = devCohorts,
    verbosePlotTable = T,
    alphas.seq = alphas.seq,
  )

ress <-
  biPDSmlv2(
    s,
    CVal = "LOSOCV",
    GlobalOp = "none",
    devCohorts = devCohorts,
    verbosePlotTable = T,
    alphas.seq = 1,
    authors = authors
  )


res$coefDf
res$cv.plot.list[[3]]
res$table.cv
res$table.val
pdf("./results/res.pdf", width = 5, height = 4)
res$cv.plot.list
res$val.plot.list
dev.off()


require(immunedeconv)
exp <- merkel.compendium$expr.list$merkel
dim(exp)

ann <- merkel.compendium$anno.list$merkel

sym <-
  mapg(rownames(exp))
length(rownames(exp)) == length(sym)
duinx <- which(duplicated(sym))
nainx <- which(is.na(sym))
inx <- unique(c(duinx, nainx))
if (length(inx) != 0)
  exp <- exp[-inx,]
rownames(exp) <- mapg(rownames(exp))

res.xcell = deconvolute(exp, "xcell")
res.ciber = deconvolute(exp, "cibersort_abs")
