library(tximport)
library(edgeR)

files <- file.path("salmon", list.files("salmon", "quant$"), "quant.sf")

txi <- tximport(files = files, type = "salmon", txOut = T)
cts <- txi$counts
colnames(cts) <- gsub("_1m_quant", "", list.files("salmon", "quant$"))
normMat <- txi$length
normMat <- normMat/exp(rowMeans(log(normMat)))

o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
y <- DGEList(cts)
y <- scaleOffset(y, t(t(log(normMat)) + o))
# filtering
keep <- filterByExpr(y)
y <- y[keep, ]
plotMDS(y, main = "MDS plot salmon quant files")
