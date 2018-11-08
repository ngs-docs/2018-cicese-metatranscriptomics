library(ggplot2)
library(ggrepel)

comp <- read.csv("sourmash-compare/tara.comp.csv")

# Label the rows
rownames(comp) <- colnames(comp)

# Transform for plotting
comp <- as.matrix(comp)

fit <- dist(comp)
fit <- cmdscale(fit)
fit <- as.data.frame(fit)

fit$lab <- gsub("_5.20_rep[1,2]_1m.khmer.pe.fq.gz", "", rownames(fit))

ggplot(fit, aes(x = V1, y = V2)) +
  geom_point() + 
  geom_label_repel(label = fit$lab) + 
  theme_minimal() +
  ggtitle("MDS plot of sourmash compare on reads")
