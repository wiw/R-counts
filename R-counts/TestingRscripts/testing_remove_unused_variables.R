rm(list=ls())
library(ggplot2)
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), ncol = cols, nrow = ceiling(numPlots/cols))
  }
 if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
    }
  }
}
scatterList <- list()
ic <- c("multiplot", "scatterList", "ic")
for (i in c("mpg", "disp")) {
	for (y in c("hp", "qsec")) {
		scattername <- paste(i, y, sep = "_")
		scatterList[[scattername]] <- ggplot(mtcars, aes(mtcars[[i]], mtcars[[y]])) + geom_point(alpha = 1/10, colour = "red", size = 3, na.rm = T) + annotate("text", x = c(max(mtcars[[i]], na.rm = T) * 0.5, max(mtcars[[i]], na.rm = T) * 0.5), y = c(max(mtcars[[y]], na.rm = T) * 0.9, max(mtcars[[y]], na.rm = T) * 0.85), label = c(paste("Pearson.Cor = ", round(cor(mtcars[[i]], mtcars[[y]], method="pearson", use="pairwise.complete.obs"), digits=2), sep = ""), paste("Spearman.Cor = ", round(cor(mtcars[[i]], mtcars[[y]], method="spearman", use="pairwise.complete.obs"), digits=2), sep = "")), size=4)
	}
}
save(list=ic, file="mtcars_test.RData")
rm(list=ls()[!(ls() %in% ic)])
print(multiplot(plotlist = scatterList, cols = 2, layout = matrix(c(1:4), ncol = 2, nrow = 2, byrow = T)))
