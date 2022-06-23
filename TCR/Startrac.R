library("Startrac")
library("tictoc")
load("TCR_all_startrac.RData")
obj.proj <- new("Startrac", TCR_all_startrac, aid = "NPC_cd8")
patient.vec <- unique(TCR_all_startrac$patient)
obj.list <- list()

for (i in c(patient.vec)) {
  require("Startrac")
  obj <- new("Startrac", subset(obj.proj@cell.data, patient ==
                                  i), aid = i)
  obj <- calIndex(obj)
  # obj <- pIndex(obj, cores = 1)
  obj.list <- append(obj.list,list(obj@cluster.data))
  print(i)
}


obj <-  do.call(rbind, obj.list)
obj$Diversity <- 1-obj$expa

obj$Diversity[is.na(obj$Diversity)] <- 0
a <- obj$Diversity[obj$majorCluster == "CD8_C9_HAVCR2"] 
b <- obj$Diversity[obj$majorCluster == "CD8_C8_CXCL13"] 
wilcox.test(a,b)

png("cd8_diver.png",width = 5, height =5, res = 400, units = "in")
ggboxplot(obj,
          x="majorCluster",y="Diversity",
          color = "#6A51A3", add = "point", outlier.colour=NULL) +
  theme(axis.text.x=element_text(angle = 45,hjust = 1))
dev.off()

