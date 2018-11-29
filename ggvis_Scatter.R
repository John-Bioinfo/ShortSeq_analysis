

setwd("D:\\python_Code\\BlastXML\\grade_test")
library(ggvis)

DataG <- read.table("test_MG.txt", sep="\t", header=TRUE)
colnames(DataG) <- c("Grade", "PN", "L", "seq","T005","T1", "GC", "BD", "DG","GB","MG")
DataG$Grade <- as.factor(DataG$Grade)
DataG %>% ggvis(~L, ~GB, fill=~Grade) %>% layer_points()
