GWASpoly.fitted <- setClass("GWASpoly.fitted",
                            slots=c(scores="list",effects="list",gene_effec = "list",params="list"),contains="GWASpoly.K")

