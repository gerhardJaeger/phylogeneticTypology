library(ape)
library(geiger)

bi <- 0.5

files <- list.files("mrbayes/")

mbfiles <- files[grepl(".mb.nex$", files)]

families <- sub(".mb.nex*", "", mbfiles)

dir.create("../data/posteriorTrees", showWarnings = F, recursive = FALSE, mode = "0777")

for (fm in families) {
    print(fm)
    trees1 <- read.nexus(paste0("../data/asjpNex/output/", fm,'.run1.t'))
    trees2 <- read.nexus(paste0("../data/asjpNex/output/", fm,'.run2.t'))
    trees <- c(trees1[-(1:floor(bi*length(trees1)))],
               trees2[-(1:floor(bi*length(trees2)))])
    if (length(trees)>1000) {
        trees <- sample(trees,1000)
    }
    write.tree(trees,paste0("../data/posteriorTrees/", fm,'.posterior.tree'))
}


##

files = list.files("revbayes")
revfiles = files[grepl(".Rev", files)]
families = sub(".Rev","",revfiles)



for (fm in families) {
    print(fm)
    trees1S = as.character(read.table(paste0("revbayes/output/",fm,"_run_1.t"), header=T)$tree)
    trees1 = read.tree(text=trees1S)
    trees2S = as.character(read.table(paste0("revbayes/output/",fm,"_run_2.t"), header=T)$tree)
    trees2 = read.tree(text=trees2S)
    trees <- c(trees1[-(1:floor(bi*length(trees1)))],
               trees2[-(1:floor(bi*length(trees2)))])
    if (length(trees)>1000) {
        trees <- sample(trees,1000)
    }
    write.tree(trees,paste0("../data/posteriorTrees/", fm,'.posterior.tree'))
}
