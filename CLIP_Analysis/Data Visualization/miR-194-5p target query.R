library(tidyverse)
mir_peaks <- readRDS("Datafiles/miRNA-merged-peaks-list-12232019-withIDs.rds")

mir_194 <-  mir_peaks$`miR-194-5p` %>% as.data.frame()

write.csv(mir_194, "ce_HEAP_miR-194-5p_targets.csv", row.names = FALSE)
