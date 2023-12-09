id1 <- paste0('6044STDY864056', 1:6)
id2 <- paste0('Pla_Camb101239', 28:35)
res <- list()
for(id in list(id1, id2)){
	sres <- list()
	for(i in 1:(length(id)-1)){
		ssres <- list()
		for(j in (i+1):length(id)){
			dd <- read.table(paste0(id[i], '_', id[j], '_correspondence.txt'), stringsAsFactors=F, skip=6, nrows=3)[, 1:2]
			colnames(dd) <- id[c(i, j)]
			if(i == 1){
				dd <- dd[order(dd[, 1]), ]
			} else {
				dd <- dd[match(sres[[1]][, i], dd[, 1]), ]
			}
			ssres[[length(ssres) + 1]] <- dd 
		}
		ssres <- do.call(cbind, ssres)
		ssres <- ssres[, !duplicated(colnames(ssres))]
		sres[[length(sres) + 1]] <- ssres
	}
	res[[length(res) + 1]] <- sres
}
saveRDS(res, file='summary_souporcell_shared_samples.RDS')

info1 <- res[[1]][[1]]
rownames(info1) <- paste0('ID', 1:3)
info2 <- res[[2]][[1]]
rownames(info2) <- paste0('ID', 4:6)
write.csv(info1, file='souporcell_shared_samples_1.csv', quote=F, row.names=T, col.names=T)
write.csv(info2, file='souporcell_shared_samples_2.csv', quote=F, row.names=T, col.names=T)
