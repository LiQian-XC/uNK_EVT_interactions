##featureCounts
library(Rsubread)
files = dir(path='mapping_res', pattern='*sortedByCoord.out.bam$')
res <- featureCounts(
		files = paste0('mapping_res/', files),
		annot.ext = '../genomes_annotations/gencode.v37.primary_assembly.annotation.gtf',
		isGTFAnnotationFile = TRUE,
		useMetaFeatures = TRUE,
		allowMultiOverlap = FALSE,
		largestOverlap = FALSE,
		ignoreDup = FALSE,
		countChimericFragments = FALSE,
		nthreads = 1,
		countMultiMappingReads = FALSE,
		strandSpecific = 0
		)

##rename colnames and rownames
gnm <- read.table('../genomes_annotations/gencode.v37.geneid.genename.txt', stringsAsFactors=F, sep='\t')
gnm[, 2] <- make.unique(gnm[, 2])
gene_len <- res$annotation$Length
names(gene_len) <- gnm[match(res$annotation$GeneID, gnm[, 1]), 2]
counts <- res$counts
colnames(counts) <- sub('.markdup.bam|Aligned.sortedByCoord.out.bam', '', colnames(counts))
colnames(counts) <- sub('^Y', '', colnames(counts))
colnames(counts) <- sub('(27_.*)P10', '\\1P11', colnames(counts))	
colnames(counts) <- sub('(27_.*)P9', '\\1P10', colnames(counts))	
colnames(counts) <- sub('(27_.*)P8', '\\1P9', colnames(counts))	
rownames(counts) <- gnm[match(rownames(counts), gnm[, 1]), 2]
cnms <- colnames(counts)
cp <- c('S1sp', 'S1L1dp', 'L1sp')
ids <- paste0('Y0', sub('_.*', '', cnms))
cds <- c('C2', 'C1')[as.integer(factor(sub('[^_]+_([^_]+)_.*', '\\1', cnms), levels=c('5', '8')))]
cps <- cp[as.integer(factor(sub('.*_', '', cnms), levels=paste0('P', 9:12)))]
colnames(counts) <- paste(ids, cps, cds, sep='_')
info <- data.frame(IDs = ids, CDs = cds, CPs = cps)
rownames(info) <- cnms

##filter genes with low read counts
ind <- apply(counts, 1, function(x) any(tapply(x, paste(cds, cps, sep='_'), function(y) sum(y > 5) > 2)))
counts <- counts[ind, ]

##calculate RPKM
rpkm <- t(t(counts)/colSums(counts))/gene_len[match(rownames(counts), names(gene_len))]*10^9
rpkm <- rpkm[, paste0(rep(paste0('Y0', c(15, 17, 18)), 6), '_', rep(cp, each=6), '_', rep(rep(c('C1', 'C2'), each=3), 3))]
ind <- apply(rpkm, 1, function(x) any(tapply(x, sub('(.*?)_', '', colnames(rpkm)), function(y) all(y > 1))))
rpkm <- rpkm[ind, ]

##PCA
pdf('PCA_for_supp.pdf', height=8, width=8)
par(mfcol=c(2, 2), xpd=T)
cols <- c('#bcbd22', '#d62728')
for(ty in cp){
    print(ty)
    mat <- rpkm[, grepl(ty, colnames(rpkm))]
    ind <- apply(mat, 1, function(x) any(tapply(x, sub('.*_', '', colnames(mat)), function(y) all(y > 1))))
    mat <- mat[ind, ]
    print(dim(mat))
    mat <- log(mat + 1)
    pca_res <- prcomp(t(mat), center = TRUE, scale. = TRUE)
    x1 <- pca_res$x[, 1]
    x2 <- pca_res$x[, 2]
    im1 <- round(summary(pca_res)$importance[2, 1]*100, 2)
    im2 <- round(summary(pca_res)$importance[2, 2]*100, 2)
    idx <- as.integer(factor(sub('.*_', '', colnames(mat)), levels=c('C1', 'C2')))
    plot(x1, x2, pch=c(22, 21)[idx], col='black', bg=cols[idx], cex=3, xlab=paste0('PC1 ', im1, '%'), ylab=paste0('PC2 ', im2, '%'), main=ty)
    legend(max(x1)*0.7, max(x2)*0.85, c('C1', 'C2'), pt.bg=cols, pch=c(22, 21), cex=1, pt.cex=1.2, bty='n')
}
dev.off()

##edgeR
library(edgeR)
de_res <- lapply(cp, function(cq) {
	print(cq)
	ind <- cps == cq
    sub_info <- info[ind, ]
    sub_count <- counts[, ind]
    sub_count <- sub_count[apply(sub_count, 1, function(x) any(tapply(x, sub_info$CDs, function(y) sum(y > 5) > 2))), ]
    y <- DGEList(counts=sub_count)
    y <- calcNormFactors(y)
    Subject <- sub_info$IDs
    Treat <- sub_info$CDs
    design <- model.matrix(~Subject + Treat)
    rownames(design) <- colnames(sub_count)
    yy <- estimateDisp(y, design)
    fit <- glmFit(yy, design)
    lrt <- glmLRT(fit)
    res <- topTags(lrt, n=nrow(yy))
    res <- as.data.frame(res)
    expr <- cpm(yy)
    expr <- expr[match(rownames(res), rownames(expr)), ]
    fres <- cbind(res, expr)
	fres
})
names(de_res) <- cp
saveRDS(de_res, file='edgeR_res.RDS')

##summarize differentially expressed genes between C2+ and C1+HLA-C for each uNK subset
genes <- list()
mat <- matrix(0, nrow=nrow(rpkm), ncol=3)
rownames(mat) <- rownames(rpkm)
colnames(mat) <- cp
for(i in 1:length(cp)){
	nm <- cp[i]
	dd <- de_res[[nm]]
    dd <- dd[rownames(dd) %in% rownames(rpkm), ]
	expr <- rpkm[match(rownames(dd), rownames(rpkm)), colnames(dd)[grepl('_', colnames(dd))]]
	ff <- sub('.*_', '', colnames(expr))
	ind <- apply(expr, 1, function(x) any(tapply(x, ff, function(y) all(y > 1))))
	dd <- dd[ind, ]
	fdr <- dd$FDR
	lfd <- dd$logFC
	gg <- rownames(dd)
	genes[[i]] <- gg

	ind <- fdr < 0.05 & abs(lfd) > log2(1.5)
	mat[match(gg[ind], rownames(mat)), i] <- sign(lfd[ind])
}

genes <- unique(unlist(genes))
mat <- mat[match(genes, rownames(mat)), ]
mat <- mat[order(apply(mat, 1, paste, collapse='_')), ]
mat[mat==0] <- 'Non-sig'
mat[mat==1] <- 'Up-reg'
mat[mat== -1] <- 'Down-reg'
mnm <- apply(mat, 1, function(stas) paste(paste0(cp, ':', stas), collapse=';'))
expr <- rpkm[match(rownames(mat), rownames(rpkm)), ]
plf <- lapply(cp, function(nm){
    dd <- de_res[[nm]][match(rownames(mat), rownames(de_res[[nm]])), c('FDR', 'logFC')]
    colnames(dd) <- paste0(nm, '_', colnames(dd))
    dd[is.na(dd[, 1]), 1] <- 1
    dd[is.na(dd[, 2]), 2] <- 0
    dd
})
plf <- do.call(cbind, plf)
eid <- sub('\\..*', '', gnm[match(rownames(mat), gnm[, 2]), 1])
ftab <- cbind(eid, mnm, plf, expr)
rownames(ftab) <- rownames(mat)
colnames(ftab)[1:2] <- c('Ensembl ID', 'DE_summary')
write.csv(ftab, file='table_all_genes_paired_comparisons.csv', quote=F, row.names=T, col.names=T)

##upset plot
library(stringr)
dd <- ftab[, 2]
mat <- matrix(0, nrow=2, ncol=7)
colnames(mat) <- 1:7
for(i in 1:length(cp)){
    cq <- cp[i]
    mat[1, i] <- sum((str_count(dd, 'Up-reg') == 1) & grepl(paste0(cq, ':Up-reg'), dd))
    mat[2, i] <- sum((str_count(dd, 'Down-reg') == 1) & grepl(paste0(cq, ':Down-reg'), dd))
    colnames(mat)[i] <- cq
}
tcps <- combn(cp, 2)
for(i in 1:ncol(tcps)){
    cp1 <- tcps[1, i]
    cp2 <- tcps[2, i]
    mat[1, i+3] <- sum(grepl(paste0(cp1, ':Up-reg'), dd) & grepl(paste0(cp2, ':Up-reg'), dd) & (str_count(dd, 'Up-reg') == 2))
    mat[2, i+3] <- sum(grepl(paste0(cp1, ':Down-reg'), dd) & grepl(paste0(cp2, ':Down-reg'), dd) & (str_count(dd, 'Down-reg') == 2))
    colnames(mat)[i+3] <- paste0(cp1, '&', cp2)
}
mat[1, 7] <- sum(grepl('S1sp:Up-reg', dd) & grepl('L1sp:Up-reg', dd) & grepl('S1L1dp:Up-reg', dd))
mat[2, 7] <- sum(grepl('S1sp:Down-reg', dd) & grepl('L1sp:Down-reg', dd) & grepl('S1L1dp:Down-reg', dd))
colnames(mat)[7] <- paste(cp, collapse='&')

pdf('DE_number_upset.pdf', height=4, width=5, useDingbats=F)
library(UpSetR)
num <- colSums(mat)
upset(fromExpression(num), keep.order=T)
dev.off()
