library(limma)
library(edgeR)
library(magrittr)
library(NMF)

##read in counts, alter/fix column names. 
fcc <- read.table("tolerance.counts.fc",stringsAsFactors = FALSE,header = TRUE)
	colnames(fcc) <- gsub("bowtie.Sample_","",colnames(fcc))
	colnames(fcc) <- gsub(".out.accepted_hits.bam","",colnames(fcc))
	colnames(fcc) <- gsub("2000","200",colnames(fcc))
	colnames(fcc) <- gsub("KC.F2","KC.F1",colnames(fcc))
	colnames(fcc) <- gsub("ER.F1","ER.F2",colnames(fcc))

##separate counts from metadata
fcc2 <- as.matrix(fcc[,7:81])

##table of factors. process them. 
tags <- read.table("../bowtie.bams.list")
	tags <- gsub("bowtie/Sample_","",tags[,1])
	tags <- gsub(".out/accepted_hits.bam","",tags)
	tags <- do.call(rbind,strsplit(x = tags,split = "-"))
	tags[20:30,2] <- "F2"
	tags[41:49,2] <- "F1"
	tags[25:27,3] <- "200"
	tags <- cbind(tags,c(rep(1,9),rep(2,10),rep(4,11),rep(2,10),rep(4,9),rep(1,8),rep(3,8),rep(3,10)))
	tags <- cbind(tags,c(rep(1,19),rep(2,11),rep(1,10),rep(2,9),rep(1,26)))
	colnames(tags) <- c("population","tolerance","dose","replicate","location","region")
	facs <- paste(tags[,1],tags[,2],tags[,3],tags[,5],tags[,6],sep = ".")

##remove samples with various issues. Weird outlier ER sample should be checked.

	###remove three lowest coverage samples. two ER undosed and one dosed BP. 

	#lowcov <- c(19:21)
	#cat("removing ", colnames(fcc2)[lowcov],"\n")
	#fcc <- fcc[,-lowcov]
	#fcc2 <- fcc2[,-lowcov]
	#tags <- tags[-lowcov,]
	#facs <- facs[-lowcov]

	
	###remove the single dropout individual. one undosed ER individual. 
	dropout <- which(colnames(fcc2)=="ER.F2.0.1")
	fcc <- fcc[,-(20+6)]
	fcc2 <- fcc2[,-20]
	tags <- tags[-20,]
	facs <- facs[-20]
	
	###test out removing outliers as well. 
	undosedhigh <- which(colnames(fcc2)%in%c("BI.F1.0.3","BP.F2.0.3"))
	fcc <- fcc[,-(undosedhigh+6)]
	fcc2 <- fcc2[,-undosedhigh]
	tags <- tags[-undosedhigh,]
	facs <- facs[-undosedhigh]

	#toleranthigh <- "BP.F2.200.1"
	#toleranthigh <- which(colnames(fcc2)%in%toleranthigh)
	
	#remove tolerant high	
	#fcc <- fcc[,-toleranthigh]
	#fcc2 <- fcc2[,-toleranthigh]
	#tags <- tags[-toleranthigh,]
	#facs <- facs[-toleranthigh]

	##remove weird outlier
	##figure out what's going on with this guy?
	undosedhigh <- which(colnames(fcc2)%in%"ER.F2.0.4")
	fcc <- fcc[,-(undosedhigh+6)]
	fcc2 <- fcc2[,-undosedhigh]
	tags <- tags[-undosedhigh,]
	facs <- facs[-undosedhigh]

##image counts
image(log(fcc2[order(rowSums(fcc2),decreasing = T),]))

##create annotation matrix
annot <- fcc[,1:6]
annot[,2] <- gsub(";.*","",annot[,2])
annot <- annot[,-(3:5)]
	
#create design matrix
TS <- factor(paste(tags[,1],tags[,2],tags[,3],tags[,6],sep = "."))
design <- model.matrix(~0+TS)
colnames(design) <- levels(TS)

###start analysis. turn counts into DGEList object. 
fcc3 <- DGEList(counts = fcc2, genes = annot, group = TS)

##toss out some low coverage genes
keep <- rowSums(cpm(fcc3))>20
fcc3 <- fcc3[keep,,keep.lib.sizes=FALSE]
dim(fcc3)

##time to figure out what the hell this is doing
##TMM normalization. creates sample weights as object in fcc3 list. 
fcc3 <- calcNormFactors(fcc3)
image(log(fcc3$counts[order(rowSums(fcc3$counts),decreasing = T),]))

##MDS plot of samples for 20 most variable genes
plotMDS(fcc3,labels = colnames(fcc2),top = 20, col = as.numeric(factor(tags[,3])), gene.selection = "common",prior.count = 5)

##estimate GLM dispersions
fcc3 <- estimateGLMCommonDisp(fcc3, design, verbose=TRUE)
fcc3 <- estimateGLMTrendedDisp(fcc3, design)
fcc3 <- estimateGLMTagwiseDisp(fcc3, design)

##fit glm model
fit <- glmFit(fcc3,design)

##contrasts to be evaluated.
cont.matrix <- makeContrasts(
    do_x_to = (((BI.F1.0.1-BI.F1.200.1)+(FLAX.F1.0.1-FLAX.F1.200.1)+(SH.F1.0.1-SH.F1.200.1)+(KC.F1.0.2-KC.F1.200.2))/4)-(((NBH.F2.0.1-NBH.F2.200.1)+(BP.F2.0.1-BP.F2.200.1)+(NWK.F2.0.1-NWK.F2.200.1)+(ER.F2.0.2-ER.F2.200.2))/4),
    do_x_to_x_re = ((((BI.F1.0.1-BI.F1.200.1)+(FLAX.F1.0.1-FLAX.F1.200.1)+(SH.F1.0.1-SH.F1.200.1))/3)-(((NBH.F2.0.1-NBH.F2.200.1)+(BP.F2.0.1-BP.F2.200.1)+(NWK.F2.0.1-NWK.F2.200.1))/3))-((KC.F1.0.2-KC.F1.200.2)-(ER.F2.0.2-ER.F2.200.2)),
    Ndo_x_tol = ((((BI.F1.0.1-BI.F1.200.1)+(FLAX.F1.0.1-FLAX.F1.200.1)+(SH.F1.0.1-SH.F1.200.1))/3)-(((NBH.F2.0.1-NBH.F2.200.1)+(BP.F2.0.1-BP.F2.200.1)+(NWK.F2.0.1-NWK.F2.200.1))/3)),
    Sdo_x_tol = (KC.F1.0.2-KC.F1.200.2)-(ER.F2.0.2-ER.F2.200.2),
	Sendo = (((BI.F1.0.1-BI.F1.200.1)+(FLAX.F1.0.1-FLAX.F1.200.1)+(SH.F1.0.1-SH.F1.200.1)+(KC.F1.0.2-KC.F1.200.2))/4),
	Toldo = (((NBH.F2.0.1-NBH.F2.200.1)+(BP.F2.0.1-BP.F2.200.1)+(NWK.F2.0.1-NWK.F2.200.1)+(ER.F2.0.2-ER.F2.200.2))/4),
# 	do_x_to_x_re2 = ((((BI.F1.0.1-BI.F1.200.1)+(FLAX.F1.0.1-FLAX.F1.200.1))/2)-(((NBH.F2.0.1-NBH.F2.200.1)+(BP.F2.0.1-BP.F2.200.1))/2))-((KC.F1.0.2-KC.F1.200.2)-(ER.F2.0.2-ER.F2.200.2)),
#	N1_do_x_to = (BI.F1.0.1-BI.F1.200.1)-(NBH.F2.0.1-NBH.F2.200.1),
#	N2_do_x_to = (FLAX.F1.0.1-FLAX.F1.200.1)-(BP.F2.0.1-BP.F2.200.1),
#	N3_do_x_to = (SH.F1.0.1-SH.F1.200.1)-(NWK.F2.0.1-NWK.F2.200.1),
#	S1_do_x_to = (KC.F1.0.2-KC.F1.200.2)-(ER.F2.0.2-ER.F2.200.2),
#	N1_do_x_to = (BI.F1.0.1-BI.F1.200.1)-(NBH.F2.0.1-NBH.F2.200.1)-((((FLAX.F1.0.1-FLAX.F1.200.1)-(BP.F2.0.1-BP.F2.200.1))+((SH.F1.0.1-SH.F1.200.1)-(NWK.F2.0.1-NWK.F2.200.1))+((KC.F1.0.2-KC.F1.200.2)-(ER.F2.0.2-ER.F2.200.2)))/3),
#	N2_do_x_to = (FLAX.F1.0.1-FLAX.F1.200.1)-(BP.F2.0.1-BP.F2.200.1)-((((BI.F1.0.1-BI.F1.200.1)-(NBH.F2.0.1-NBH.F2.200.1))+((SH.F1.0.1-SH.F1.200.1)-(NWK.F2.0.1-NWK.F2.200.1))+((KC.F1.0.2-KC.F1.200.2)-(ER.F2.0.2-ER.F2.200.2)))/3),
#	N3_do_x_to = (SH.F1.0.1-SH.F1.200.1)-(NWK.F2.0.1-NWK.F2.200.1)-((((BI.F1.0.1-BI.F1.200.1)-(NBH.F2.0.1-NBH.F2.200.1))+((FLAX.F1.0.1-FLAX.F1.200.1)-(BP.F2.0.1-BP.F2.200.1))+((KC.F1.0.2-KC.F1.200.2)-(ER.F2.0.2-ER.F2.200.2)))/3),
#	S1_do_x_to = (KC.F1.0.2-KC.F1.200.2)-(ER.F2.0.2-ER.F2.200.2)-((((BI.F1.0.1-BI.F1.200.1)-(NBH.F2.0.1-NBH.F2.200.1))+((FLAX.F1.0.1-FLAX.F1.200.1)-(BP.F2.0.1-BP.F2.200.1))+((SH.F1.0.1-SH.F1.200.1)-(NWK.F2.0.1-NWK.F2.200.1)))/3),
	utol = ((BI.F1.0.1+FLAX.F1.0.1+SH.F1.0.1+KC.F1.0.2)/4)-((NBH.F2.0.1+BP.F2.0.1+NWK.F2.0.1+ER.F2.0.2)/4),
	Nu_tol = ((BI.F1.0.1+FLAX.F1.0.1+SH.F1.0.1)/3)-((NBH.F2.0.1+BP.F2.0.1+NWK.F2.0.1)/3),
	Su_tol = KC.F1.0.2-ER.F2.0.2,
	utol_re = (((BI.F1.0.1+FLAX.F1.0.1+SH.F1.0.1)/3)-((NBH.F2.0.1+BP.F2.0.1+NWK.F2.0.1)/3))-(KC.F1.0.2-ER.F2.0.2),
	levels = design
 	)

##evaluate a single contrast, examine results. 
lrt <- glmLRT(fit, contrast=cont.matrix[,5])
topTags(lrt)
results <- decideTestsDGE(lrt)
detags <- rownames(fcc3)[as.logical(results)]
plotSmear(lrt, de.tags=detags)
abline(h=c(-1, 1), col="blue")

subC <- cpm(fcc3$counts[detags,],log=TRUE,prior.count=1)
aheatmap(subC)
plotMDS(subC,
	labels = colnames(fcc2),top = 400, 
	col = as.numeric(factor(tags[,3])), 
	gene.selection = "common")

##evaluate all contrasts. store results. 
ncont <- 10
lrt <- list()
for(i in 1:ncont){
	lrt[[i]] <- glmLRT(fit, contrast=cont.matrix[,i])
	}

outvals <- list()
outvals[["FDR"]] <- numeric(0)
outvals[["logFC"]] <- numeric(0)
outvals[["logCPM"]] <- numeric(0)

for(i in 1:ncont){
	tmp <- topTags(lrt[[i]],sort.by="none",n=24441)
	outvals[["FDR"]] <- cbind(outvals[["FDR"]],tmp$table$FDR)
	outvals[["logFC"]] <- cbind(outvals[["logFC"]],tmp$table$logFC)
	outvals[["logCPM"]] <- cbind(outvals[["logCPM"]],tmp$table$logCPM)
	}

outvals[["sig"]] <- rowSums(outvals[["FDR"]]<0.05)>0

resum <- table(apply(outvals[["FDR"]]<0.05,MAR=1,FUN=function(x){paste(as.numeric(x),collapse=".")}))
resum <- cbind(do.call(rbind,strsplit(split="\\.",names(resum))),as.numeric(resum))
class(resum)<-"numeric"
resum<-resum[order(resum[,ncont+1],decreasing=TRUE),]
colSums(outvals[["FDR"]]<0.03)



##plot some results. 

subC <- cpm(fcc3,log=TRUE,prior.count=1)

subi <- !grepl("200",colnames(fcc3$counts))
#subg <- which(outvals[["FDR"]][,7]<0.05&rowSums(outvals[["FDR"]]<0.05)==1)
subg <- which(outvals[["FDR"]][,7]<0.03)
#subg <- which((outvals[["FDR"]][,1]<0.05)&(outvals[["FDR"]][,3]>0.05))

plotMDS(fcc3$counts[subg,subi],
labels = colnames(fcc2)[subi],top = 20, 
col = as.numeric(factor(tags[,3][subi])), 
gene.selection = "common")

aheatmap(subC[subg,subi])

#fold change from CPM average
subC2 <- subC[subg,subi]-outvals[["logCPM"]][subg,1]
subC2[subC2>2] <- 2
subC2[subC2<(-2)] <- -2

ann <- tags[subi,c(2,3,6)]
aheatmap(subC2,annCol=ann)


vennDiagram(outvals[["FDR"]][,c(1,3,7)]<0.05)


