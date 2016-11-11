########################################################################
# A script to analyse intron-retention with chromatin marks in gut development
########################################################################

require(gdata)
require(GenomicRanges)	# for overlap calculations
require(rtracklayer)	# to parse .bed files as GRanges
require(ggplot2)		# for nice plots
require(BSgenome.Mmusculus.UCSC.mm9) # for proper seqinfo object
require(RColorBrewer)	# to choose nice colors
require(plyr)			# for transformation and summary of data.frames

#-----------------------------------------------------------------------
# Define some parameters and file names
#-----------------------------------------------------------------------
COL_CELL_TYPE=brewer.pal(8, "Pastel1")
COL_ALL=brewer.pal(8, "Dark2")[8]

CELL_TYPE_LEVELS <- c("e125", "ISC", "Paneth", "enterocyte")
MARK_LEVELS <- c("k4me3", "k27ac")
SCORE_TH=0.1

VERSION="v04"

# design / dimensions

# mark
#	cell_type
#		 time
#			replicate
# 
#/data/Chromatin/k4me3_peaks/K4me3_called_peaks_e125_1.csv

chromMarksDesignFile <- "data/chromatin_marks_design.csv"
#~ intronDesignFile <- "data/retained_intron_design.csv"
#~ exonDesignFile <- "data/skipped_exon_design.csv"
designFile <- "data/splicing_design.csv"

allIntronFiles <- "data/Intron_Retention_bed_files/Intron_Retention_Annotated/intron_retained_in_all_genome_direction_annotated_miso_ordered.bed.mm10.BED.mm9"

allExonsFiles <- "data/Events_With_Signals/Jonas_Exon_Skipping/SE.events.annotated.mm9.bed"

outPrefix <- paste0("results/", VERSION, "_GutDev")
dir.create(dirname(outPrefix), recursive=TRUE, showWarnings = FALSE)

#-----------------------------------------------------------------------
# Parse input data
#-----------------------------------------------------------------------

# read design and meta-data of marks
chromMarkDesign <- read.table(chromMarksDesignFile, header=TRUE)
# take only a subest with existing RNA-seq data
chromMarkDesign$cell_type <- factor(chromMarkDesign$cell_type, levels=CELL_TYPE_LEVELS)
chromMarkDesign$mark <- factor(chromMarkDesign$mark, MARK_LEVELS)

# read meta-data of introns and exons
#~ intronDesign <- read.table(intronDesignFile, header=TRUE)
#~ intronDesign$cell_type <- factor(intronDesign$cell_type, levels=CELL_TYPE_LEVELS)

#~ exonDesign <- read.table(exonDesignFile, header=TRUE)
#~ exonDesign$cell_type <- factor(exonDesign$cell_type, levels=CELL_TYPE_LEVELS)

design <- read.table(designFile, header=TRUE)
design$cell_type <- factor(design$cell_type, levels=CELL_TYPE_LEVELS)

########################################################################
# DEBUG: take only exons into account
design <- design[which(design$type == "exon"),]
########################################################################

# get seqinfo object from annotation package
genome <- BSgenome.Mmusculus.UCSC.mm9
seqInfo <- seqinfo(genome)

# parse chromatin mark peaks as GenomicRanges objects
markGRL <- lapply(chromMarkDesign$file, function(f){
	gr <- GRanges(read.delim(as.character(f), header=TRUE), seqinfo=seqInfo)
})

# parse all annotated introns and exons
allIntrons <- import(allIntronFiles, format="bed", seqinfo=seqInfo)
allExons <- import(allExonsFiles, format="bed", seqinfo=seqInfo)

# make list of all inrons and exon
allReg = List(
	intron=allIntrons,
	exon=allExons
)

# parse all retained introns and skipped exons
regGRL <- lapply(design$file, function(f){
	df <- read.delim(as.character(f), header=FALSE)
	gr <- GRanges(df[,1], IRanges(df[,2], df[,3]), score=df[,4])
})

regGRL <- GRangesList(regGRL)

# add number of regions to design table
design$n <- sapply(regGRL, length)

#~ # parse all skipped exons
#~ exonGRL <- lapply(exonDesign$file, function(f){
#~ 	df <- read.delim(as.character(f), header=FALSE)
#~ 	gr <- GRanges(df[,1], IRanges(df[,2], df[,3]), PSI=df[,4])
#~ })
#~ exonGRL <- GRangesList(exonGRL)


#-----------------------------------------------------------------------
# Compare length of introns and exons
#-----------------------------------------------------------------------

# create a data.frame with all intron sizes from retained introns
widthDF <- do.call("rbind", lapply(1:nrow(design), function(i){
  data.frame(design[rep(i, length(regGRL[[i]])),], width=width(regGRL[[i]]))
}))

# add all introns
allIntronsWidthDF <- data.frame(
    name="All", 
    cell_type="all",
    type="intron",
    RNAseq_replicate=1,
    file=allIntronFiles,
    n=length(allIntrons),
    width=width(allIntrons)
  )


# add all exons
allExonsWidthDF <- data.frame(
    name="All", 
    cell_type="all",
    type="exon",
    RNAseq_replicate=1,
    file=allExonsFiles,
    n=length(allExonsFiles),
    width=width(allExons)
)
  
    
widthDF <- rbind(
    allIntronsWidthDF,
    allExonsWidthDF,
    widthDF
)

#-----------------------------------------------------------------------
# plot both
g <- ggplot(widthDF, aes(x=as.factor(RNAseq_replicate), y=width, fill=cell_type)) + 
    geom_boxplot()+
    facet_grid(type~cell_type, scales="free_x") +
    theme_bw() + scale_fill_manual(values=c(COL_ALL, COL_CELL_TYPE)) + 
    xlab("RNA-seq replicates") + ylab("Size [bp]")
ggsave(paste0(outPrefix, ".region_size.boxplot.pdf"))

g +  scale_y_log10() 
ggsave(paste0(outPrefix, ".region_sizes.boxplot_log10.pdf"))


#-----------------------------------------------------------------------
# plot only introns
intronWidthDF <- subset(widthDF, type == "intron")

g <- ggplot(intronWidthDF, aes(x=as.factor(RNAseq_replicate), y=width, fill=cell_type)) + 
    geom_boxplot()+
    facet_grid(.~cell_type, scales="free_x") +
    theme_bw() + scale_fill_manual(values=c(COL_ALL, COL_CELL_TYPE)) + 
    xlab("RNA-seq replicates") + ylab("Intron size [bp]")
ggsave(paste0(outPrefix, ".intron_sizes.boxplot.pdf"))

g +  scale_y_log10() 
ggsave(paste0(outPrefix, ".intron_sizes.boxplot_log10.pdf"))


#-----------------------------------------------------------------------
# plot only exons
exonWidthDF <- subset(widthDF, type == "exon")

g <- ggplot(exonWidthDF, aes(x=as.factor(RNAseq_replicate), y=width, fill=cell_type)) + 
    geom_boxplot()+
    facet_grid(.~cell_type, scales="free_x") +
    theme_bw() + scale_fill_manual(values=c(COL_ALL, COL_CELL_TYPE)) + 
    xlab("RNA-seq replicates") + ylab("Exon size [bp]")
ggsave(paste0(outPrefix, ".exon_sizes.boxplot.pdf"))

g +  scale_y_log10() 
ggsave(paste0(outPrefix, ".exon_sizes.boxplot_log10.pdf"))


#-----------------------------------------------------------------------
# Plot number of retained introns and skipped exons per cell and across replicates
#-----------------------------------------------------------------------

g <- ggplot(design, aes(x=RNAseq_replicate, y=n, fill=cell_type)) + 
  geom_bar(stat="identity", color="black") + 
  geom_text(aes(label=n), vjust=-1.25) + 
  facet_grid(type~cell_type) + 
  theme_bw() + scale_fill_manual(values=COL_CELL_TYPE) + 
  xlab("RNA-seq replicates") + ylab("Events")

ggsave(paste0(outPrefix, ".regions.barplot.pdf"))

#-----------------------------------------------------------------------
# plot only introns
g <- ggplot(subset(design, type=="intron"), aes(x=RNAseq_replicate, y=n, fill=cell_type)) + 
  geom_bar(stat="identity", color="black") + 
  geom_text(aes(label=n), vjust=-1.25) + 
  facet_grid(~cell_type) + 
  theme_bw() + scale_fill_manual(values=COL_CELL_TYPE) + 
  xlab("RNA-seq replicates") + ylab("Retained Introns")

ggsave(paste0(outPrefix, ".retained_introns.barplot.pdf"))


#-----------------------------------------------------------------------
# plot only exons

g <- ggplot(subset(design, type=="exon"), aes(x=RNAseq_replicate, y=n, fill=cell_type)) + 
  geom_bar(stat="identity", color="black") + 
  geom_text(aes(label=n), vjust=-1.25) + 
  facet_grid(~cell_type) + 
  theme_bw() + scale_fill_manual(values=COL_CELL_TYPE) + 
  xlab("RNA-seq replicates") + ylab("Skipped exons")

ggsave(paste0(outPrefix, ".skipped_exons.barplot.pdf"))

#-----------------------------------------------------------------------
# Plot number and size distribution of chromatin mark peaks
#-----------------------------------------------------------------------

# create a data.frame with all chrmatin peaks
peakDF <- do.call("rbind", lapply(1:nrow(chromMarkDesign), function(i){
  data.frame(chromMarkDesign[rep(i, length(markGRL[[i]])),], width=width(markGRL[[i]]))
}))

# reformat as factor with order
#~ peakDF$mark = factor(peakDF$mark, c("k4me3", "k27ac"))

# count peaks per condition
countPeakDF <- ddply(peakDF, c("mark", "cell_type", "replicate"), summarise, 
                     count=length(width))

# plot number of peaks
g <- ggplot(countPeakDF, aes(x=as.factor(replicate), y=count, fill=cell_type)) + 
  geom_bar(stat="identity", color="black") + 
  geom_text(aes(label=count), vjust=1.25, size=4) + 
  facet_grid(mark~cell_type) + 
  theme_bw() + scale_fill_manual(values=COL_CELL_TYPE) + 
  xlab("ChIP-seq replicates") + ylab("Number of chromatin mark peaks")

ggsave(paste0(outPrefix, ".mark_peaks.barplot.pdf"))

# plot length distribution of chromatin marks
g <- ggplot(peakDF, aes(x=as.factor(replicate), y=width, fill=cell_type)) + 
  geom_boxplot() + scale_y_log10() +
  facet_grid(mark~cell_type, scales="free_x") +
  theme_bw() + scale_fill_manual(values=COL_CELL_TYPE) + 
  xlab("ChIP-seq replicates") + ylab("Chromatin mark peaks size [bp]")
ggsave(paste0(outPrefix, ".mark_peak.size.boxplot.pdf"))  
  
#-----------------------------------------------------------------------
# calculate overlaps
#-----------------------------------------------------------------------

# initialize data.frame
allDF <- data.frame()
allScoreDF <- data.frame()

# iterate over chrom mark types
for (MARK in MARK_LEVELS){
	
	# iteratre over cell types
	for (CELL_TYPE in CELL_TYPE_LEVELS){
		
		# iterate over RNA-seq replicates
		for (RNA_REP in unique(design$RNAseq_replicate)) {

			# iterate over type (retained intron / skipped exon )
#~ 			for (REG_TYPE in levels(design$type)) {
			for (REG_TYPE in c("exon")){
				
				# choose index of region
				i <- which(
					design$cell_type == CELL_TYPE & 
					design$RNAseq_replicate == RNA_REP &
					design$type == REG_TYPE)
				
				# iterate chrom-mark ChIP-seq replicates
				for (MARK_REP in unique(chromMarkDesign$replicate)) {
					
					# choose index of marks
					m <- which(chromMarkDesign$mark == MARK 
						& chromMarkDesign$cell_type == CELL_TYPE
						& chromMarkDesign$replicate == MARK_REP)
					
					# debug message
					message(paste("INFO:", MARK, CELL_TYPE, RNA_REP, MARK_REP, i, m))
					
					# compute overlap with marks
					ol <- countOverlaps(regGRL[[i]], markGRL[[m]], ignore.strand=TRUE) > 0
										
					# construct data frame with all regions as rows and column if overlaps
					scoreDF <- data.frame(
						mark=MARK, 
						cell_type=CELL_TYPE, 
						RNAseq_replicate=RNA_REP, 
						mark_replicate=MARK_REP, 
						type=REG_TYPE,
						mark_overlap=ol,
						n=sum(ol), 
						percent=100* sum(ol) / length(ol),
						score=as.numeric(score(regGRL[[i]]))
						)
					
					# add data fame to allScoreDF
					allScoreDF <- rbind(allScoreDF, scoreDF)
					
					
				} 
			}
		}
	}
}


#allScoreDF$cell_type <- factor(allScoreDF$cell_type, CELL_TYPE_LEVELS)

# Run the functions length, mean, and sd on the value of "score" for each group, 
sumData <- ddply(allScoreDF, c("type", "mark", "mark_replicate", "cell_type", "RNAseq_replicate", "mark_overlap"), summarise,
   N    = length(score),
   mean = mean(score),
   median = median(score),
   sd   = sd(score),
   se   = sd / sqrt(N)
)


# calculate p-values
pvalDF <- ddply(allScoreDF, c("type", "mark", "mark_replicate", "cell_type", "RNAseq_replicate"), summarise, pval  = wilcox.test(score~mark_overlap)$p.value
)

g <- ggplot(allScoreDF, aes(x=factor(mark_overlap, c(TRUE, FALSE), c("overlap", "no overlap")), y=score, fill=cell_type))+
	geom_boxplot() + #geom_jitter() +
	facet_grid(mark*mark_replicate~type*cell_type*RNAseq_replicate) + 
	geom_text(data=pvalDF, aes(x=1.5, y=1.1, label=paste0("p=", signif(pval,3))), size=2) +
	geom_text(data=sumData, aes(y=-0.1, label=paste0("n=", N)), size=2) +
	theme_bw() + theme(text=element_text(size=10), axis.text.x = element_text(angle = 60, hjust = 1)) + 
	scale_fill_manual(values=COL_CELL_TYPE) + 
	xlab("Overlap with mark") + ylab("Score (PSI/PRI)")

ggsave(paste0(outPrefix, ".regions.score_by_overlap.boxplot.pdf"), w=10.5, h=7)

#-----------------------------------------------------------------------
# group introns/exons in score by retained/skipped
#-----------------------------------------------------------------------
allScoreDF$event <- allScoreDF$score <= SCORE_TH
allScoreDF$event <- factor(allScoreDF$event, c(TRUE, FALSE), c("excluded", "included"))

# calculate p-values by fisher test on event groups
pvalDF <- ddply(allScoreDF, c("type", "mark", "mark_replicate", "cell_type", "RNAseq_replicate"), summarise, pval  = fisher.test(event,mark_overlap)$p.value
)

summaryDF <- ddply(allScoreDF, c("type", "mark", "mark_replicate", "cell_type", "RNAseq_replicate", "event"), summarise, 
   N    = length(n),
   mean = mean(n),
   sd   = sd(n),
   se   = sd / sqrt(N),
   percent_mean = mean(percent),
   percent_sd   = sd(percent),
   percent_se   = percent_sd / sqrt(N),
   mean_p = mean(p_val)
)



g <- ggplot(allScoreDF, aes(x=event, y=mark_overlap, fill=cell_type)) + 
	geom_bar(color="black") + 
	geom_text(aes(label=signif(percent_mean,3)), vjust=-1.25) + 
  geom_text(aes(x=1.5, label=ifelse(is.na(mean_p), "", paste0("p=", signif(mean_p,3)))), vjust=-2.5) + 
  #~ 	geom_text(aes(label=paste("n=",round(mean,1)), y=0), angle=90) + 
	ylim(0,1.4*max(cdata$percent_mean, na.rm=TRUE)) +
	facet_grid(mark*mark_replicate~cell_type) + 
	theme_bw() + scale_fill_manual(values=COL_CELL_TYPE) + 
	xlab("Introns") + ylab("Introns overlapping chromatin mark [%]")
	
ggsave(paste0(outPrefix, ".introns_overlapping_marks.barplot_percent.pdf"))

#----------



# Run the functions length, mean, and sd on the value of "score" for each group, 
sumData <- ddply(allScoreDF, c("type", "mark", "mark_replicate", "cell_type", "RNAseq_replicate", "mark_overlap"), summarise,
   N    = length(score),
   mean = mean(score),
   median = median(score),
   sd   = sd(score),
   se   = sd / sqrt(N),
)


# calculate p-values
pvalDF <- ddply(allScoreDF, c("type", "mark", "mark_replicate", "cell_type", "RNAseq_replicate"), summarise, pval  = wilcox.test(score~mark_overlap)$p.value
)



#~ # plot length distribution of chromatin marks
#~ g <- ggplot(peakDF, aes(x=as.factor(replicate), y=width, fill=cell_type)) + 
#~   geom_boxplot() + scale_y_log10() +
#~   facet_grid(mark~cell_type, scales="free_x") +
#~   theme_bw() + scale_fill_manual(values=COL_CELL_TYPE) + 
#~   xlab("ChIP-seq replicates") + ylab("Chromatin mark peaks size [bp]")
#~ ggsave(paste0(outPrefix, ".mark_peak.size.boxplot.pdf"))  
  	
#~ g <- ggplot(cdata, aes(x=group, y=percent_mean, fill=cell_type)) + 
#~ 	geom_bar(stat="identity", color="black") + 
#~ 	geom_errorbar(aes(ymax = percent_mean + percent_sd, ymin=percent_mean - percent_sd), width=0.25) + 
#~ 	geom_text(aes(label=signif(percent_mean,3)), vjust=-1.25) + 
#~   geom_text(aes(x=1.5, label=ifelse(is.na(mean_p), "", paste0("p=", signif(mean_p,3)))), vjust=-2.5) + 
#~   	#geom_text(aes(label=paste("n=",round(mean,1)), y=0), angle=90) + 
#~ 	ylim(0,1.4*max(cdata$percent_mean, na.rm=TRUE)) +
#~ 	facet_grid(mark*mark_replicate~cell_type) + 
#~ 	theme_bw() + scale_fill_manual(values=COL_CELL_TYPE) + 
#~ 	xlab("Introns") + ylab("Introns overlapping chromatin mark [%]")
	

					#---------------------------------------------------
					# split regions in two groups
					#---------------------------------------------------
					
					olAll <- countOverlaps(allReg[[REG_TYPE]], markGRL[[m]], ignore.strand=TRUE) > 0
					# 
					#f.test <- fisher.test(ol, olAll)
					
					olFactor <- factor(as.character(ol), levels=c("TRUE", "FALSE"))
					olAllFactor <- factor(as.character(olAll), levels=c("TRUE", "FALSE"))
	
					t <- cbind(table(olFactor), table(olAllFactor))
					ft <- fisher.test(t)
					p <- ft$p.value
					
					# add info to allDF 
					row1 <- data.frame(
						mark=MARK, 
						cell_type=CELL_TYPE, 
						RNAseq_replicate=RNA_REP, 
						mark_replicate=MARK_REP, 
						group="retained",
						n=sum(ol), 
						percent=100* sum(ol) / length(ol),
						p_val=p
						)
					row2 <- data.frame(
						mark=MARK, 
						cell_type=CELL_TYPE, 
						RNAseq_replicate=RNA_REP, 
						mark_replicate=MARK_REP, 
						group="all",
						n=sum(olAll), 
						percent=100* sum(olAll) /length(olAll),
						p_val=NA
						)
	
					allDF <- rbind(allDF, row1)
					allDF <- rbind(allDF, row2)					





					#---------------------------------------------------
					# old stuff
					#---------------------------------------------------
					

					
				} 
			}
		}
	}
}


# reformat as factor with order
allDF$cell_type = factor(allDF$cell_type, c("e125", "e145", "ISC", "Paneth",  "enterocyte", "ae"))
allDF$mark = factor(allDF$mark, c("k4me3", "k27ac"))



# Run the functions length, mean, and sd on the value of "change" for each group, 
# broken down by sex + condition
cdata <- ddply(allDF, c("mark", "cell_type", "mark_replicate", "group"), summarise,
               N    = length(n),
               mean = mean(n),
               sd   = sd(n),
               se   = sd / sqrt(N),
               percent_mean = mean(percent),
               percent_sd   = sd(percent),
               percent_se   = percent_sd / sqrt(N),
               mean_p = mean(p_val)
)

 
g <- ggplot(cdata, aes(x=group, y=percent_mean, fill=cell_type)) + 
	geom_bar(stat="identity", color="black") + 
	geom_errorbar(aes(ymax = percent_mean + percent_sd, ymin=percent_mean - percent_sd), width=0.25) + 
	geom_text(aes(label=signif(percent_mean,3)), vjust=-1.25) + 
  geom_text(aes(x=1.5, label=ifelse(is.na(mean_p), "", paste0("p=", signif(mean_p,3)))), vjust=-2.5) + 
  #~ 	geom_text(aes(label=paste("n=",round(mean,1)), y=0), angle=90) + 
	ylim(0,1.4*max(cdata$percent_mean, na.rm=TRUE)) +
	facet_grid(mark*mark_replicate~cell_type) + 
	theme_bw() + scale_fill_manual(values=COL_CELL_TYPE) + 
	xlab("Introns") + ylab("Introns overlapping chromatin mark [%]")
	
ggsave(paste0(outPrefix, ".introns_overlapping_marks.barplot_percent.pdf"))



g <- ggplot(cdata, aes(x=group, y=mean, fill=cell_type)) + 
	geom_bar(stat="identity", color="black") + 
	geom_errorbar(aes(ymax = mean + sd, ymin=mean - sd), width=0.25) + 
	geom_text(aes(label=round(mean,1)), vjust=-1.25) + 
	ylim(0,1.2*max(cdata$mean, na.rm=TRUE)) +
	facet_grid(mark*mark_replicate~cell_type) + 
	theme_bw() + scale_fill_manual(values=COL_CELL_TYPE) + 
	xlab("Introns") + ylab("Introns overlapping chromatin mark")
	
ggsave(paste0(outPrefix, ".introns_overlapping_marks.barplot.pdf"))




