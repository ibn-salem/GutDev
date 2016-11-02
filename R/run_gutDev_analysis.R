########################################################################
# A script to analyse intron-retention with chromatin marks in gut development
########################################################################

require(gdata)
require(GenomicRanges)
require(rtracklayer)
require(ggplot2)
require(BSgenome.Mmusculus.UCSC.mm9) # for proper seqinfo object
require(RColorBrewer)
require(plyr)

COL_CELL_TYPE=brewer.pal(8, "Pastel1")
COL_ALL=brewer.pal(8, "Dark2")[8]

CELL_TYPE_LEVELS <- c("e125", "ISC", "Paneth", "enterocyte")

VERSION="v02"

# design / dimensions

# mark
#	cell_type
#		 time
#			replicate
# 
#/data/Chromatin/k4me3_peaks/K4me3_called_peaks_e125_1.csv

chromMarksDesign <- "data/chromatin_marks_design.csv"
intronsDesign <- "data/retained_intron_design.csv"

allIntronFiles <- "data/Intron_Retention_bed_files/Intron_Retention_Annotated/intron_retained_in_all_genome_direction_annotated_miso_ordered.bed.mm10.BED.mm9"

# read metadata of marks
chromMarkFiles <- read.table(chromMarksDesign, header=TRUE)
# take only a subest with existing RNA-seq data
chromMarkFiles <- subset(chromMarkFiles, cell_type %in% CELL_TYPE_LEVELS)
chromMarkFiles$cell_type <- factor(chromMarkFiles$cell_type, levels=CELL_TYPE_LEVELS)

# use only replicate 1
#chromMarkFiles <- subset(chromMarkFiles, replicate==1)

# read metadata of introns
intronFiles <- read.table(intronsDesign, header=TRUE)
intronFiles <- intronFiles[!is.na(intronFiles$file),]
intronFiles$cell_type <- factor(intronFiles$cell_type, levels=CELL_TYPE_LEVELS)
# remove e145 

outPrefix <- paste0("results/", VERSION, "_GutDev")
dir.create(dirname(outPrefix), recursive=TRUE, showWarnings = FALSE)
# use only replicate 1
#intronFiles <- subset(intronFiles, intron_replicate==1)

# get seqinfo object from annotation package
genome <- BSgenome.Mmusculus.UCSC.mm9
seqInfo <- seqinfo(genome)

# parse chromatin mark peaks as GenomicRanges objects
markGRL <- lapply(chromMarkFiles$file, function(f){
	if(is.na(f)){
		GRanges()
	}else{
		gr <- GRanges(read.delim(as.character(f), header=TRUE), seqinfo=seqInfo)
	}
})


# parse all introns
allIntrons <- import(allIntronFiles, format="bed", seqinfo=seqInfo)

# parse all retained introns
# parse introns as GenomicRanges objects
intronGRL <- lapply(intronFiles$file, function(f){
	if(is.na(f)){
		GRanges()
	}else{
		gr <- import(as.character(f), format="bed", seqinfo=seqInfo)
	}
})
#~ intronGRL <- c(list("all"=allIntrons), intronGRL)
# intronGRL <- GRangesList(intronGRL)


#-----------------------------------------------------------------------
# Compare length of introns
#-----------------------------------------------------------------------

# create a data.frame with all intron sizes from retained introns
widthDF <- do.call("rbind", lapply(1:nrow(intronFiles), function(i){
  data.frame(intronFiles[rep(i, length(intronGRL[[i]])),], width=width(intronGRL[[i]]))
}))

# add all introns
allIntronsWidthDF <- data.frame(
    name="All Introns", 
    cell_type="all introns",
    intron_replicate=1,
    file=allIntronFiles,
    width=width(allIntrons)
  )
  
widthDF <- rbind(
    allIntronsWidthDF,
    widthDF
)


g <- ggplot(widthDF, aes(x=as.factor(intron_replicate), y=width, fill=cell_type)) + 
    geom_boxplot()+
    facet_grid(.~cell_type, scales="free_x") +
    theme_bw() + scale_fill_manual(values=c(COL_ALL, COL_CELL_TYPE)) + 
    xlab("Retained introns replicates") + ylab("Intron size [bp]")
ggsave(paste0(outPrefix, ".intron_sizes.boxplot.pdf"))

g +  scale_y_log10() 
ggsave(paste0(outPrefix, ".intron_sizes.boxplot_log10.pdf"))

#-----------------------------------------------------------------------
# Plot number of retained introns per cell and across replicates
#-----------------------------------------------------------------------
intronFiles$n <- sapply(intronGRL, length)

g <- ggplot(intronFiles, aes(x=intron_replicate, y=n, fill=cell_type)) + 
  geom_bar(stat="identity", color="black") + 
  geom_text(aes(label=n), vjust=-1.25) + 
  facet_grid(~cell_type) + 
  theme_bw() + scale_fill_manual(values=COL_CELL_TYPE) + 
  xlab("RNA-seq replicates") + ylab("Retained Introns")

ggsave(paste0(outPrefix, ".retained_introns.barplot.pdf"))

#-----------------------------------------------------------------------
# Plot number and size distribution of chromatin mark peaks
#-----------------------------------------------------------------------

# create a data.frame with all chrmatin peaks
peakDF <- do.call("rbind", lapply(1:nrow(chromMarkFiles), function(i){
  data.frame(chromMarkFiles[rep(i, length(markGRL[[i]])),], width=width(markGRL[[i]]))
}))

# reformat as factor with order
peakDF$cell_type = factor(peakDF$cell_type, c("e125", "e145", "ISC", "Paneth",  "enterocyte", "ae"))
peakDF$mark = factor(peakDF$mark, c("k4me3", "k27ac"))

# count peaks per condition
countPeakDF <- ddply(peakDF, c("mark", "cell_type", "replicate"), summarise, 
                     count=length(width))

# plot number of peaks
g <- ggplot(countPeakDF, aes(x=replicate, y=count, fill=cell_type)) + 
  geom_bar(stat="identity", color="black") + 
  geom_text(aes(label=count), hjust=-1.25, angle=90) + 
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

# iterate over chrom mark types
for (MARK in as.character(unique(chromMarkFiles$mark))){
	
	# iteratre over cell types
	for (CELL_TYPE in as.character(unique(intronFiles$cell_type))){
		
		# iterate over intron retation replicates
		for (INTRON_REP in unique(intronFiles$intron_replicate)) {

			i <- which(intronFiles$cell_type == CELL_TYPE & intronFiles$intron_replicate == INTRON_REP)
			
			# iterate chrom-mark ChIP-seq replicates
			for (MARK_REP in unique(chromMarkFiles$replicate)) {

				m <- which(chromMarkFiles$mark == MARK 
					& chromMarkFiles$cell_type == CELL_TYPE
					& chromMarkFiles$replicate == MARK_REP)
				
				# debug message
				message(paste(MARK, CELL_TYPE, INTRON_REP, MARK_REP, i, m))
				
				ol <- countOverlaps(intronGRL[[i]], markGRL[[m]], ignore.strand=TRUE) > 0
				olAll <- countOverlaps(allIntrons, markGRL[[m]], ignore.strand=TRUE) > 0
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
					intron_replicate=INTRON_REP, 
					mark_replicate=MARK_REP, 
					group="retained",
					n=sum(ol), 
					percent=100* sum(ol) / length(ol),
					p_val=p
					)
				row2 <- data.frame(
					mark=MARK, 
					cell_type=CELL_TYPE, 
					intron_replicate=INTRON_REP, 
					mark_replicate=MARK_REP, 
					group="all",
					n=sum(olAll), 
					percent=100* sum(olAll) /length(olAll),
					p_val=NA
					)

				allDF <- rbind(allDF, row1)
				allDF <- rbind(allDF, row2)
				
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




