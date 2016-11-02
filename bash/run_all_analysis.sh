

#=======================================================================
# Prepare data
#=======================================================================

#-----------------------------------------------------------------------
# download liftOver tool from USCS and chain files to convert mm10 cooridnates to mm9

#~ # UCSC liftover chains
#~ mkdir -p data/UCSC
#~ wget -P data/UCSC http://hgdownload.cse.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm9.over.chain.gz
#~ gunzip data/UCSC/*.gz

#~ # download liftOver tool from UCSC:
#~ mkdir bin
#~ wget -P bin http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
#~ chmod u+x bin/liftOver

#-----------------------------------------------------------------------
# convert introns from mm10 to mm9

for F in data/Intron_Retention_bed_files/*/*.bed ; do
	
	echo "File -> $F"
	
		#~ | awk -F"\t" '{if ($4 == "+") {print $1,$2,$3,".","0",$4} else {print $1,$3,$2,".","0",$4} }' OFS="\t" \

	# convert into standard .bed file
	cat ${F} \
		| awk -F"\t" '{if ($3 >= $2) {print $1,$2,$3} else {print $1,$3,$2} }' OFS="\t" \
		> ${F}.mm10.BED
	
	# liftover to mm9
	./bin/liftOver \
		${F}.mm10.BED \
		data/UCSC/mm10ToMm9.over.chain \
		${F}.mm10.BED.mm9 \
		${F}.mm10.BED.mm9_unmapped
done



#-----------------------------------------------------------------------
# Run R script to convert .xlsx tables to standart bed files
#~ Rscript l


#=======================================================================
# Analyse Overlap
#=======================================================================
