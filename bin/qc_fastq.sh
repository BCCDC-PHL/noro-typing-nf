FORWARD=$1
REVERSE=$2
SAMPLE=$(echo $FORWARD | cut -d_ -f1)
OUTFILE="${SAMPLE}.QC_fastq.tsv"

if [[ ${FORWARD} == *.gz ]]; then
	gunzip -c ${FORWARD} | fastx_quality_stats -o ${SAMPLE}.R1.raw.tsv -Q 33 &&
	gunzip -c ${REVERSE} | fastx_quality_stats -o ${SAMPLE}.R2.raw.tsv -Q 33

elif [[ ${FORWARD} == *fastq ]]; then
	fastx_quality_stats -i ${FORWARD} -o ${SAMPLE}.R1.raw.tsv -Q 33 &&
	fastx_quality_stats -i ${REVERSE} -o ${SAMPLE}.R2.raw.tsv -Q 33

fi

# write header line
echo -e -n 'sample\ttotal_reads_r1\tmax_read_len_r1\t' > $OUTFILE
echo -e -n 'num_max_len_reads_r1\tmean_base_quality_r1\t' >> $OUTFILE
echo -e -n 'total_reads_r2\tmax_read_len_r2\t' >> $OUTFILE
echo -e 'num_max_len_reads_r2\tmean_base_quality_r2' >> $OUTFILE

echo -e -n ${SAMPLE} >> $OUTFILE

for i in {1,2}; do 
	BASES=`awk '{s+=$2}END{print s}' ${SAMPLE}.R$i.raw.tsv`
	PERQUALITY=`awk '{s+=$5}END{print s}' ${SAMPLE}.R$i.raw.tsv`
	AVGQUALITY=`echo "scale=2; $PERQUALITY/$BASES" | bc`
	READTOTAL=`sed -n -e 2p ${SAMPLE}.R$i.raw.tsv | awk '{ print $2 }'`
	MAXLENGTH=`sed '$!d' ${SAMPLE}.R$i.raw.tsv | awk '{ print $1 }'`
	NOMAXREADS=`sed '$!d' ${SAMPLE}.R$i.raw.tsv | awk '{ print $2 }'`

	echo -e -n '\t'$READTOTAL'\t'$MAXLENGTH'\t'$NOMAXREADS'\t'$AVGQUALITY >> $OUTFILE
done 

echo "" >> $OUTFILE