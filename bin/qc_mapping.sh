FILENAME=$1
SAMPLE=$(echo $FILENAME | cut -d. -f1)
TMPFILE=${SAMPLE}.raw.tsv
OUTFILE=${SAMPLE}.coverage.tsv

echo -e 'sample\tpercent_coverage\tmean_coverage' > $OUTFILE

bedtools genomecov -ibam $FILENAME -g REFERENCE > $TMPFILE

PERCOV=`head -1 $TMPFILE | awk '{ print $5}'`
PERCOV=`echo $PERCOV | sed -e 's/[eE]+*/\\*10\\^/'`
PERCOV=`echo "scale=10; 1 - $PERCOV" | bc`
PERCOV=`echo "$PERCOV * 100" | bc`

COVCOL=(`awk '{ print $2 }' $TMPFILE`)
FREQCOL=(`awk '{ print $5 }' $TMPFILE`)
MULTIPLICATION=0
LENGTH=${#COVCOL[@]}
LENGTH=`expr $LENGTH / 2 - 1`

for i in `seq 1 $LENGTH`; do
	VALUE1=`echo ${FREQCOL[$i]} | sed -e 's/[eE]+*/\\*10\\^/'`
	VALUE2=`echo "scale=10; ${COVCOL[$i]} * ($VALUE1)" | bc | awk '{printf "%.02f", $0}'`
	MULTIPLICATION=`echo "scale=10; $MULTIPLICATION + $VALUE2" | bc | awk '{printf "%.02f", $0}'`
done

echo -e $SAMPLE'\t'$PERCOV'\t'$MULTIPLICATION >> $OUTFILE