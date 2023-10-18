name=$1
rn1='-'
id=${name#*$rn1}
id=${id%'/'}
echo $id

for L in L001 L002 L003
do 
for R in R1 R2 I1 
do
#mv $id*$L*$R* ${id}_${L}_${I}_001.fastq.gz
echo $id*$L*$R* ${id}_${L}_${I}_001.fastq.gz
done
done

