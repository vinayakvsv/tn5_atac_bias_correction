module load seq/samtools/1.2
#
all=("$@")
#
bamfile=${all[0]}
echo $bamfile
#
bamfilename=$(basename $bamfile)
bamname=${bamfilename%.bam}
echo $bamname
#positive strand
samtools view -h -F 16 $bamfile | samtools sort - $bamname.pos.sorted
samtools index $bamname.pos.sorted.bam
#negative strand
samtools view -h -f 16 $bamfile | samtools sort - $bamname.neg.sorted
samtools index $bamname.neg.sorted.bam
