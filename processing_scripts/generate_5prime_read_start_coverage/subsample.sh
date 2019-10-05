all=("$@")
bamfile=${all[0]}
fraction=${all[1]}
#
bamfilename=${bamfile%.bam}
#samtools view -s $fraction -b $bamfile | samtools sort - $bamfilename.$fraction.sorted
samtools view -h $bamfile "$fraction" | samtools sort - $bamfilename.$fraction.sorted
#samtools reheader $bamfilename.$fraction.bam
samtools index $bamfilename.$fraction.sorted.bam
