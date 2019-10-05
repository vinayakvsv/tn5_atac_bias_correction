all=("$@")
indir=${all[0]}
#
pillai="/n/data2/mgh/ragon/pillai"
strsplit="$pillai/pipelines/epigenetics_pipeline/motif_footprints/atacseq_footprinting/transpos_site_cov/split_bam_strand.sh"
#
bams=$(find $indir -name *.bam)
#
echo $bams
echo $strsplit
#
for i in $bams; do
	name=$(basename $i)
	outname=${name%.bam}
	echo $outname
	echo $i
	bsub -q mcore -n 6 -W 12:00 -o out.$outname.txt -e err.$outname.txt bash $strsplit $i
done

