#!/bin/bash

# Gets full path of analysis folders where trimmed fasta files and name files are found. Descend into each directory and run a mothur batch script in each.

# Could be made interactive to read in a top directory from which to search. For now, must be started from a terminal open to the correct folder.

regions=$(echo v4)

for region in $regions
do

prefix=combined_seqs
prefixReset=$prefix
reference=/data/results/home/mmeier/dbs/mothur/silva/silva.nr_v128.$region.align
refname=${reference##*/}
trainset=/data/results/home/mmeier/dbs/mothur/greengenes/gg_13_8_99.fasta
taxonomy=/data/results/home/mmeier/dbs/mothur/greengenes/gg_13_8_99.gg.tax
#trainset=/data/results/home/mmeier/dbs/mothur/rdp/trainset16_022016.rdp/trainset16_022016.rdp.fasta
#taxonomy=/data/results/home/mmeier/dbs/mothur/rdp/trainset16_022016.rdp/trainset16_022016.rdp.tax
rootdir=/data/results/home/mmeier/data/2017Consortia_project/Renuka_region_subsampled/
#rootdir=/data/results/home/mmeier/data/2017Consortia_project/Renuka_redo/Auto_user_SN2-5-AgNP_V6_VSL_soil_48_025


for directory in $(find $rootdir -wholename "*/analysis" -type d )
#for directory in $(echo /data/results/home/mmeier/data/2017Consortia_project/Renuka_redo/Auto_user_SN2-5-AgNP_V6_VSL_soil_48_025/analysis)
	 	
	do echo "Starting analysis for $directory ..."
	cd $directory
	echo "Current working directory:"
	echo $PWD


	prefix=$prefixReset
	
	doneChecker=$(ls | grep $region | grep biom)
	if [ -f "$doneChecker" ]
	then echo "Biom file already present. Skipping directory ${directory}."
	
	else

	fileChecker=$(ls | grep .good.unique.align$ | grep $region)

	if  [ -f "$fileChecker" ]

		then echo "Alignment previously complete. Skipping step to save computational resources."

		else echo "Alignment not previously done; continuing with preprocessing and alignment."
		
		echo "Starting summary.seqs in $directory"
		mothur "#summary.seqs(fasta=$prefix.fasta, processors=15)"
		logfile=$(ls -t | grep log | head -n 1)
		min=$(cat $logfile  | grep 2.5%-tile | awk '{print $4}')
		max=$(cat $logfile  | grep 97.5%-tile | awk '{print $4}')
		echo $logfile
		echo $min
		echo $max
		echo "Starting first round of screen.seqs in $directory"
		mothur "#screen.seqs(fasta=$prefix.fasta, group=stability.contigs.groups, maxambig=0, minlength=$min, maxlength=$max, summary=combined_seqs.summary, maxhomop=6, processors=15)"

		# mothur "#screen.seqs(fasta=$prefix.fasta, group=stability.contigs.groups, maxambig=0, optimize=start-end-minlength, criteria=50, summary=combined_seqs.summary, maxhomop=6, processors=12)"
		echo "Getting unique sequencess in $directory"
		mothur "#unique.seqs(fasta=$prefix.good.fasta)"
		echo "Counting sequences in $directory"
	        mothur "#count.seqs(name=$prefix.good.names, group=stability.contigs.good.groups)"

		### This block could be repeated with different reference files.
		
			echo "Running alignment in $directory"
		        mothur "#align.seqs(fasta=$prefix.good.unique.fasta, reference=$reference, flip=t, processors=15)"
	
			# Rename output files according to the reference to which they were aligned!
			mv $prefix.good.unique.align $prefix.$refname.good.unique.align
			mv $prefix.good.unique.align.report $prefix.$refname.good.unique.align.report
			cp $prefix.good.count_table $prefix.$refname.good.count_table
	
	fi

	# reset prefix to $prefix.$refname for remaining steps
	
		prefixReset=$prefix
		prefix=$prefix.$refname

	
	if [ -f "$fileChecker" ]
		then echo "Also skipping summary sequences since we assume it is done after alignment."
		else echo "Getting summary of sequences for $directory"

		mothur "#summary.seqs(fasta=$prefix.good.unique.align, count=$prefix.good.count_table, processors=15)"

	fi

	logfile=$(ls -t | grep log | head -n 1)
	start=$(cat $logfile | grep 2.5%-tile | awk '{print $2}')
	end=$(cat $logfile  | grep 75%-tile | awk '{print $3}')

	# Alternate code for screen.seqs:
	# mothur "#screen.seqs(fasta=$prefix.good.unique.align, count=$prefix.good.count_table, start=$start, end=$end, processors=12)"
	
	#fileChecker=$(ls | grep .good.unique.good.align$ | grep $region)
	#if [ -f "$fileChecker" ]
                #then echo "Skipping screen sequences"
		#else echo "Running screen.seqs on ${directory}"
	mothur "#screen.seqs(fasta=$prefix.good.unique.align, count=$prefix.good.count_table, optimize=start, criteria=75, minlength=150, maxlength=200, processors=15)"
	#fi

	echo $directory
       	mothur "#filter.seqs(fasta=$prefix.good.unique.good.align, vertical=T, processors=14)"
	echo $directory
        mothur "#unique.seqs(fasta=$prefix.good.unique.good.filter.fasta, count=$prefix.good.good.count_table)"
	echo $directory
        mothur "#pre.cluster(fasta=$prefix.good.unique.good.filter.unique.fasta, count=$prefix.good.unique.good.filter.count_table, diffs=2, processors=14)"
	echo $directory
        mothur "#chimera.vsearch(fasta=$prefix.good.unique.good.filter.unique.precluster.fasta, count=$prefix.good.unique.good.filter.unique.precluster.count_table, dereplicate=t, processors=14)"
	echo $directory
	mothur "#remove.seqs(fasta=$prefix.good.unique.good.filter.unique.precluster.fasta, accnos=$prefix.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)"
	echo $directory
        mothur "#classify.seqs(fasta=$prefix.good.unique.good.filter.unique.precluster.pick.fasta, count=$prefix.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference=$trainset, taxonomy=$taxonomy, cutoff=80, processors=14)"
	echo $directory
        mothur "#remove.lineage(fasta=$prefix.good.unique.good.filter.unique.precluster.pick.fasta, count=$prefix.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy=$prefix.good.unique.good.filter.unique.precluster.pick.gg.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)"
	echo $directory
	mothur "#dist.seqs(fasta=$prefix.good.unique.good.filter.unique.precluster.pick.pick.fasta, cutoff=0.03, processors=14)"
	echo $directory
        mothur "#cluster(column=$prefix.good.unique.good.filter.unique.precluster.pick.pick.dist, count=$prefix.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, processors=14)"
	echo $directory
        mothur "#make.shared(list=$prefix.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=$prefix.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03)"
	echo $directory
        mothur "#classify.otu(list=$prefix.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=$prefix.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=$prefix.good.unique.good.filter.unique.precluster.pick.gg.wang.pick.taxonomy, label=0.03)"
	echo $directory
        mothur "#phylotype(taxonomy=$prefix.good.unique.good.filter.unique.precluster.pick.gg.wang.pick.taxonomy)"
	echo $directory
        mothur "#make.shared(list=$prefix.good.unique.good.filter.unique.precluster.pick.gg.wang.pick.tx.list, count=$prefix.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=1)"
	echo $directory
        mothur "#classify.otu(list=$prefix.good.unique.good.filter.unique.precluster.pick.gg.wang.pick.tx.list, count=$prefix.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=$prefix.good.unique.good.filter.unique.precluster.pick.gg.wang.pick.taxonomy, label=1)"
	echo $directory
	echo "The following shared file was used:"
	currentshared=$(ls -t | grep shared$ | head -n 1)
	echo $currentshared
        mothur "#count.groups(shared=$currentshared)"
	echo $directory
        mothur "#rarefaction.single(shared=$currentshared, calc=sobs, freq=100)"

	##### TREES AND BIOM FILES
	echo $directory
	mothur "#dist.seqs(fasta=$prefix.good.unique.good.filter.unique.precluster.pick.pick.fasta, output=lt, processors=14)"
	echo $directory
	mothur "#cluster(phylip=$prefix.good.unique.good.filter.unique.precluster.pick.pick.phylip.dist, count=$prefix.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, processors=14)"
	
	### Bin sequences here - creates FASTA file with OTU name for each entry in the fasta file, as well as group information

	mothur "#bin.seqs(list=$prefix.good.unique.good.filter.unique.precluster.pick.pick.phylip.opti_mcc.list, fasta=$prefix.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=$prefix.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03)"

	### OUTPUT IS FASTA

	# Calculate distance
#	echo $directory
#	mothur "#dist.seqs(fasta=$prefix.good.unique.good.filter.unique.precluster.pick.pick.phylip.opti_mcc.0.03.fasta, output=lt, processors=12)"
	# Cluster and make tree
#	echo $directory
#	mothur "#cluster(phylip=$prefix.good.unique.good.filter.unique.precluster.pick.pick.phylip.opti_mcc.0.03.phylip.dist, count=$prefix.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, processors=12)"


	echo $directory
	mothur "#clearcut(phylip=$prefix.good.unique.good.filter.unique.precluster.pick.pick.phylip.dist)"

	echo $directory
	## Last FASTA with info: $prefix.align.good.unique.good.filter.unique.precluster.pick.pick.phylip.opti_mcc.0.03.fasta
	
# Get group names of those that are sampled well

goodSamples=()
total=$(wc -l *map | grep total | awk '{print $1}')
for i in *map
do
lines=$(wc -l $i | awk '{print $1}')
name=$(wc -l $i | awk '{print $2}')
pct=$(echo "100*$lines/$total" | bc)
echo $name
echo $lines
echo $total
echo $pct

if [ "$pct" -lt 1 ]
then
echo "$name is NOT part of the dataset"
else
group=$(echo $name | rev | cut -d. -f 2  | rev)
goodSamples+=($group)
fi
done

goodGroups=$(echo "${goodSamples[@]}" | sed 's/\s/-/g')

echo ${goodSamples[@]}
echo $goodGroups

	# Make new shared file
	echo $directory
	mothur "#make.shared(list=$prefix.good.unique.good.filter.unique.precluster.pick.pick.phylip.opti_mcc.list, count=$prefix.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)"
	echo $directory
	mothur "#normalize.shared(shared=$prefix.good.unique.good.filter.unique.precluster.pick.pick.phylip.opti_mcc.shared, groups=$goodGroups)"
	# Outputs: combined_seqs.silva.nr_v128.v6.align.good.unique.good.filter.unique.precluster.pick.pick.phylip.opti_mcc.0.03.norm.shared

	# Make biom format file
	echo $directory
	mothur "#make.biom(shared=$prefix.good.unique.good.filter.unique.precluster.pick.pick.phylip.opti_mcc.0.03.norm.shared, constaxonomy=$prefix.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy, label=0.03)"
	echo "Finished with $directory, proceeding..."
	sleep 3

fi # Finish if Biom file present, go to next directory.

done
done
