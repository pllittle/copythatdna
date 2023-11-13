#!/bin/sh

[ ! -z $srcPL_dnaGetCNT ] && [ $srcPL_dnaGetCNT -eq 1 ] && return 0
[ -z "$git_dir" ] && git_dir=$(cd $(dirname $BASH_SOURCE)/../..; pwd)

for repo in baSHic; do
	repo_dir=$git_dir/$repo
	check_array $repo baSHic && tmp_url=https://github.com/pllittle/$repo.git
	
	while true; do
		if [ ! -d "$repo_dir" ]; then
			cd "$git_dir"
			git clone "$tmp_url" >&2
			[ $? -eq 0 ] && break
		else
			cd "$repo_dir"
			git pull >&2
			[ $? -eq 0 ] && break
		fi
		echo -e "Some error in cloning $repo, contact pllittle" >&2 && return 1
	done
	
done

for fn in base colors getEnv install linux_git linux_perl; do
	. $git_dir/baSHic/scripts/$fn.sh
	[ ! $? -eq 0 ] && echo -e "Some error in src-ing baSHic's $fn script" >&2 \
		&& return 1
done

unset repo_dir fn

# DNA Copy Number: Get bin read counts
make_contig_sizes(){
	local sam_dir bam_fn output_fn
	local rm_contigs cmd status
	
	while [ ! -z $1 ]; do
		case $1 in
			-s | --sam_dir )
				shift
				sam_dir="$1"
				;;
			-b | --bam_fn )
				shift
				bam_fn="$1"
				;;
			-o | --output_fn )
				shift
				output_fn="$1"
				;;
			-r | --rm_contigs )
				shift
				rm_contigs="$1"
				;;
		esac
		shift
	done
	
	[ -z $sam_dir ] 			&& echo "Add -s <samtools dir>" >&2 && return 1
	[ -z $bam_fn ] 				&& echo "Add -b <bam filename>" >&2 && return 1
	[ -z $output_fn ] 		&& echo "Add -o <output filename>" >&2 && return 1
	[ -z "$rm_contigs" ] 	&& echo "Not removing any contigs" >&2
	
	cmd="$sam_dir/bin/samtools view -H $bam_fn | grep '^@SQ'"
	cmd="$cmd | cut -f2-3 | sed 's|SN:||g' | sed 's|LN:||g'"
	[ ! -z "$rm_contigs" ] && cmd="$cmd	| grep -v -E '$rm_contigs'"
	cmd="$cmd > $output_fn"
	eval $cmd
	
}
make_genome_windows(){
	local bt_dir input_fn output_fn windsize
	local status
	
	while [ ! -z $1 ]; do
		case $1 in
			-b | --bt_dir )
				shift
				bt_dir="$1"
				;;
			-i | --input_fn )
				shift
				input_fn="$1"
				;;
			-o | --output_fn )
				shift
				output_fn="$1"
				;;
			-w | --windsize )
				shift
				windsize="$1"
				;;
		esac
		shift
	done
	
	[ -z $bt_dir ] 		&& echo "Add -b <bedtools dir>" >&2 && return 1
	[ -z $input_fn ] 	&& echo "Add -i <input contig sizes>" >&2 && return 1
	[ -z $output_fn ] && echo "Add -o <output bed file>" >&2 && return 1
	[ -z $windsize ] 	&& echo "Add -w <window size>" >&2 && return 1
	
	[ -f $output_fn ] && return 0
	
	echo -e "`date`: Start makewindows" >&2
	$bt_dir/bin/bedtools makewindows \
		-g $input_fn -w $windsize \
		> $output_fn
	status=$?
	[ ! $status -eq 0 ] && new_rm $output_fn \
		&& echo -e "`date`: Error makewindows" >&2 \
		&& return 1
	echo -e "`date`: Finish makewindows" >&2
	
	return $status
}
get_GCcontent(){
	local bt_dir fasta_fn input_fn output_fn
	local orig_ncols status
	
	while [ ! -z $1 ]; do
		case $1 in
			-b | --bt_dir )
				shift
				bt_dir="$1"
				;;
			-f | --fasta_fn )
				shift
				fasta_fn="$1"
				;;
			-i | --input_fn )
				shift
				input_fn="$1"
				;;
			-o | --output_fn )
				shift
				output_fn="$1"
				;;
		esac
		shift
	done
	
	[ -z $bt_dir ] 		&& echo "Add -b <bedtools dir>" >&2 && return 1
	[ -z $fasta_fn ] 	&& echo "Add -f <reference fasta>" >&2 && return 1
	[ -z $input_fn ] 	&& echo "Add -i <input bed file>" >&2 && return 1
	[ -z $output_fn ] && echo "Add -o <output bed file>" >&2 && return 1
	
	[ -s $output_fn ] && return 0
	orig_ncols=$(head -n 1 $input_fn | sed 's|\t|\n|g' | wc -l)
	
	echo -e "`date`: GC content started" >&2
	$bt_dir/bin/bedtools nuc -fi $fasta_fn \
		-bed $input_fn | tail -n +2 \
		| cut -f 1-${orig_ncols},$((orig_ncols+2)) \
		> $output_fn
	status=$?
	[ ! $status -eq 0 ] && new_rm $output_fn \
		&& echo -e "`date`: Error with GC content" >&2 \
		&& return 1
	echo -e "`date`: GC content finished" >&2
	
	return $status
}
get_windowCounts(){
	local bam_fn fasta_fn bt_dir sam_dir status min_mapq
	local rm_contigs out_dir windsize clean_dir tmp_fn ncores
	local bed_fn rm_bam
	
	clean_dir=0; ncores=1; rm_bam=0
	while [ ! -z $1 ]; do
		case $1 in
			-b | --bam_fn )
				shift
				bam_fn="$1"
				;;
			-e | --bed_fn )
				shift
				bed_fn="$1"
				;;
			-f | --fasta_fn )
				shift
				fasta_fn="$1"
				;;
			-i | --bt_dir )
				shift
				bt_dir="$1"
				;;
			-j | --sam_dir )
				shift
				sam_dir="$1"
				;;
			-k | --clean_dir )
				clean_dir=1
				;;
			-n | --name )
				shift
				name="$1"
				;;
			-o | --out_dir )
				shift
				out_dir="$1"
				;;
			-q | --min_mapq )
				shift
				min_mapq="$1"
				;;
			-r | --rm_contigs )
				shift
				rm_contigs="$1"
				;;
			-t | --ncores )
				shift
				ncores="$1"
				;;
			-w | --windsize )
				shift
				windsize="$1"
				;;
			-m | --rm_bam )
				rm_bam=1
				;;
		esac
		shift
	done
	
	# Check inputs
	[ -z $bam_fn ] 				&& echo "Add -b <bam fn>" >&2 && return 1
	[ -z $fasta_fn ] 			&& echo "Add -f <fasta fn>" >&2 && return 1
	[ -z $bt_dir ] 				&& echo "Add -i <bedtools install dir>" >&2 && return 1
	[ -z $sam_dir ] 			&& echo "Add -j <samtools install dir>" >&2 && return 1
	[ -z $name ] 					&& echo "Add -n <extra name tag>" >&2 && return 1
	[ -z $out_dir ] 			&& echo "Add -o <out dir>" >&2 && return 1
	[ -z "$rm_contigs" ] 	&& echo "Add -r <remove contigs>" >&2 && return 1
	[ -z $windsize ] 			&& echo "Add -w <window size>" >&2 && return 1
	[ -z $min_mapq ]			&& echo "Add -q <minimum mapping quality>" >&2 && return 1
	[ $windsize -eq 0 ] && [ -z $bed_fn ] && echo "Add -e <bed fn>" >&2 && return 1
	
	if [ ! -f $bt_dir/bin/bedtools ] \
		|| [ ! $($bt_dir/bin/bedtools --version &> /dev/null; echo $?) -eq 0 ]; then
		echo "Install bedtools or setup environment!" >&2 && return 1
	fi
	if [ ! -f $sam_dir/bin/samtools ] \
		|| [ ! $($sam_dir/bin/samtools --version &> /dev/null; echo $?) -eq 0 ]; then
		echo "Install samtools or setup environment!" >&2 && return 1
	fi
	
	new_mkdir $out_dir
	
	# Clear directory
	[ $clean_dir -eq 1 ] && [ `ls $out_dir | wc -l` -gt 0 ] \
		&& new_rm $out_dir/windows.$name.$windsize.* \
			$out_dir/$name.*
	
	# If final file exists, done
	tmp_fn=$out_dir/windows.$name.$windsize.gc.mMAPQ${min_mapq}.bed.gz
	if [ -f $tmp_fn ]; then
		[ $rm_bam -eq 1 ] && new_rm $out_dir/$name.sort.bam \
			$out_dir/$name.sort.bam.bai
		echo -e "`date`: Job completed" >&2 && return 0
	fi
	
	# If sort.bam files exist, clean up
	new_rm $out_dir/$name.bam $out_dir/$name.sort.bam.tmp.*.bam \
		$out_dir/$name.sort.bam.tmp.*.bam.bai \
		$out_dir/windows.$name.$windsize.bed \
		$out_dir/windows.$name.$windsize.gc.bed \
		$out_dir/windows.$name.$windsize.gc.mMAPQ${min_mapq}.bed
	
	# Make genome file
	make_contig_sizes -s $sam_dir -b $bam_fn \
		-o $out_dir/$name.chrom.sizes -r "$rm_contigs"
	
	# Subset chromosomes from bam
	tmp_fn=$out_dir/$name.sort.bam
	if [ ! -f $tmp_fn ]; then
		export OMP_NUM_THREADS=$ncores
		
		# Get chr1-22 X Y MT basically ...
		echo -e "`date`: Subset contigs w/ samtools" >&2
		$sam_dir/bin/samtools view -h -o $out_dir/$name.bam \
			$bam_fn `cut -f1 $out_dir/$name.chrom.sizes` >&2
		status=$?
		if [ ! $status -eq 0 ]; then
			new_rm $out_dir/$name.bam
			echo -e "`date`: Error with subset contigs" >&2
			return 1
		fi
		
		# Sort subsetted bam file
		echo -e "`date`: Sort subsetted bam w/ samtools" >&2
		$sam_dir/bin/samtools sort -@ $((ncores-1)) \
			-o $out_dir/$name.sort.bam $out_dir/$name.bam >&2
		status=$?
		if [ ! $status -eq 0 ]; then
			new_rm $out_dir/$name.sort.bam
			echo -e "`date`: Error with sorting bam" >&2
			return 1
		fi
		new_rm $out_dir/$name.bam
		
		# Create bam index file
		echo -e "`date`: Create index file w/ samtools" >&2
		$sam_dir/bin/samtools index -b \
			-@ $((ncores-1)) $out_dir/$name.sort.bam >&2
		[ ! $? -eq 0 ] && echo "Error in index" >&2 && return 1
	fi
	
	# Make windows bed
	if [ -z "$bed_fn" ] && [ $windsize -gt 0 ]; then
		export OMP_NUM_THREADS=1
		make_genome_windows -b $bt_dir \
			-i $out_dir/$name.chrom.sizes \
			-o $out_dir/windows.$name.$windsize.bed \
			-w $windsize
		status=$?; [ ! $status -eq 0 ] && return 1
	else
		new_rm $out_dir/windows.$name.$windsize.bed
		cp $bed_fn $out_dir/windows.$name.$windsize.bed
	fi
	
	# Append GC content
	get_GCcontent -b $bt_dir -f $fasta_fn \
		-i $out_dir/windows.$name.$windsize.bed \
		-o $out_dir/windows.$name.$windsize.gc.bed
	status=$?; [ ! $status -eq 0 ] && return 1
	
	# Get window counts/gc/coverage/num bases covered
	tmp_fn=$out_dir/windows.$name.$windsize.gc.mMAPQ${min_mapq}.bed
	if [ ! -f $tmp_fn ]; then
		echo -e "`date`: Start coverageBed" >&2
		$sam_dir/bin/samtools view -h -b -q $min_mapq \
			$out_dir/$name.sort.bam | $bt_dir/bin/coverageBed \
			-split -a $out_dir/windows.$name.$windsize.gc.bed -b stdin \
			-sorted > $out_dir/windows.$name.$windsize.gc.mMAPQ${min_mapq}.bed
		status=$?; [ ! $status -eq 0 ] \
			&& echo -e "`date`: Error in coverage" >&2 && return 1
		echo -e "`date`: Finish coverageBed" >&2
	fi
	
	# Clean up
	echo -e "`date`: Cleaning up intermediate files" >&2
	new_rm $out_dir/windows.$name.$windsize.bed \
		$out_dir/windows.$name.$windsize.gc.bed
	[ $rm_bam -eq 1 ] \
		&& new_rm $out_dir/$name.sort.bam $out_dir/$name.sort.bam.bai
	
	# Compress
	echo -e "`date`: Compressing final file" >&2
	gzip $out_dir/windows.$name.$windsize.gc.mMAPQ${min_mapq}.bed
	
	return 0
	
}

srcPL_dnaGetCNT=1

## EOF
