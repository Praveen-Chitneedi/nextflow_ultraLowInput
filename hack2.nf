folderIN  = "$baseDir/data/0_fastq"

params.reads = "$folderIN/*_R{1,2}_*.fastq" 
reads_ch = Channel.fromFilePairs(params.reads)

process fastqc {

	echo true
	publishDir "$baseDir/data/1_fastqc", mode: 'copy'

	input:
	tuple val(sample_id), file(reads_file) from reads_ch

    	output:
    	file("fastqc_${sample_id}_logs") into fastqc_ch

	script:

	"""
#	echo $sample_id $reads_file;

	mkdir -p "fastqc_${sample_id}_logs"

	fastqc -o fastqc_${sample_id}_logs -f fastq -q $folderIN/${reads_file}
	
	"""
}


reads_ch = Channel.fromFilePairs(params.reads)
folderOUT = "$baseDir/data/2_cutadapt"
process cutadapt {
	
	echo true
//	publishDir folderOUT, mode: 'copy'

	input:
	tuple val(sample_id), file(reads_file) from reads_ch
	
	output:
//    	tuple val(sample_id), file ('${sample_id}*.fastq') into res1

	exec:
	file1 = reads_file[0]
	file2 = reads_file[1]
//	println("process cutadapt ${sample_id} $file1 $file2")

	script:
	"""
	echo process cutadapt : $sample_id $file1 $file2;	
	mkdir -p $folderOUT
	cutadapt -o $baseDir/data/2_cutadapt/${sample_id}_R1.fastq $folderIN/$file1
	cutadapt -o $baseDir/data/2_cutadapt/${sample_id}_R2.fastq $folderIN/$file2 
#	cutadapt -o ${sample_id}_R1.fastq $folderIN/$file1
#	cutadapt -o ${sample_id}_R2.fastq $folderIN/$file2 
	"""

}

