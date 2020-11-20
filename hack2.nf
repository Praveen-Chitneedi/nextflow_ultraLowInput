/* 
 * pipeline input parameters 
 */
folderIN  = "$baseDir/data/0_fastq"
reads     = "$folderIN/*_R{1,2}_*.fastq" 

log.info """\
         R N A S E Q - N F   P I P E L I N E    
         ===================================
         dir          : ${params.dir}
         genome       : ${params.genome}
         """
         .stripIndent()


folderOUT = "$baseDir/data/1_fastqc"

Channel.fromFilePairs(reads).into{reads_ch1; reads_ch2}

process fastqc {

	publishDir folderOUT, mode: 'copy'

	input:
	tuple val(sample_id), file(reads_file) from reads_ch1

    	output:
    	file("fastqc_${sample_id}_logs") into fastqc_ch

	script:
	"""
	mkdir -p "fastqc_${sample_id}_logs"
	fastqc -o fastqc_${sample_id}_logs -f fastq -q $folderIN/${reads_file}
	"""
}


folderOUT = "$baseDir/data/2_cutadapt"

process cutadapt {
	container 'dceoy/cutadapt' //pegi3s/cutadapt'
	echo true
	publishDir "$folderOUT", pattern: '*.fastq', mode: 'copy'

	input:
	tuple val(sample_id), file(reads_file) from reads_ch2
	
	output:
    	tuple val(sample_id), file ("${sample_id}_*.fastq") into res_cutadapt_ch1

	script:
	"""
	echo process cutadapt : $sample_id ${reads_file[0]} ${reads_file[1]}	
	mkdir -p $folderOUT
	cutadapt --quiet -o ${sample_id}_R1.fastq ${reads_file[0]} 
	cutadapt --quiet -o ${sample_id}_R2.fastq ${reads_file[1]}  
	"""
}
res_cutadapt_ch2 = Channel.fromFilePairs("$folderOUT/*_R{1,2}.fastq")


folderOUT = "$baseDir/data/3_fastqc"
process fastqc2 {

	publishDir folderOUT, mode: 'copy'

	input:
	tuple val(sample_id), file(reads_file) from res_cutadapt_ch1

    	output:
    	file("fastqc_${sample_id}_logs")

	script:
	"""
	mkdir -p "fastqc_${sample_id}_logs"
	fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads_file}
	"""
}

index_ch = Channel.fromPath("$folderIN/idx/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex")
folderOUT = "$baseDir/data/4_bwa"
process bwa {
	container 'dceoy/bwa' 
	publishDir folderOUT, mode: 'copy'

	input:
	file(index) from index_ch
	tuple val(sample_id), file(reads_file) from res_cutadapt_ch2

//    	output:
//    	file("${sample_id}.sam") 

	script:
	"""
	bwa mem -M -t 5 $index ${reads_file[0]} ${reads_file[1]} > ${sample_id}.sam
	"""
}

/*
process index {


    input:
    path transcriptome from params.transcriptome

    

    output:
    path 'index' into index_ch

 

    script:

    if( params.aligner == 'salmon' )

    """

    salmon index --threads $task.cpus -t $transcriptome -i index

    """

    else if( params.aligner == 'kallisto' )

    """

    mkdir index

    kalisto index -i index/transcriptome.idx $transcriptome

    """

    else

    throw new IllegalArgumentException("Unknown aligner $params.aligner")

}

params.aligner = "kalisto"
*/
