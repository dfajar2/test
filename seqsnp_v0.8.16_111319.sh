#!/bin/bash
set -eo pipefail

echo $@"
" > date_command.txt


export PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:${PATH}
export PATH=/media/software/cutadapt/1.18/bin:${PATH}
#export PATH=/media/software/bowtie2/2.2.3/bin:${PATH}
export PATH=/media/software/bowtie2/2.3.4.3/bin:${PATH}
export PATH=/media/software/samtools/1.9/bin:${PATH}
export PATH=/media/software/ts/1.0/bin:${PATH}

#Freebayes exported below depending on Version to use
#export PATH= /media/software/freebayes/1.0.2/bin:${PATH}
#export PATH=/media/software/freebayes/1.2.0-2-g29c4002/bin:${PATH}

export PATH=/media/software/bedtools/2.27.1/bin:${PATH}
export PATH=/media/software/bedops/2.4.35/bin:${PATH}
export PATH=/media/software/vcftools/0.1.17/bin:/media/software/vcflib/bin:${PATH}
export PATH=/media/software/vcflib/bin:${PATH}

export PATH=/media/software/bwa/0.7.17/bin:$PATH

# Export PERL5LIB for vcftools
export PERL5LIB=/media/software/vcftools/0.1.17/share/perl/5.26.1:${PERL5LIB}



echo Running seqsnp pipeline
echo ==================
beg=`date +%s`

usage="

USAGE:

    $ /path/to/SeqSNP_script.sh --ngs NGS<number> --reads path/to/fastq_samples_reads(.gz)/directory/ --reference path/to/Reference_genome --targets path/to/SNP_target_file ] ...

    Required:
    --ngs           NGS<number> . It'll be the Output directory name.
    --reads         full_path/to/fastq_samples_reads(.gz)/
                    Compressed reads files directory. If uncompressed files, 
                    also use -u.   
    --reference     path/to/Reference_genome 
                    Reference Genome must be in fasta format (extensions .fa .fasta or .fna)
    --targets       full_path/to/SNP_target_file
                    SNP target file must be in bed format (.bed) specific for the
                    Reference genome.
    
    Optional: 
     --jobs         Num_jobs_for_ts [60] 
     --ploidy       Ploidy [2]
     --step         Selected steps to run [all]:
            all             - Default. Includes trimming, reads alignment, freebayes 
	                      variant calling and report generation, 
	    Trimming        - Trim reads
	    Alignment	    - Bowtie2/BWA alignment
	    VariantAnalysis - Freebayes variant calling 
	    AddUncalled     - Check for SNPs in target list and include
	                      uncalled SNPs 
	    Report          - Generate report from FULL vcf . 
     --only         Use if you ONLY want to run the step decladred with --step. 
                    e.g. '--step Alignment --only' will run the alignment step and exit.
     --stop         Use if you want so stop the analysis after a specific step \(--step\). The 
                    pipeline will start from the beginning of the process.
                    e.g. '--step Alignment --stop' will stop after the sample alignments are done.
     --uncompressed Declare for uncompressed input read files (fastq format). 
                    Default .gz files.
     --use_bwa      Use BWA instead of Bowtie2 for read alignment. Default: Bowtie2
     --large_genome Change index algorithm to bwtsw. Use only with --use_bwa
     --fb_version   Freebayes version.
            100 - Freebayes v1.0.0
	    102 - Freebayes v1.0.2 (Default)
            120 - Freebayes v1.2.0
     --offtarget    Calculate number of offtarget reads not hitting SNPs locations. 
     --S_name_delim Decrapted - It is assumed the name is delimited by the underscore. Sample names
                    will be splitted by '_S'. Default delimiter to extract sample names is '_'.
     --verbose      verbose. Default FALSE.
"
for arg in "$@"; do
     shift
     case "$arg" in
	     "--ngs") 
		     set -- "$@" "-n" ;;
	     "--reads") 
		     set -- "$@" "-f" ;;
	     "--reference") 
		     set -- "$@" "-R" ;;
	     "--targets") 
		     set -- "$@" "-t" ;;
	     "--jobs") 
		     set -- "$@" "-j" ;;
	     "--ploidy") 
		     set -- "$@" "-p" ;;
	     "--step") 
		     set -- "$@" "-s" ;;
	     "--fb_version") 
		     set -- "$@" "-F" ;;
	     "--use_bwa")
		     set -- "$@" "-a" ;;
	     "--large_genome")
		     set -- "$@" "-l" ;;
	     "--only") 
		     set -- "$@" "-o" ;;
	     "--stop")
		     set -- "$@" "-S" ;;
	     "--uncompressed") 
		     set -- "$@" "-u" ;;
	     "--offtarget")
		     set -- "$@" "-g" ;;
	     "--verbose") 
		     set -- "$@" "-v" ;;
             "--S_name_delim") 
		     set -- "$@" "-N" ;;
	     *) 
		     set -- "$@" "$arg"
     esac
done

#n:f:R:t:j:p:s:F:eoudvN
l=0
j=60
p=2
s=all
o=0
S=0
u=0
F=102
a=0
g=0
N=0

while getopts "n:f:R:t:j:p:s:F:oSudglvaN" options; do
	case "${options}" in
		n)
			n=${OPTARG} ;;
		f)
			f=${OPTARG} ;;
		R)
			R=${OPTARG} ;;
		t)
			t=${OPTARG} ;;
		j)
			j=${OPTARG} ;;
		p)
			p=${OPTARG} ;;
		s)
			s=${OPTARG} ;;
		o)
			o=1 ;;
		S)
			S=1 ;;
		u)
		        u=1 ;;
		g)
			g=1 ;;
		F)
		        F=${OPTARG} ;;
		a)
			a=1 ;;
		l)
			l=1 ;;
		v)
			set -x ;;
		N)
			N=1 ;;
		*)
			echo "${usage}" ; exit 1;
	esac
done

shift $((OPTIND-1))
if [ -z "${n}" ] || [ -z "${f}" ] || [ -z "${R}" ] || [ -z "${t}" ] ; then
	echo ; echo "ERROR - Missing arguments"; echo "$usage"; exit 1
fi

echo "Log file: " | tee ${n}.log

if [ "$F" != "100" ] && [ "$F" != "102" ] && [ "$F" != "120" ];
then
	echo ERROR - Invalid Freebayes Version \(-F\). Version ${F} not supported. | tee -a ${n}.log
	echo Valid versions: 100 \(for v1.0.0\), 102 \(for v1.0.2\), 120 \(for v1.2.0\) | tee -a ${n}.log
	echo ; echo -------- | tee -a ${n}.log; echo "$usage"; exit 1; 
fi

if [ "$F" == "100" ];
then
    echo Freebayes v1.0.0 selected | tee -a ${n}.log
    export PATH=/media/software/freebayes/1.0.2/bin:${PATH}
    add_param=" --no-mnps --no-complex "
elif [ "$F" == "102" ];
then
    echo Freebayes v1.0.2 selected | tee -a ${n}.log
    export PATH=/media/software/anaconda2/envs/freebayes1.0.2/bin:${PATH}
    add_param=" --no-mnps --no-complex "
elif [ "$F" == "120" ];
then
    echo Freebayes v1.2.0 selected | tee -a ${n}.log
    export PATH=/media/software/freebayes/1.2.0-2-g29c4002/bin:${PATH}
    add_param=" "
else
	echo ERROR - Invalid Freebayes Version \(-F\). Version ${F} not supported. | tee -a ${n}.log
	echo Valid versions: 100 \(for v1.0.0\), 102 \(for v1.0.2\), 120 \(for v1.2.0\) | tee -a ${n}.log
	echo ; echo -------- | tee -a ${n}.log; echo "$usage"; exit 1; 
fi

NGS=${n}
reads=${f}
reference=${R}
targets=${t}
numjobs=${j}
base_dir=$(pwd) 
ploidy=${p}
ref_path=$(echo $reference | awk -F"/" '{OFS="/"; $(NF--)="";print}')
ref_complete_name=$(echo $reference | awk -F"/" '{print $NF}')
ref_basename=$(echo $reference | awk -F"/" '{print $NF}'| awk -F".fa" '{print $1}')
targets_path=$(echo $targets | awk -F"/" '{OFS="/"; $(NF--)="";print}')
targets_SNP_basename=$(echo $targets | awk -F"/" '{print $NF}'| awk -F".bed" '{print $1}')
echo Files and paths
echo ---------------
echo "base directory       = "${base_dir} | tee -a ${n}.log
echo "NGS experiment       = "${NGS} | tee -a ${n}.log
echo "reads directory      = "${reads} | tee -a ${n}.log
echo "reference genome     = "${reference} | tee -a ${n}.log
echo "targets file         = "${targets} | tee -a ${n}.log
echo "Reference_path       = "${ref_path} | tee -a ${n}.log
echo "Reference_base_name  = "${ref_basename} | tee -a ${n}.log
echo "Targets_path         = "${targets_path} | tee -a ${n}.log
echo "Targets_SNP_basename = "${targets_SNP_basename} | tee -a ${n}.log
echo
echo Additional parameters
echo ---------------------
echo "TS_cores             = "${numjobs} | tee -a ${n}.log
echo "Ploidy               = "${ploidy} | tee -a ${n}.log
#echo -----
#if [ $m -eq 1 ] ;
#then
#	echo Samples will be merged | tee -a ${n}.log
#else
#	echo Samples will not be merged | tee -a ${n}.log
#fi
#echo ----- | tee -a ${n}.log
echo
echo Steps selected to run = ${s} | tee -a ${n}.log

if [[ ${o} -eq 0 ]];
then
	echo Running steps starting from SymLinksData process | tee -a ${n}.log
elif [[ ${o} -eq 1 ]];
then
	echo Running ONLY step ${s} | tee -a ${n}.log
else
	echo ERROR - Wrong -o selection | tee -a ${n}.log
	echo "$usage" | tee -a ${n}.log
	exit 1
fi

if [[ ${S} -eq 1 ]];
then
	echo Analysis WILL stop after running \"${s}\" step | tee -a ${n}.log
elif [[ ${S} -gt 1 ]];
then
	echo ERROR - Wrong -S \(--stop\) selection | tee -a ${n}.log
	echo "$usage" | tee -a ${n}.log
	exit 1
fi

#Check paths
if [ ! -f $reference ] || [ ! -f $targets ] || [ ! -d ${reads} ];
then
	echo ----- | tee -a ${n}.log; echo ERROR !!! File or directories do not exist. Check files and paths. | tee -a ${n}.log; exit 1
fi

if [ "$s" != "all" ] && [ "$s" != "Report" ] && [ "$s" != "AddUncalled" ] && [ "$s" != "VariantAnalysis" ] && [ "$s" != "Alignment" ] && [ "$s" != "Trimming" ];
then
	echo ERROR - Invalid Step \(-s\) Selection. ${s} step do not exist. | tee -a ${n}.log
	echo Valid options: Trimming, Alignment, VariantAnalysis, AddUncalled, Report | tee -a ${n}.log
	echo | tee -a ${n}.log; echo -------- | tee -a ${n}.log; echo "$usage" | tee -a ${n}.log; exit 1; 
fi


# Get sample names from data folder 
echo 
echo Samples
echo -------

if [ "$N" == "0" ]; #If separated by _
then
	echo Samples name delimiter = _ | tee -a ${n}.log
	samples=$(for f in `ls --color=never -1 "${reads}" | grep _R1 | grep -v Unde | grep  -v _I | grep -v config.xml |  grep -v  FastqSummaryF1L1.txt | awk -F"_" '{print $1}' | sort -u` ; do echo $f; done)
elif [ "$N" == "1" ]; #If separated by _S
then
	echo Samples separated by _S | tee -a ${n}.log
	samples=$(for f in `ls --color=never -1 "${reads}" | grep _R1 | grep -v Unde | grep  -v _I | grep -v config.xml |  grep -v  FastqSummaryF1L1.txt | awk -F"_S" '{print $1}' | sort -u` ; do echo $f; done)  
fi
numsamp=$(printf '%s\n' $samples | wc -w)
echo Number of Samples: $numsamp | tee -a ${n}.log
echo Names of the first 10 samples: `for f in {1..10} ; do echo $samples |  awk -v v=${f} '{print $v}' ; done` | tee -a ${n}.log

#Check if SE or PE

if [ "$u" == "0" ];
then
	echo --------------------- | tee -a ${n}.log
	echo Read files are compressed \(.gz\) | tee -a ${n}.log
else
	echo --------------------- | tee -a ${n}.log
	echo Read files are uncompressed \(.fastq\) | tee -a ${n}.log
fi

if ls -1 $reads | grep -q _R2 ;
then
	echo --------------------- | tee -a ${n}.log
	echo "This study contains PE reads" | tee -a ${n}.log
	run=PE
	echo --------------------- | tee -a ${n}.log
else
	echo --------------------- | tee -a ${n}.log
	echo "This study contains SE reads" | tee -a ${n}.log
	run=SE
	echo --------------------- | tee -a ${n}.log
fi

if [ "$a" == "0" ];
then
	echo IMPORTANT - Bowtie will be used for read mapping | tee -a ${n}.log
elif [ "$a" == "1" ] && [ "$l" == "0" ];
then
	echo IMPORTANT - BWA will be used for read mapping | tee -a ${n}.log
	echo IMPORTANT - IS indexing algorithm \(for small genomes\) | tee -a ${n}.log
elif [ "$a" == "1" ] && [ "$l" == "1" ];
then
	echo IMPORTANT - BWA will be used for read mapping | tee -a ${n}.log
	echo IMPORTANT - bwtsw indexing algorithm \(for large genomes\) | tee -a ${n}.log
fi


echo --------------------- | tee -a ${n}.log
echo IMPORTANT - Verify target file, please | tee -a ${n}.log
num_snps=$(wc -l ${targets})
echo SNPs in file: $num_snps | tee -a ${n}.log
echo First lines from ${targets_SNP_basename} file | tee -a ${n}.log
echo --------------------- | tee -a ${n}.log
echo ${targets_path}${targets_SNP_basename} | tee -a ${n}.log
head -n3 ${targets} | tee -a ${n}.log

# Check reference file name
echo
if [[ "$reference" == *fasta  ]]
then
	echo Reference is in fasta format | tee -a ${n}.log
elif [[ "$reference" == *fa  ]]
then
        echo Reference is in fasta format | tee -a ${n}.log
elif [[ "$reference" == *fna  ]]
then
	echo Reference is in fasta format | tee -a ${n}.log
else
        echo ERROR - What the hell??? File does not have .fasta, .fa or .fna suffix | tee -a ${n}.log
        exit 1
fi

#Check if bed positions are present in reference genome
if [ `bedtools getfasta -fi $reference -bed $targets | grep -c ">"` -eq `wc -l $targets | awk '{print $1}'` ];
then
	echo Valid bed file. All positions are present in reference genome | tee -a ${n}.log
else
	echo ERROR - Invalid bed file. Coordinate\(s\) do not exist on reference genome | tee -a ${n}.log
	exit 1
fi

echo 
offtarget=${g}
#echo offtarget: ${g}

if [ "$offtarget" -eq 0 ]
then
	echo IMPORTANT - Offtarget reads counts will be not included in analysis. Offtarget value: ${offtarget} | tee -a ${n}.log
elif [ "$offtarget" -eq 1 ]
then
	echo IMPORTANT - Offtarget reads counts will be included in analysis. Offtarget value: ${offtarget} | tee -a ${n}.log
fi


#Safe stop
echo --------------- | tee -a ${n}.log
echo Check files and paths | tee -a ${n}.log
echo Press Y to continue, any other key to exit. | tee -a ${n}.log
echo --------------- | tee -a ${n}.log
read input
if [ "$input" != "Y" ] && [ "$input" != "y" ]; 
then
	echo Exiting... | tee -a ${n}.log;  
	echo
	exit 0
fi

### # Check for ${NGS} folder - to continue
### if [ -d ${NGS} ]
### then
###	echo ; echo ---------------; echo WARNING; echo; echo ${NGS} folder already exist!!
###	echo Do you want to delete the \"${NGS}\" folder and run the analysis again? \<Type Y to continue\>
###	read input2
###	if [ "$input2" != "Y" ] && [ "$input2" != "y" ];
###	then
###		echo Exiting...
###		exit 1
###	else
###		rm -r $NGS
###	fi
### fi

# Check for ${NGS} folder and create or cd
if [ -d ${NGS} ]
then
	cd $NGS
else
	mkdir $NGS
	cd $NGS
fi

function SymLinksData() {
cd ${base_dir}/${NGS}
	
# Create symbolic links of read files

# Check for ${NGS} folder - to continue
if [ -d ${base_dir}/${NGS}/data ]
then
	echo ; echo ---------------; echo SymLinksData Step; echo WARNING; echo ${base_dir}/${NGS}/data folder already exist!!
	echo Do you want to delete the \"${base_dir}/${NGS}/data\" folder and run \"${s}\" step again? \<Type Y to continue\>
	read input2
	if [ "$input2" != "Y" ] && [ "$input2" != "y" ];
	then
		echo Exiting...
		exit 1
	else
		rm -r ${base_dir}/${NGS}/data
	fi
fi

mkdir data ; cd data
for  file in $(ls --color=never -1 ${reads}); do ln -s ${reads}$file ; done

data_dir=`echo ${base_dir}/${NGS}/data/`
echo Data directory: $data_dir
lst=$(ls -1 $data_dir)

if [ "$u" == "1" ] && [ -d ${base_dir}/${NGS}/trimmed ]
then
	echo Reads already trimmed
	reads=${base_dir}/${NGS}/trimmed/
elif [ $u == "1" ];
then	
	mkdir ../trimmed
	source /media/software/python3_env/bin/activate
	
	if [ "$run" == "PE" ]
	then
		for f in `ls *_R1*q | grep -v Unde | awk -F "_R1" '{print $1}'`; do echo Trimming Sample ${f}; /media/software/cutadapt/2.3/bin/cutadapt -j 24 --nextseq-trim=20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT  ${f}*R1*q ${f}*R2*q --format=fastq -o ../trimmed/${f}_R1_trimmed.fastq -p ../trimmed/${f}_R2_trimmed.fastq > ../trimmed/${f}.log; done
	elif [ "$run" == "SE" ]
	then
		for f in `ls *_R1*q | grep -v Unde | awk -F "_R1" '{print $1}'`; do echo Trimming Sample ${f}; /media/software/cutadapt/2.3/bin/cutadapt -j 24 --nextseq-trim=20 -g AGATCGGAAGAGCACACGTCTGAACTCCAGTCA ${f}*q --format=fastq -o ../trimmed/${f}_R1_trimmed.fastq >../trimmed/${f}.log; done
	fi
	deactivate
	#Define $reads from trimmed folder
	cd ${base_dir}/${NGS}/trimmed/
	reads=${base_dir}/${NGS}/trimmed/
	echo Reads $reads
fi

if [ "$u" == "0" ] && [ -d ${base_dir}/${NGS}/trimmed ]
then
	echo Reads already trimmed
	reads=${base_dir}/${NGS}/trimmed/
elif [ $u == "0" ];
then	
	mkdir ../trimmed
	source /media/software/python3_env/bin/activate
	
	if [ "$run" == "PE" ]
	then
		for f in `ls *_R1*gz | grep -v Unde | awk -F "_R1" '{print $1}'`; do echo Trimming Sample ${f}; /media/software/cutadapt/2.3/bin/cutadapt -j 24 --nextseq-trim=20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT  ${f}*R1*.fastq.gz ${f}*R2*.fastq.gz --format=fastq -o ../trimmed/${f}_R1_trimmed.fastq -p ../trimmed/${f}_R2_trimmed.fastq > ../trimmed/${f}.log; done
	elif [ "$run" == "SE" ]
	then
		for f in `ls *_R1*gz | grep -v Unde | awk -F "_R1" '{print $1}'`; do echo Trimming Sample ${f}; /media/software/cutadapt/2.3/bin/cutadapt -j 24 --nextseq-trim=20 -g AGATCGGAAGAGCACACGTCTGAACTCCAGTCA ${f}*.fastq.gz --format=fastq -o ../trimmed/${f}_R1_trimmed.fastq >../trimmed/${f}.log; done
	fi
	deactivate
	#Define $reads from trimmed folder
	cd ${base_dir}/${NGS}/trimmed/
	reads=${base_dir}/${NGS}/trimmed/
	echo Reads $reads
fi
#echo debug!!
#exit 1
# Get sample names from trimmed data folder
samples=$(for f in `ls --color=never -1 "${reads}" | grep _R1 | grep -v Unde | grep  -v _I | grep -v config.xml |  grep -v  FastqSummaryF1L1.txt | awk -F"_" '{print $1}' | sort -u` ; do echo $f; done) 
echo Samples list: $samples
cd ..
}

function CheckReference(){
cd ${base_dir}/${NGS}

# mkdir reference; cd reference # not needed??? 

# Check reference file name
if [[ "$reference" == *fasta  ]]
then
	echo Reference is in fasta format
elif [[ "$reference" == *fa ]]
then
	echo Reference is in fasta format
elif [[ "$reference" == *fna ]]
then
	echo Reference is in fasta format
else
	echo What the hell??? File does not have .fasta, .fa nor .fna suffix
	exit 1
fi

if [ "$a" == "0" ]; then
	#Check Reference Index
	if [ -f $reference ] && [ -f ${ref_path}${ref_basename}.1.bt2 ]
	then
		echo Reference and Bowtie index already exist.
	elif [ -f $ref_name ]
	then
		echo Reference exist. Building Bowtie2 reference index.
		cd ${ref_path}
		bowtie2-build --threads 32  ${ref_complete_name} ${ref_basename}
		cd -
	else
		echo Fasta file does not exist to build index ; exit 1
	fi
elif [ "$a" == "1" ] && [ "$l" == "0" ]; then
	# Check Reference and BWA index
	if [ -f $reference ] && [ -f ${ref_path}${ref_basename}*.bwt ]
	then
		echo Reference and BWA index already exist.
	elif [ -f $reference ]
	then
		echo Reference exist. Building BWA reference index.
		cd ${ref_path}
		bwa index ${ref_complete_name}
		cd -
	else
		echo Fasta file does not exist to build index ; exit 1
	fi
elif [ "$a" == "1" ] && [ "$l" == "1" ]; then
	# Check Reference and BWA index
	if [ -f $reference ] && [ -f ${ref_path}${ref_basename}*.bwt ]
	then
		echo Reference and BWA index already exist.
	elif [ -f $reference ]
	then
		echo Reference exist. Building BWA reference index.
		cd ${ref_path}
		bwa index -a bwtsw ${ref_complete_name}
		cd -
	else
		echo Fasta file does not exist to build index ; exit 1
	fi
fi
cd ${base_dir}
}

function Bowtie2Alignment() {
cd ${base_dir}/${NGS}

if [ -d ${base_dir}/${NGS}/bowtie2 ] 
then
	echo ; echo ---------------; echo Alignment Step; echo WARNING; echo ${base_dir}/${NGS}/bowtie2 folder already exist!!
	echo Do you want to delete the \"${base_dir}/${NGS}/bowtie2\" folder and run \"${s}\" step again? \<Type Y to continue\>
	read input2
	if [ "$input2" != "Y" ] && [ "$input2" != "y" ];
	then
		echo Exiting...
		exit 1
	else
		rm -r ${base_dir}/${NGS}/bowtie2
	fi 
fi

# Bowtie2 alignment
mkdir bowtie2 ; cd bowtie2
echo Files list
lst=`ls --color=never -1 ${reads}`
echo Checking if SE or PE
if ls -1 $reads | grep -q _R2 ;
then
	echo "Paired data"
	echo set -eo pipefail > Bowtie_bash.sh
	for f in $samples ; do echo echo Aligning Sample ${f} >> Bowtie_bash.sh; echo bowtie2 --threads 48 --rg-id ${f} --rg SM:${f} --rg LB:${f} --rg CN:LGC_Genomics --rg PL:Illumina --dovetail --minins 0 --maxins 1000 -x ${ref_path}${ref_basename} -1 `ls ${reads}${f}*_R1*fastq* | tr "\n" "," | sed 's/,$//g'` -2 `ls ${reads}${f}*_R2*fastq* | tr "\n" "," | sed 's/,$//g'` 2\>${f}.log \| samtools view -Sbu - \| samtools sort -@32 -m 8G \> ${f}_sorted.bam; done >> Bowtie_bash.sh ; 
	echo "ls *.bam > orig_bam.lst ; samtools merge -b orig_bam.lst --threads 60 merged_orig_bam_sorted.bam ; samtools index merged_orig_bam_sorted.bam" >> Bowtie_bash.sh
else
	echo "Unpaired data"
	echo set -eo pipefail > Bowtie_bash.sh
	echo Reads: $reads
	for f in $samples ; do echo echo Aligning Sample ${f} >> Bowtie_bash.sh; echo bowtie2 --threads 48 --rg-id ${f} --rg SM:${f} --rg LB:${f} --rg CN:LGC_Genomics --rg PL:Illumina --dovetail --minins 0 --maxins 1000 -x ${ref_path}${ref_basename} -U `ls ${reads}${f}*_R1*fastq* | tr "\n" "," | sed 's/,$//g'` 2\>${f}.log \| samtools view -Sbu - \| samtools sort -@32 -m 8G \> ${f}_sorted.bam ; done >> Bowtie_bash.sh
	echo "ls *.bam > orig_bam.lst ; samtools merge -b orig_bam.lst --threads 60 merged_orig_bam_sorted.bam ; samtools index merged_orig_bam_sorted.bam" >> Bowtie_bash.sh
#uncomment line below for debuging
#echo; echo STOP!!!!; cat Bowtie_bash.sh; exit 1
fi


if [ -f Bowtie_bash.sh ]                
then
	echo Aligning samples with bowtie2.
	sh Bowtie_bash.sh
	if [ $? -eq 1 ]
	then
		echo "Bowtie alignment error. Exit error 1."
		exit 1
	fi
else
	echo No bowtie file \(Bowtie_bash.sh\) generated.; exit 1
fi
#Exit bowtie2
#uncomment line below for debuging 
#echo; echo STOP!!!!; cat Bowtie_bash.sh; exit 1
cd ..
}

function BWAAlignment() {

cd ${base_dir}/${NGS}
if [ -d ${base_dir}/${NGS}/bwa ]
then
	echo ; echo ---------------; echo Alignment Step; echo WARNING; echo ${base_dir}/${NGS}/bwa folder already exist!!
	echo Do you want to delete the \"${base_dir}/${NGS}/bwa\" folder and run \"${s}\" step again? \<Type Y to continue\>
	read input2
	if [ "$input2" != "Y" ] && [ "$input2" != "y" ];
	then
		echo Exiting...
		exit 1
	else
		rm -r ${base_dir}/${NGS}/bwa
	fi
fi

#BWA alignment
mkdir bwa; cd bwa
# BWA alignment
echo Files list
lst=`ls --color=never -1 ${reads}`

ls ${reads} 
echo Checking if SE or PE
if ls -1 $reads | grep -q _R2 ;
then
	echo "Paired data"
	echo set -eo pipefail > BWA_bash.sh
	for f in $samples ; do echo echo Aligning Sample ${f} >> BWA_bash.sh ; echo bwa mem -U 17 -M -t 56 ${reference} `ls ${reads}${f}*_R1*fastq | tr "\n" "," | sed 's/,$//g'` `ls ${reads}${f}*_R2*fastq | tr "\n" "," | sed 's/,$//g'` 2\>${f}.log \| samtools view -Sbu - \| samtools sort -@32 -m 8G \> ${f}_sorted.bam ; done  >> BWA_bash.sh
	echo "ls *.bam > orig_bam.lst ; samtools merge -b orig_bam.lst --threads 60 merged_orig_bam_sorted.bam ; samtools index merged_orig_bam_sorted.bam" >> BWA_bash.sh
else
	echo "Unpaired data"
	echo set -eo pipefail > BWA_bash.sh
	for f in $samples ; \
	do echo Aligning Sample ${f} >> BWA_bash.sh ; echo bwa mem -U 17 -M -t 56 ${reference} `ls ${reads}${f}*_R1*fastq |  tr "\n" "," | sed 's/,$//g'` 2\>${f}.log \| samtools view -Sbu - \| samtools sort -@32 -m 8G \> ${f}_sorted.bam ; done  >> BWA_bash.sh;
	echo "ls *.bam > orig_bam.lst ; samtools merge -b orig_bam.lst --threads 60 merged_orig_bam_sorted.bam ; samtools index merged_orig_bam_sorted.bam" >> BWA_bash.sh
fi

echo Running BWA scripts

if [ -f BWA_bash.sh ]
then
	echo Aligning samples with BWA
	sh BWA_bash.sh
else
	echo No BWA file \(BWA_bash.sh\) generated.; exit 1
fi
#Exit bwa
cd ..
}

function Alignment(){
if [ "$a" == "0" ]; then
	Bowtie2Alignment
elif [ "$a" == "1" ]; then
	BWAAlignment
fi
}

function VariantAnalysis() {
cd ${base_dir}/${NGS}

#Check if old analysis exists
if [ -d ${base_dir}/${NGS}/freebayes ]
then
	echo ; echo ---------------; echo Alignment Step; echo WARNING; echo ${base_dir}/${NGS}/freebayes folder already exist!!
	echo Do you want to delete the \"${base_dir}/${NGS}/freebayes\" folder and run \"${s}\" step again? \<Type Y to continue\>
	read input2
	if [ "$input2" != "Y" ] && [ "$input2" != "y" ];
	then
		echo Exiting...
		exit 1
	else
		rm -r ${base_dir}/${NGS}/freebayes
	fi
fi

#check if required alignment files exist
if [ "$a" == "0" ]; then
	if [ `ls ${base_dir}/${NGS}/bowtie2/*.bam | wc -l ` -lt 1 ]
	then
		echo No alignment \(BOWTIE2\) files found
		exit 1
	fi
elif [ "$a" == "1" ]; then
	if [ `ls ${base_dir}/${NGS}/bwa/*.bam | wc -l ` -lt 1 ]
	then
		echo No alignment \(BWA\) files found
		exit 1
	fi
fi

# ts Task Spooler setup  
# Setup to run simultaneously $numjobs jobs 
ts -S $numjobs
pwd
mkdir freebayes ; cd freebayes
pwd
mkdir targets ; cd targets
# Split targets
if [ -f $targets ] 
then
	echo Copying target file only SNPs
	cp ${targets} ${targets_SNP_basename}.bed
	cp ${targets_SNP_basename}.bed ${targets_basename}.bed
	head ${targets_basename}.bed       	
	split -l 50 ${targets_basename}.bed
else
	echo No target bed files; exit 1
fi


# Run freebayes
if [ "$a" == "0" ]; then
	if [ -f xaa ]; then
		min_alt_frac=$(echo  1/"${ploidy}"/3 | bc -l )
		echo --min-alternate-fraction ${min_alt_frac}
		for f in `ls x* `; do echo freebayes --targets ${f} -f  ${reference} ${base_dir}/${NGS}/bowtie2/merged_orig_bam_sorted.bam --min-base-quality 20 --min-supporting-allele-qsum 10 --read-mismatch-limit 4 --min-coverage 4 --mismatch-base-quality-threshold 10 --min-alternate-count 2 --report-genotype-likelihood-max --exclude-unobserved-genotypes --genotype-qualities --ploidy ${ploidy} --min-alternate-fraction ${min_alt_frac} --report-monomorphic ${add_param} ${d} \> ../${f}.vcf  \&\& touch ../${f}.done > ${f}_freebayes.sh ; done
		for freebayes_sh in `ls x*_freebayes.sh` ; do ts -E sh $freebayes_sh; done
	else
		echo No split regions files present \(e.g. xaa\) ; exit 1 
	fi
elif [ "$a" == "1" ]; then
	if [ -f xaa ]; then
		min_alt_frac=$(echo  1/"${ploidy}"/3 | bc -l )
		echo --min-alternate-fraction ${min_alt_frac}
		for f in `ls x* `; do echo freebayes --targets ${f} -f  ${reference} ${base_dir}/${NGS}/bwa/merged_orig_bam_sorted.bam --min-base-quality 20 --min-supporting-allele-qsum 10 --read-mismatch-limit 4 --min-coverage 4 --mismatch-base-quality-threshold 10 --min-alternate-count 2 --report-genotype-likelihood-max --exclude-unobserved-genotypes --genotype-qualities --ploidy ${ploidy} --min-alternate-fraction ${min_alt_frac} --report-monomorphic ${add_param} ${d} \> ../${f}.vcf  \&\& touch ../${f}.done > ${f}_freebayes.sh ; done
		for freebayes_sh in `ls x*_freebayes.sh` ; do ts -E sh $freebayes_sh; done
	else
		echo No split regions files present \(e.g. xaa\); exit 1
	fi
fi



#Exit targets
cd ..

# Merging vcf files
while [ `ls -1 *vcf | wc -l` -gt `ls -1 *.done 2>/dev/null | wc -l` ]
do
	echo -en "\rFreebayes running, please wait. `ls -1 *.done 2>/dev/null | wc -l` / `ls -1 *vcf | wc -l` jobs completed < "$((`date +%s` - beg))" sec >" 
done

#Concatenate vcf and filter for SNPs if necessary
if [ -f xaa.done ] 
then
	ls *vcf > list
	vcf-concat -f list  | vcf-sort -p 60 > ${NGS}_vcf_concatenated_sorted_target_filtered.vcf
else
	echo No vcf from subsets \(i.e. xaa.vcf\); exit 1 
fi
}

function AddUncalled() {

if [ -f ${base_dir}/${NGS}/freebayes/${NGS}_concatenated_sorted_target_filtered_uncalled.vcf ]
then
        echo ; echo ---------------; echo Alignment Step; echo WARNING; echo VCF file with uncalled locations already exist!! ${base_dir}/${NGS}/freebayes/${NGS}_concatenated_sorted_target_filtered_uncalled.vcf
        echo Do you want to delete vcf and intermediate files and run \"${s}\" step again? \<Type Y to continue\>
        read input2
        if [ "$input2" != "Y" ] && [ "$input2" != "y" ];
	then
		echo Exiting...
		exit 1
	else
		rm ${base_dir}/${NGS}/freebayes/${NGS}_samples.tsv
		rm ${base_dir}/${NGS}/freebayes/${NGS}_uncalled_positions.bed
		rm ${base_dir}/${NGS}/freebayes/uncalled_bed.tmp
		rm ${base_dir}/${NGS}/freebayes/${NGS}_concatenated_sorted_target_filtered_uncalled.vcf
		rm ${base_dir}/${NGS}/${NGS}_concatenated_sorted_target_filtered_FULL*
	fi
fi	
	
	
cd ${base_dir}/${NGS}/freebayes

#sample names
if [ -f ${NGS}_vcf_concatenated_sorted_target_filtered.vcf ]
then
	grep -m 1 "^#CHROM" ${NGS}_vcf_concatenated_sorted_target_filtered.vcf | cut -f 10- > ${NGS}_samples.tsv
	bedtools intersect -v -a ${targets} -b ${NGS}_vcf_concatenated_sorted_target_filtered.vcf > ${NGS}_uncalled_positions.bed
	bedtools getfasta -fi ${reference} -fo ${NGS}_uncalled_positions.fasta -bed ${NGS}_uncalled_positions.bed
else
	echo Uncalled Step - No vcf file ; exit 1
fi

vcf=${NGS}_vcf_concatenated_sorted_target_filtered.vcf

if [ -f ${NGS}_uncalled_positions.bed ]
then
	echo Uncalled positions bed file present
	echo Intersect
	bedtools intersect -v -a ${targets} -b ${vcf} > uncalled_bed.tmp
	uncalled_bed_file="uncalled_bed.tmp"
	echo create fasta file from uncalled locations
	bedtools getfasta -fi ${reference} -bed uncalled_bed.tmp -fo uncalled_bed.fasta.tmp
	uncalled_ref_bases="uncalled_bed.fasta.tmp"
	# Number_of_samples
	num_samples=`awk '{print NF}' <(grep -m 1 ^#CHROM ${vcf} | cut -f 10-)`
	echo Number of samples $num_samples
	uncalled_ref_bases="uncalled_bed.fasta.tmp"
	first=`cat $uncalled_ref_bases | awk 'BEGIN {RS=""}{if (/^>/) gsub(/\n/," ",$0) ;  print $0 }' | tr ">" "\n" | grep -v ^$ | sed 's/[:,-]/ /g'  | awk '{OFS="\t"; print $ 1,$3,$1"_"$3,$4,".",0,".","DP=0;DPB=0;EPPR=0;GTI=0;MQMR=0;NS=0;NUMALT=0;ODDS=0;PAIREDR=0;PQR=0;PRO=0;QR=0;RO=0;RPPR=0","GT:GQ:DP:DPR:RO:QR:AO:QA"}'`
	# SL3.0ch09       30588048   SL3.0ch09_30588048     A       .       0       .       DP=0;DPB=0;EPPR=0;GTI=0;MQMR=0;NS=0;NUMALT=0;ODDS=0;PAIREDR=0;PQR=0;PRO=0;QR=0;RO=0;RPPR=0      GT:GQ:DP:DPR:RO:QR:AO:QA
	# SL3.0ch09       60566973   SL3.0ch09_60566973     T       .       0       .       DP=0;DPB=0;EPPR=0;GTI=0;MQMR=0;NS=0;NUMALT=0;ODDS=0;PAIREDR=0;PQR=0;PRO=0;QR=0;RO=0;RPPR=0      GT:GQ:DP:DPR:RO:QR:AO:QA
	ploid_tmp=$(coco=$(awk -v ploidy_minus1=$ploidy-1 'BEGIN {while (c++<ploidy_minus1-1) printf "./"}') ; echo $coco)
	second=`seq $num_samples | awk -v PLOID=$ploid_tmp '{print "\t"PLOID".:.:.:.:.:.:.:."}'`
	#        ././././.:.:.:.:.:.:.:.        ././././.:.:.:.:.:.:.:. ././././.:.:.:.:.:.:.:.
	# Create complete vcf lines of uncalled locations
	IFS=$'\n'       # make newlines the only separator
	for j in $first; do echo $j,$second ; done | sed 's/,/\t/g' | sed 's/ //g' | sed 's/\t\t/\t/g' > ${NGS}_concatenated_sorted_target_filtered_uncalled.vcf
	# Join both vcfs and comments
	# Comments
	grep "^#" ${vcf} > ${NGS}_new_vcf.tmp
	# Grep lines from called and uncalled sites with no comments
	cat <(grep -v "^#" ${vcf}) <(cat ${NGS}_concatenated_sorted_target_filtered_uncalled.vcf) | sort -V -k 1,1 -k 2,2n >> ${NGS}_new_vcf.tmp
    	# Rename file
	cp ${NGS}_new_vcf.tmp ${NGS}_concatenated_sorted_target_filtered_FULL.vcf
	vcfbreakmulti ${NGS}_concatenated_sorted_target_filtered_FULL.vcf | vcfallelicprimitives -k -g > ../${NGS}_concatenated_sorted_target_filtered_FULL_allelicprim_breakmulti.vcf
else
	echo No uncalled positions bed file
fi

#Exit freebayes
cd ..
}

function Report() {
fullvcf=${NGS}_concatenated_sorted_target_filtered_FULL_allelicprim_breakmulti.vcf

fullvcf_basename=$(echo ${NGS}_concatenated_sorted_target_filtered_FULL_allelicprim_breakmulti.vcf | awk -F".vcf" '{print $1}')

if [ -f ${fullvcf_basename}.vcf ]
then
	mv ${fullvcf_basename}.vcf temp.recode.vcf
	grep "^#" temp.recode.vcf > ${fullvcf_basename}.vcf
	bedtools intersect -a temp.recode.vcf -b ${targets} >> ${fullvcf_basename}.vcf
	bgzip --threads 8 ${fullvcf_basename}.vcf
	tabix -p vcf ${fullvcf_basename}.vcf.gz
	# Hard filter for at least a depth of eight reads, otherwise "."
	vcftools --gzvcf ${fullvcf_basename}.vcf.gz --minDP 8 --recode --recode-INFO-all --out ${fullvcf_basename}
	mv ${fullvcf_basename}.recode.vcf temp.recode.vcf
	grep "^#" temp.recode.vcf > ${fullvcf_basename}.recode.vcf
	bedtools intersect -a temp.recode.vcf -b ${targets} >> ${fullvcf_basename}.recode.vcf
	paste <(paste <(vcf-query -f 'CHROM:POS\tREF\tALT[\t%SAMPLE"GT"]\n' ${fullvcf_basename}.recode.vcf | head -n1) <(vcf-query -f '[AC %SAMPLE\t]NS\tAF\n' ${fullvcf_basename}.recode.vcf | head -n1)) <(echo; paste <(vcf-query -f '%CHROM:%POS\t%REF\t%ALT[\t%GT]\n' ${fullvcf_basename}.recode.vcf ) <(vcf-query  -f '[%RO\/%AO\t]%INFO/NS\t%INFO/AF\n' ${fullvcf_basename}.recode.vcf )) | sed 's/^\t//g' | sed 's/"\t/\t/g' | sed 's/"/_/g' | sed 's/|/\//g' > ${NGS}_spreadsheet.tab #Added
else
	echo No full vcf file; exit 1
fi


#if [ -f ${fullvcf_basename}.vcf ]
#then
#	bgzip --threads 8 ${fullvcf_basename}.vcf
#	tabix -p vcf ${fullvcf_basename}.vcf.gz
#	# Hard filter for at least a depth of eight reads, otherwise "."
#	vcftools --gzvcf ${fullvcf_basename}.vcf.gz --minDP 8 --recode --recode-INFO-all --out ${fullvcf_basename}
#	paste <(paste <(vcf-query -f 'CHROM:POS\tREF\tALT[\t%SAMPLE"GT"]\n' ${fullvcf_basename}.recode.vcf | head -n1) <(vcf-query -f '[AC %SAMPLE\t]NS\tAF\n' ${fullvcf_basename}.recode.vcf | head -n1)) <(echo; paste <(vcf-query -f '%CHROM:%POS\t%REF\t%ALT[\t%GT]\n' ${fullvcf_basename}.recode.vcf ) <(vcf-query  -f '[%RO\/%AO\t]%INFO/NS\t%INFO/AF\n' ${fullvcf_basename}.recode.vcf )) | sed 's/^\t//g' | sed 's/"\t/\t/g' | sed 's/"/_/g' > ${NGS}_spreadsheet.tab
#else
#	echo No full vcf file; exit 1
#fi
}

function offtarget(){
if [ "$a" == "0" ]; then
	cd ${base_dir}/${NGS}/bowtie2
	mkdir offtarget
elif [ "$a" == "1" ]; then
	cd ${base_dir}/${NGS}/bwa
	mkdir offtarget
fi 

echo Stats merged bam file
samtools flagstat merged*bam > ../merged_bam_stats.txt
	
echo Generating Regions plus / minus 100 bases around SNP
sort-bed ${targets} | bedops --range 100 --everything - > offtarget/${targets_SNP_basename}_paddedSortedRegions.bed

echo Extracting offtarget reads
for f in `ls *bam | grep -v merged | awk -F".bam" '{print $1}'`; do echo Sample ${f} ; bam2bed < ${f}.bam | bedops --not-element-of -100% - offtarget/${targets_SNP_basename}_paddedSortedRegions.bed > offtarget/${f}_offtarget_reads.bed ; done

echo Copying offtarget counts 
cd offtarget
wc -l *_offtarget_reads.bed > offtarget_reads_counts.txt
cp offtarget_reads_counts.txt ${base_dir}/${NGS}
cd ${base_dir}/${NGS}
}


if [ "$s" == "all" ];
then
	SymLinksData
	CheckReference
	Alignment
	VariantAnalysis
	AddUncalled
	Report
elif [ "$s" == "Trimming" ] && [ "$o" -eq 1 ];
then
	SymLinksData
elif [ "$s" == "Trimming" ] && [ "$S" -eq 1 ];
then
	SymLinksData
elif [ "$s" == "Trimming" ] && [ "$o" -eq 0 ];
then
	SymLinksData
	CheckReference
	Alignment
	VariantAnalysis
	AddUncalled
	Report
elif [ "$s" == "Alignment" ] && [ "$o" -eq 1 ]; 
then
	#SymLinksData
	CheckReference
	Alignment
elif [ "$s" == "Alignment" ] && [ "$S" -eq 1 ];
then
	SymLinksData
	CheckReference
	Alignment
elif [ "$s" == "Alignment" ] && [ "$o" -eq 0 ];
then
	SymLinksData
	CheckReference
	Alignment
	VariantAnalysis
	AddUncalled
	Report
elif [ "$s" == "VariantAnalysis" ] && [ "$o" -eq 1 ];
then
	VariantAnalysis
elif [ "$s" == "VariantAnalysis" ] && [ "$S" -eq 1 ];
then
        SymLinksData
	CheckReference
	Alignment
	VariantAnalysis
elif [ "$s" == "VariantAnalysis" ] && [ "$o" -eq 0 ];
then
	VariantAnalysis
	AddUncalled
	Report    
elif [ "$s" == "AddUncalled" ] && [ "$o" -eq 1 ];
then
        AddUncalled
elif [ "$s" == "AddUncalled" ] && [ "$S" -eq 1 ];
then
	SymLinksData
	CheckReference
	Alignment
	VariantAnalysis
        AddUncalled
elif [ "$s" == "AddUncalled" ] && [ "$o" -eq 0 ];
then
	AddUncalled
	Report
elif [ "$s" == "Report" ] && [ "$o" -eq 1 ];
then
	Report
elif [ "$s" == "Report" ] && [ "$S" -eq 1 ]; 
then
	SymLinksData
	CheckReference
	Alignment
	VariantAnalysis
	AddUncalled
	Report
elif [ "$s" != "all" ] && [ "$s" != "Report" ] && [ "$s" != "AddUncalled" ] && [ "$s" != "VariantAnalysis" ] && [ "$s" != "Alignment" ] ;
then
	echo Invalid Step \(-s\) Selection
	exit 1;
fi

if [ "$offtarget" -eq 1 ] ;
then
	offtarget
	paste <(paste <(vcf-query -f 'CHROM:POS\tREF\tALT[\t%SAMPLE"GT"]\n' ${fullvcf_basename}.recode.vcf | head -n1) <(vcf-query -f '[AC_%SAMPLE\t]NS\tDP\tAF\n' ${fullvcf_basename}.recode.vcf | head -n1)) <(echo; paste <(vcf-query -f '%CHROM:%POS\t%REF\t%ALT[\t%GT]\n' ${fullvcf_basename}.recode.vcf ) <(vcf-query  -f '[%RO\/%AO\t]%INFO/NS\t%INFO/DP\t%INFO/AF\n' ${fullvcf_basename}.recode.vcf )) | sed 's/^\t//g' | sed 's/"\t/\t/g' | sed 's/"/_/g' | sed 's/|/\//g' > ${NGS}_spreadsheet_DP.tab
fi

echo
echo - 
echo
cowsay -f $(echo bunny default dragon moose sheep skeleton stegosaurus turkey | sed 's/ /\n/g' | shuf | head -n1) " Done !!! -- Time: $((`date +%s` - beg)) sec)" 

exit 0

