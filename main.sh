#!/usr/bin/env bash
set -e

map=''
input=''
databasefna='/home/alexey/tax_n_refs/silva_132_97_16S.fna'
databasetax='/home/alexey/tax_n_refs/taxonomy_7_levels.txt'
databasetree='/home/alexey/tax_n_refs/97_otus.tre'
trim_sw='4:12'
trim_minlen='180'
threads='10'


# Parameters module
while [ -n "$1" ]; do
case $1 in
    -h) echo ""
        echo "This script perform 16s Illumina pipeline in a QIIME-manner. It take raw data from remote folder, \
        make trimming (via Trimmomatic) and joining (fastq-join), filter chimeric sequences (vsearch)."
        echo ""
        echo "vsearch, Trimmomatic and fastq-join REQUIRED"
        echo ""
        echo "script parameters:"
        echo ""
        echo "-i - input directory with raw .fastq.gz files. You can catch them from remote source. Make sure, that your \
        file named as <sample_name>_L001_R1_001.fastq.gz (forward) and <sample_name>_L001_R2_001.fastq.gz (reverse)"
        echo "-m - map file in QIIME format. Second column must contain filename prefix (just <sample_name>, without _L001..)"
	echo
        echo "-p - optional -parameter file. It must content full path for reference database and trimming parameters.\
         By default, they are:"
        echo "          Database:   databasefna='/home/alexey/tax_n_refs/silva_132_97_16S.fna'"
        echo "          Taxonomy for database:  databasetax='/home/alexey/tax_n_refs/taxonomy_7_levels.txt'"
        echo "          Database tree:  databasetree='/home/alexey/tax_n_refs/97_otus.tre'"
        echo "          Trimmomatic sliding window:  trim_sw='4:12'"
        echo "          Trimmomatic minimum length: trim_minlen='180'"
        echo "          Threads for chimera finding:    threads='10'"

        exit 0;;
    -i) input=$2
     shift;;
    -m) map=$2
     shift;;
    -p) param=$2
     shift;;
    *) echo "$1 is not a valid parameter. Please, read -h annotation"
      exit 1
esac
shift
done


#export other variables

if [[ -f ${param} ]]; then
	source $param
fi



#check a map as a file

map=`pwd $map`/$map
if [[ ${map} == "" ]]; then
    echo "Please, enter a map name, use -m parameter"
    exit 1
elif [[ ! -f "${map}" ]]; then
    echo "Error: can't find a map - ${map}";
    exit 2
fi

#check an input
if [[ ${input} == "" ]]; then
    echo "Please, enter an input directory, use -i parameter"
    exit 1
elif [[ ! -d "${input}" ]]; then
    echo "Error: wrong input directory - ${input}";
    exit 2
fi



# check a programs accessibility

#if ! `trimmomatic -h &>/dev/null`; then
#    echo "Error: can't find Trimmomatic"
#    exit 3
#fi

if ! `fastq-join -h &>/dev/null`; then
    echo "Error: can't find fastq-join"
    exit 3
fi

if ! `vsearch -h &>/dev/null`; then
    echo "Error: can't find vsearch"
    exit 3
fi



#check files in <input>, and copy them

echo -n 'Copy raw files...'

if [[ -d raw_data ]]; then
    rm -r raw_data/*.*
    else
    mkdir raw_data
fi


if [[ -f "samples.txt" ]]; then
    rm samples.txt
fi

lenMap=`wc -l ${map} | cut -f 1 -d " "`
for (( i=2; i <= ${lenMap}; i++ )); do
    line=`sed -n -e ${i}p ${map}`;
    sampleID=`echo ${line} | awk '{print $1}'`
    prefix=`echo ${line} | awk '{print $2}'`
    echo ${sampleID} >> samples.txt


    if [[ -f "${input}/${prefix}_L001_R1_001.fastq.gz" ]] && [[ -f "${input}/${prefix}_L001_R2_001.fastq.gz" ]]; then
       cp ${input}/${prefix}_L001_R1_001.fastq.gz raw_data/${sampleID}_L001_R1_001.fastq.gz
       cp ${input}/${prefix}_L001_R2_001.fastq.gz raw_data/${sampleID}_L001_R2_001.fastq.gz
       else
       echo
       echo "Error in ${prefix} pair in ${input}: please, check a table or files in folder"
       exit 2
    fi
    done
echo 'done.'




#trimming and joining

echo 'Trimming...'

if [[ -d trimmed ]]; then
    rm -r trimmed/*
    else
    mkdir trimmed
fi
if [[ -d tmp_trimmed ]]; then
    rm -r tmp_trimmed/*
    else
    mkdir tmp_trimmed
fi

for sample in $(cat samples.txt); do
echo ${sample};
dirname="tmp_trimmed/${sample}";
forward="raw_data/${sample}_L001_R1_001.fastq.gz";
reverse="raw_data/${sample}_L001_R2_001.fastq.gz";
mkdir ${dirname};
trimmomatic PE -phred33 2>>trimmed/trimming.log \
                1>/dev/null \
                ${forward} \
                ${reverse} \
                ${dirname}/trimed_paired_forward.fastq.gz \
                ${dirname}/trimed_upaired_forward.fastq.gz \
                ${dirname}/trimed_paired_reverse.fastq.gz \
                ${dirname}/trimed_unpaired_reverse.fastq.gz \
                SLIDINGWINDOW:${trim_sw} MINLEN:${trim_minlen};
fastq-join 1>>trimmed/trimming.log \
                ${dirname}/trimed_paired_forward.fastq.gz \
                ${dirname}/trimed_paired_reverse.fastq.gz \
                -o ${dirname}/trimmed_%.fastq.gz;
cp ${dirname}/trimmed_join.fastq.gz \
    trimmed/${sample}.fastq.gz;
done

rm tmp_trimmed -r
echo 'Trimming done'
echo




#chimera removing

echo 'Chimera removing...'
if [[ -d dechimered ]]; then
    rm -r dechimered/*
    else
    mkdir dechimered
fi
if [[ -d tmp_dechimered ]]; then
    rm -r tmp_dechimered/*
    else
    mkdir tmp_dechimered
fi

for sample in $(cat samples.txt); do
    echo ${sample};
    dirname="tmp_dechimered/${sample}";
    mkdir ${dirname}
    cp trimmed/${sample}.fastq.gz ${dirname}/${sample}.fastq.gz
    gzip -d ${dirname}/${sample}.fastq.gz;
    sed -n '1~4s/^@/>/p;2~4p' ${dirname}/${sample}.fastq > ${dirname}/${sample}.fna
#    convert_fastaqual_fastq.py \    # attention to output filepath - in this way, it will include fastaqual/
#                -c fastq_to_fastaqual\
#                -f ${dirname}/${sample}.fastq\
#                -o ${dirname}/fastaqual;
    vsearch 2>>dechimered/chimera.log --uchime_ref ${dirname}/${sample}.fna\
        --nonchimeras ${dirname}/dechim_${sample}.fna\
        --db ${databasefna} --threads ${threads};
    mv ${dirname}/dechim_${sample}.fna\
        dechimered/${sample}.fna;

done

rm -r tmp_dechimered
echo 'Chimera removing done'


exit 0