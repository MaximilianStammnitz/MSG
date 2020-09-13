#!/bin/bash

# MSG: streamlined (1) Manta, (2) svimmer and (3) GraphTyper2-based SV calling and genotyping
# Max Stammnitz, Transmissible Cancer Group, University of Cambridge
# 2020

# MSG

# INPUT
# -d: path to list with absolute paths to input BAMs for SV discovery (accompanied by BAI or CSI indices)
# -g: path to list with absolute paths to additional input BAMs for SV genotyping (accompanied by BAI or CSI indices)
# -r: path to reference genome FASTA
# -o: path to global output folder
# -w: path to file of region windows in chr:start-end format
# -c: number of CPUs (processes) for Manta and Graphtyper2

VERSION=1.1


####################################### FUNCTIONS #######################################

# AUXILIARY FUNCTIONS
# print_help()
# Prints a help guide if -h, or no arguments, are input
print_help() {
    echo
    echo
    echo "| MSG"
    echo "| (1) Manta, (2) svimmer and (3) GraphTyper2-based SV calling and genotyping"
    echo "| Version $VERSION"
    echo "|"
    echo "| Required input:"
    echo "|    -d  List with absolute paths to input BAMs for SV calling (accompanied by BAI or CSI indices)."
    echo "|    -g  List with absolute paths to additional input BAMs for SV genotyping (accompanied by BAI or CSI indices)."
    echo "|    -r  Absolute path to reference genome FASTA file (accompanied by FAI index)."
    echo "|    -o  Absolute path to the output folder."
    echo "|    -w  Absolute path to file of regions to use, one per line in CHR:START-END format."
    echo "|    -c  Number of CPUs (processes) for Manta and GraphTyper2."
    echo "|"
    echo "| Options:"
    echo "|    -h  Print this usage information and exit."
    echo "|    -v  Print version and exit."
    echo "|"
    echo "| Usage:"
    echo "|    MSG -d /path/to/discovery_cohort.txt -g /path/to/genotyping_cohort.txt -r /path/to/reference_genome.fa -o /path/to/output_directory -w /path/to/regions.txt -c cores <INT>"
    echo
    echo
}

# check_file()
# Checks if a file exists and is not empty. In that case, it displays an error message and exits
# Used for checking the output of each step
check_file() {

    if [ ! -s $1 ]; then
        echo -e "\nERROR: Output file $1 was not correctly generated. Please check the logs folder for more information.\n" >&2
        exit 1
    fi

}

# PIPELINE STEPS
# (1) individual_calling()
# Calls variants using Manta from each tumour sample, individually
individual_calling() {

    # Create directories for output files and logs
    mkdir -p $OUTDIR/1_individual_calls
    while read FILE; do

        NAME=`basename $FILE`
        NAME=`echo ${NAME} | cut -d'.' -f1`
        echo -e "\nCalling on $NAME\n"

        # Create sample-specific directory
        mkdir -p $OUTDIR/1_individual_calls/$NAME

        # Configure Manta run with the default settings

	## make a BED file from the regions list
	mkdir -p $OUTDIR/regions
	cp $REGIONS $OUTDIR/regions/regions.txt
	cut -f 1 -d ':' $OUTDIR/regions/regions.txt > $OUTDIR/regions/regions1
	cut -f 2 -d ':' $OUTDIR/regions/regions.txt | cut -f 1 -d '-' > $OUTDIR/regions/regions2
	cut -f 2 -d ':' $OUTDIR/regions/regions.txt | cut -f 2 -d '-' > $OUTDIR/regions/regions3
	paste $OUTDIR/regions/regions1 $OUTDIR/regions/regions2 $OUTDIR/regions/regions3 > $OUTDIR/regions/regions.bed
	bgzip -c $OUTDIR/regions/regions.bed > $OUTDIR/regions/regions.bed.gz
	tabix -C -p bed $OUTDIR/regions/regions.bed.gz
	rm $OUTDIR/regions/regions1 $OUTDIR/regions/regions2 $OUTDIR/regions/regions3

	## run Manta with the regions list
        configManta.py \
        --tumorBam $FILE \
        --referenceFasta $REFERENCE \
        --runDir $OUTDIR/1_individual_calls/$NAME \
	--callRegions $OUTDIR/regions/regions.bed.gz

        # Run Manta with the default settings
        $OUTDIR/1_individual_calls/$NAME/runWorkflow.py -j $CPUS

    done < $CALL

}

# (2) merge_calls()
# Merge Manta SV from all calls, using svimmer
merge_calls() {

    # Create directory for output split files
    mkdir -p $OUTDIR/2_merge

    # Copy, rename, INV-convert, gzip and index each tumour_SV.vcf file within the same folder
    for FILE in `ls -1 $OUTDIR/1_individual_calls`; do

        cp $OUTDIR/1_individual_calls/$FILE/results/variants/tumorSV.vcf.gz $OUTDIR/2_merge/${FILE}_manta.vcf.gz
	gunzip $OUTDIR/2_merge/${FILE}_manta.vcf.gz

	## convert INV format from Manta v1.6 outputs to that of previous versions
	SAMT=$(which samtools)
	convertInversion.py $SAMT $REFERENCE $OUTDIR/2_merge/${FILE}_manta.vcf > $OUTDIR/2_merge/${FILE}_manta_inv.vcf
	mv $OUTDIR/2_merge/${FILE}_manta_inv.vcf $OUTDIR/2_merge/${FILE}_manta.vcf

	bgzip -c $OUTDIR/2_merge/${FILE}_manta.vcf > $OUTDIR/2_merge/${FILE}_manta.vcf.gz
        tabix -C -p vcf $OUTDIR/2_merge/${FILE}_manta.vcf.gz

    done

    # Create input file list for svimmer
    ls -1 $OUTDIR/2_merge/*.gz > $OUTDIR/2_merge/input_merging.txt

    # Provide an input list of chromosomes to merge
    cut -f 1 -d ':' $OUTDIR/regions/regions.txt | uniq |  tr '\n' ' ' > $OUTDIR/regions/chroms.txt

    # Merge SV calls with svimmer, only using the specified genome regions
    svimmer $OUTDIR/2_merge/input_merging.txt \
    --threads $CPUS \
    --ids \
    --output $OUTDIR/2_merge/SVs_manta_merged.vcf \
    --max_distance 500 \
   $(<$OUTDIR/regions/chroms.txt)

}

# (3) genotype_calls()
# Genotype the svimmer output VCF across all tumours and normals, with GraphTyper2
genotype_calls() {

    # Create directory
    mkdir -p $OUTDIR/3_genotyped

    # Sort, gzip and index input VCF
    cat $OUTDIR/2_merge/SVs_manta_merged.vcf | vcf-sort > $OUTDIR/2_merge/SVs_manta_merged.sorted.vcf
    bgzip -c $OUTDIR/2_merge/SVs_manta_merged.sorted.vcf > $OUTDIR/2_merge/SVs_manta_merged.sorted.vcf.gz
    bcftools index -c $OUTDIR/2_merge/SVs_manta_merged.sorted.vcf.gz

    # Launch GraphTyper2
    graphtyper genotype_sv \
    $REFERENCE \
    $OUTDIR/2_merge/SVs_manta_merged.sorted.vcf.gz \
    --sams=$GENO \
    --output=$OUTDIR/3_genotyped \
    --region_file=$REGIONS \
    --csi \
    --verbose \
    --threads=$CPUS

    # Summarise chromosome-wise results
    cut -f 1 -d ':' $OUTDIR/regions/regions.txt | uniq > $OUTDIR/regions/chroms.list
    wc $OUTDIR/regions/chroms.list | cut -f 1 -d ' ' > $OUTDIR/regions/nchroms

    if [ $OUTDIR/regions/nchroms = 1 ]; then

        while read chrom; do
	   zcat $OUTDIR/3_genotyped/${chrom}/*.gz > $OUTDIR/4_SVs_final.vcf
	done < $OUTDIR/regions/chroms.list
        rm -r $OUTDIR/regions

   else

	head -n 1 $OUTDIR/regions/chroms.list > $OUTDIR/regions/chroms1
        while read chrom; do
	   zcat $OUTDIR/3_genotyped/${chrom}/*.gz > $OUTDIR/4_SVs_final.vcf
        done < $OUTDIR/regions/chroms1

	sed '1d' $OUTDIR/regions/chroms.list > $OUTDIR/regions/chroms2
	while read chrom; do
	    bcftools concat $OUTDIR/4_SVs_final.vcf $OUTDIR/3_genotyped/${chrom}/*.gz >> $OUTDIR/4_SVs_final_n.vcf
	done < $OUTDIR/regions/chroms2

	mv $OUTDIR/4_SVs_final_n.vcf $OUTDIR/4_SVs_final.vcf
	rm -r $OUTDIR/regions
    fi

}


################################### END OF FUNCTIONS ####################################


# Check that dependencies (Manta, svimmer, Graphtyper, tabix, vcf-sort, bgzip) are installed
hash configManta.py 2>/dev/null || { echo -e "\nERROR: configManta.py command not found. Please install Manta (https://github.com/Illumina/manta) and add its ../bin directory to your PATH.\n" >&2; exit 1; }
hash svimmer 2>/dev/null || { echo -e "\nERROR: svimmer command not found. Please install svimmer (https://github.com/DecodeGenetics/svimmer) and add its directory to your PATH.\n" >&2; exit 1; }
hash graphtyper 2>/dev/null || { echo -e "\nERROR: graphtyper command not found. Please install GraphTyper2 (https://github.com/DecodeGenetics/graphtyper) and add its directory to your PATH.\n" >&2; exit 1; }
hash tabix 2>/dev/null || { echo -e "\nERROR: tabix: command not found. Please install the SAMtools tabix package and add its directory to your PATH.\n" >&2; exit 1; }
hash bgzip 2>/dev/null || { echo -e "\nERROR: bgzip: command not found. Please install the SAMtools tabix package and add its directory to your PATH.\n" >&2; exit 1; }
hash vcf-sort 2>/dev/null || { echo -e "\nERROR: vcf-sort: command not found. Please install VCFtools and add its directory to your PATH.\n" >&2; exit 1; }

# If no arguments (or -h): print help
if [ "$#" -eq 0 ]; then
    print_help
    exit 0
fi

# Parse input
CALL=""
GENO=""
REFERENCE=""
REGIONS="no"
OUTDIR=""
CPUS=1
EXTRA=""
while getopts "d:g:r:w:o:c:p:hv?" OPT; do
  case $OPT in
    d)
      CALL=$OPTARG
      ;;
    g)
      GENO=$OPTARG
      ;;
    r)
      REFERENCE=$OPTARG
      ;;
    w)
      REGIONS=$OPTARG
      ;;
    o)
      OUTDIR=$OPTARG
      ;;
    c)
      CPUS=$OPTARG
      ;;
    p)
      EXTRA=$OPTARG
      ;;
    h)
      print_help
      exit 0
      ;;
    v)
      echo "MSG $VERSION"
      exit 0
      ;;
    \?)
      print_help
      echo -e "Invalid option: -$OPTARG\n" >&2
      exit 1
      ;;
  esac
done

# Check that all mandatory inputs are present
if [ -z "$CALL" ]; then
   print_help
   echo -e "Input BAM files list for SV discovery (-d) is required\n" >&2
   exit 1
fi

if [ -z "$GENO" ]; then
   print_help
   echo -e "Input BAM files list for SV genotyping (-g) is required\n" >&2
   exit 1
fi

if [ -z "$REFERENCE" ]; then
   print_help
   echo -e "Reference genome FASTA file (-r) is required\n" >&2
   exit 1
fi

if [ -z "$OUTDIR" ]; then
   print_help
   echo -e "Path to the output folder (-o) is required\n" >&2
   exit 1
fi

# Sanity checks:
# Check that: Input files exist; FASTA is indexed; there are input BAMs; all BAMs are indexed; CPUS >0;
# extra Platypus options do not contain options used within the pipeline
if [ ! -s $REFERENCE ]; then
    echo -e "\nERROR: Reference FASTA file not found or empty. Please check the path.\n" >&2
    exit 1
fi

if [ ! -s ${REFERENCE}.fai ]; then
    echo -e "\nERROR: Reference FASTA file has no accompanying .fai index file. Please index the FASTA file using 'samtools faidx' or equivalent.\n" >&2
    exit 1
fi

if  [ "$REGIONS" != "no" ] && [ ! -s $REGIONS ]; then
    echo -e "\nERROR: Regions file not found or empty. Please check the path.\n" >&2
    exit 1
fi

if ! [[ $CPUS =~ ^[1-9]+[0-9]*$ ]]; then
    echo -e "\nERROR: Number of CPUs must be greater than 0\n" >&2
    exit 1
fi

if [ ! -s $CALL ]; then
    echo -e "\nERROR: $CALL list not found. Please check the path.\n" >&2
    exit 1
fi

if [ ! -s $GENO ]; then
    echo -e "\nERROR: $GENO list not found. Please check the path.\n" >&2
    exit 1
fi

FILES=`grep '.bam$' $CALL 2> /dev/null | wc -l`
if [ "$FILES" -eq 0 ]; then
    echo -e "\nERROR: $CALL contains no files with .bam extension. Please check the path.\n" >&2
    exit 1
fi

FILES=`grep '.bam$' $GENO 2> /dev/null | wc -l`
if [ "$FILES" -eq 0 ]; then
    echo -e "\nERROR: $GENO contains no files with .bam extension. Please check the path.\n" >&2
    exit 1
fi

while read FILE; do
    if [ ! -s ${FILE}.bai ] && [ ! -s ${FILE}.csi ]; then
        echo -e "\nERROR: The file $FILE has no accompanying .bai or .csi index file. Please index the file using 'samtools sort' and 'samtools index' or equivalent.\n" >&2
        exit 1
    fi
done < $CALL

while read FILE; do
    if [ ! -s ${FILE}.bai ] && [ ! -s ${FILE}.csi ]; then
        echo -e "\nERROR: The file $FILE has no accompanying .bai or .csi index file. Please index the file using 'samtools sort' and 'samtools index' or equivalent.\n" >&2
        exit 1
    fi
done < $GENO

# START RUNNING
# Copy all standard out and standard error to log file
mkdir -p $OUTDIR/logs
exec &> >(tee -ia $OUTDIR/logs/MSG_`date +"%y%m%d%H%M"`.log)

echo -e "\nThis is MSG $VERSION\n"
echo "Input BAM list for SV calling:      $CALL"
echo "Input BAM list for SV genotyping:   $GENO"
echo "Input reference genome:             $REFERENCE"
echo "Input regions file:                 $REGIONS"
echo "Output directory:                   $OUTDIR"
echo "Number of CPUs to use:              $CPUS"

# Check if there is a checkpoint file from a previous run in the output folder
STEP=0
if [ -s $OUTDIR/logs/CHECKPOINT ]; then
    CHK=`tail -1 $OUTDIR/logs/CHECKPOINT`
    STEP=`echo $CHK | cut -f1 -d" "`
    STEPNAME=`echo $CHK | cut -f2 -d" "`
    echo -e "\n*CHECKPOINT FILE FOUND*"
    echo -e "Resuming execution after last completed step: $STEPNAME"
fi

echo -e "\nExecution started on `date`"

# Each step is performed only if its index is higher than STEP (last finished step index)
# 1. Run Manta individually on every tumour sample
if [ "$STEP" -lt 1 ]; then

    echo -e "\n(1) RUNNING MANTA INDIVIDUALLY ON EVERY SAMPLE OF THE DISCOVERY COHORT"
    individual_calling

    # Check successful execution
    for FILE in `ls -1 $OUTDIR/1_individual_calls`; do
        check_file $OUTDIR/1_individual_calls/$FILE/results/variants/tumorSV.vcf.gz
    done < $CALL

    # Update checkpoint file
    echo -e "\nSuccess"
    echo "1 individual_calling" >> $OUTDIR/logs/CHECKPOINT

fi

# 2. Merge all Manta tumour-specific SV calls, with svimmer
if [ "$STEP" -lt 2 ]; then

    echo -e "\n(2) RUNNING SVIMMER TO MERGE ALL DISCOVERED SVS"
    merge_calls

    # Check successful execution
    check_file $OUTDIR/2_merge/SVs_manta_merged.vcf

    # Update checkpoint file
    echo -e "\nSuccess"
    echo "2 merge_calls" >> $OUTDIR/logs/CHECKPOINT

fi

# 3 Genotype the svimmer output VCF across all tumours and normals, with GraphTyper2
if [ "$STEP" -lt 3 ]; then

    echo -e "\n(3) RUNNING GRAPHTYPER2 TO GENOTYPE SVS ACROSS THE GENOTYPING COHORT"
    genotype_calls

    # Check successful execution
    check_file $OUTDIR/4_SVs_final.vcf

    # Update checkpoint file
    echo -e "\nSuccess"
    echo "3 genotype_calls" >> $OUTDIR/logs/CHECKPOINT

fi

echo -e "\nExecution finished on `date`"
echo -e "\nALL DONE!"
echo -e "Output is in: $OUTDIR/4_SVs_final.vcf\n"
