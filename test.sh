#!/bin/bash

run_umi=Y
run_demulti=N

source activate dmq

batch_num=$1    
input_read1=$2   #inputs/FACS1_R1.fastq
input_read2=$3   #inputs/FACS1_R2.fastq
barcodes=$4      #inputs/celseq_oddbatch.txt

work_dir=results/$(date +%Y%m%d_%H%M)/$batch_num/
mkdir $work_dir

echo "Umi-tools"

# CCCCCCXXXXXX by itself does not do anything, needs Ns
if [ "$run_umi" = "Y" ]; then
    extracted_1=$work_dir/`basename $input_read1 .fastq`.extracted.fastq
    extracted_2=$work_dir/`basename $input_read2 .fastq`.extracted.fastq

    corrected_bcs=$work_dir/"corrected_bcs.txt"
    cat $barcodes | awk '{print $2}' > $corrected_bcs
    
    umi_tools\
        extract\
        --extract-method='string'\
        --bc-pattern='CCCCCCNNNNNN'\
        --filter-cell-barcode\
        --whitelist=$corrected_bcs\
        --stdin=$input_read1   --read2-in=$input_read2\
        --stdout $extracted_1  --read2-out=$extracted_2\
        --log=$logdate/umi_tools.log

    #echo "Running stats"
    #./stats.sh $input_read1 $extracted_1 $work_dir
fi



if [ "$run_demulti" = "Y" ]; then
    echo ""
    echo "Je-suite"

    je demultiplex\
       F1=$extracted_1 F2=$extracted_2\
       BARCODE_FILE=$barcodes\
       BARCODE_READ_POS=READ_1\
       BARCODE_FOR_SAMPLE_MATCHING=READ_1\
       SAME_HEADERS=false\
       C=true ADD=true\
       MM=1 MMD=1 Q=10\
       QUALITY_FORMAT=null\
       XT=0 ZT=0 RCHAR=: GZ=false\
       OUTPUT_DIR=results\
       KEEP_UNASSIGNED_READ=true\
       STATS_ONLY=false\
       METRICS_FILE_NAME=$work_dir/jesuite.log

    echo ""
    echo "FINIT"
fi
    
source deactivate dmq

