#BATCHES = [x for x in range(1,8)]
BATCHES = ["1"]

#WHITELIST_USE=TRUE
WHITELIST_EXPECTCELLS=96
WHITELIST_CORRECT_THRESH=2


GENOMEDIR="/storage/static_data/mm10/"
MAX_THREAD=4
REF_GTF="/storage/static_data/mm10/mm10.gtf"

#BARCODE_FORMAT="CCCCCCNNNNNN"
BARCODE_FORMAT="NNNNNNCCCCCC"
BARCODE_EXTRACT="string"
BARCODES_EVEN="input_dir/barcodes001-096.txt"
BARCODES_ODD ="input_dir/barcodes097-192.txt"

results_dir="results_"+BARCODE_FORMAT
results_sample=results_dir+"/{sample}"

rule all:
   input:
       results_dir + "/complete_count_matrix.tsv"


rule mergeCounts:
    input:
        expand(results_sample + "/6_renamedmatrix/relabelled.matrix", sample = BATCHES)
    output:
        results_dir + "/complete_count_matrix.tsv"
    shell:
        "join -j 1 {input} > {output}"


rule renameCountHeaders:
    input:
        results_sample + "/5_countmatrix/counts.matrix"
    output:
        results_sample + "/6_renamedmatrix/relabelled.matrix"
    run:
        shell('''
        cat {input} | sed 1 s|\b[0-9]+_([ACTG]+)[_0-9.tx]+\b|F{sample}_\1|g > {output}
        ''')


rule countGenesPerCell:
    input:
        results_sample + "/4_sortedbam/sorted.bam"
    output:
        results_sample + "/5_countmatrix/counts.matrix"
    shell:
        "umi_tools counts --per-gene --gene-tag=XT --per-cell -I {input} -S {output}"


rule reindexBam:
    input:
        results_sample + "/3_featurecounts/Aligned.sortedByCoord.out.bam.featureCounts.bam"
    output:
        bam=results_sample + "/4_sortedbam/sorted.bam",
        bai=results_sample + "/4_sortedbam/sorted.bam.bai"
    run:
        shell("samtools sort {input} -o {output.bam}"),
        shell("samtools index {output.bam}")


rule reads2Genes:
    input:
        results_sample + "/2_starmap/Aligned.sortedByCoord.out.bam"
    output:
        bam=results_sample + "/3_featurecounts/Aligned.sortedByCoord.out.bam.featureCounts.bam",
        cm=results_sample + "/3_featurecounts/counts.tsv"
    params:
        threads=MAX_THREAD,
        gtf=REF_GTF,
    run:
        shell("featureCounts -a {params.gtf} -o {output.cm}  -R BAM {input} -T {params.threads}"),
        shell("mv Aligned.sortedByCoord.out.bam.featureCounts.bam {output}")


rule mapReads:
    input:
        results_sample + "/1_umi_extract/R2_extracted.fastq"
    output:
        results_sample + "/2_starmap/Aligned.sortedByCoord.out.bam"
    params:
        gendir=GENOMEDIR,
        threads=MAX_THREAD,
        outfilterMMMax=1,
        readcommand="zcat",
        outSAMtype="BAM"
    run:
        shell('''
        STAR
        --readFilesIn {input}
        --runThreadN {params.threads}
        --genomeDir {params.gendir}
        --readFilesCommand {params.readcommand}
        --outFilterMultimapNmax {params.outfilterMMMax}
        --outSAMtype {params.outSAMtype} SortedByCoordinate
        '''),
        shell("mv Aligned.sortedByCoord.out.bam {output}")



rule whitelistBCsAndUmis:
    input:
        read1="input_dir/FACS{sample}_R1.fastq",
        read2="input_dir/FACS{sample}_R2.fastq",
        bar_even=results_sample + "/barcode_even.extracted",
        bar_odd =results_sample + "/barcode_odd.extracted"
    output:
        whitelist=results_sample + "/0_umi_whitelist/guessed_barcodes",
        log=      results_sample + "/0_umi_whitelist/whitelist.log"
    params:
        plotprefix=results_sample + "/plots",
        expect_cells=WHITELIST_EXPECTCELLS,
        correct_thresh=WHITELIST_CORRECT_THRESH,
        bc_pattern=BARCODE_FORMAT,
        extract_method=BARCODE_EXTRACT
    shell:
        '''
        umi_tools whitelist\
            --extract-method='{params.extract_method}'\
            --bc-pattern='{params.bc_pattern}'\
            --set-cell-number={params.expect_cells}\
            --error-correct-threshold={params.correct_thresh}\
            --plot-prefix='{params.plotprefix}'\
            --method=umis\
            --stdin='{input.read1}'\
            --log='{output.log}'\
           --stdout {output.whitelist}
        '''


rule measureBarcodeOverlapWithWhitelist:
    input:
        guessed=results_sample + "/0_umi_whitelist/guessed_barcodes",
        bar_even=results_sample + "/__preproc/barcode_even.extracted",
        bar_odd =results_sample + "/__preproc/barcode_odd.extracted"
    output:
        results_sample + "/0_umi_whitelist/stats.txt"
    run:
        bar_used=input.bar_even if int(wildcards.sample)%2==0 else input.bar_odd

        designed_barcodes = {}

        MATCH_MAIN='matches_ḿain'
        MATCH_OTHER='matches_other'
        
        with open(bar_used, 'r') as f:
            for line in f:
                #num, barcode = line.split()
                barcode = line.strip()
                if barcode not in designed_barcodes:
                    designed_barcodes[barcode] = {
                        MATCH_MAIN : [],
                        MATCH_OTHER: []
                    }
                else:
                    print("Duplicate:", barcode, file=sys.stderr)


        with open(input.guessed, 'r') as f:
            for line in f:
                tokens = line.split()
                barcode_main = tokens[0].strip()
                barcode_others = map(lambda x: x.strip(), tokens[1].split(','))

                # Direct match
                if barcode_main in designed_barcodes:
                    designed_barcodes[barcode_main][MATCH_MAIN].append(barcode_main)

                # Matches with others
                for bar_other in barcode_others:
                    if bar_other in designed_barcodes:
                        # This other barcode exists in the designed barcodes, heres the
                        # main barcode this other barcode it belongs to
                        designed_barcodes[bar_other][MATCH_OTHER].append(barcode_main)


        with open(output[0], 'w') as fout:
            # Process map
            print("Wanted\tDirectMatch\tMergedWithOther", file=fout)
            
            for wanted_barcode in designed_barcodes:
                matches_main  = designed_barcodes[wanted_barcode][MATCH_MAIN]
                match_w_other = designed_barcodes[wanted_barcode][MATCH_OTHER]

                print("%s\t%s\t%s" % (
                    wanted_barcode, ','.join(matches_main), ','.join(match_w_other)
                ), file=fout)
        fout.close()


rule extractBCsAndUmis:
    input:
        stats=results_sample + "/0_umi_whitelist/stats.txt",
        # stats are not a required input, but good to be generated here
        read1="input_dir/FACS{sample}_R1.fastq",
        read2="input_dir/FACS{sample}_R2.fastq",
        whitelist=results_sample + "/0_umi_whitelist/guessed_barcodes"
    output:
        read1=results_sample + "/1_umi_extract/R1_extracted.fastq",
        read2=results_sample + "/1_umi_extract/R2_extracted.fastq",
        log=results_sample + "/1_umi_extract/log.txt"
    params:
        bc_pattern=BARCODE_FORMAT,
        extract_method=BARCODE_EXTRACT
    shell:
        '''
        umi_tools extract
        --error-correct-cell
        --filter-cell-barcode
        --whitelist='{input.whitelist}'
        --extract-method='{params.extract_method}'
        --bc-pattern='{params.bc_pattern}'
        --stdin='{input.read1}'
        --read2-in='{input.read2}'
        --stdout='{output.read1}'
        --read2-out='{output.read2}'
        --log='{output.log}'
        '''


rule extractBarcodes:
    input:
        odd=BARCODES_ODD,
        even=BARCODES_EVEN
    output:
        odd= results_sample + "/__preproc/barcode_odd.extracted",
        even=results_sample + "/__preproc/barcode_even.extracted"
    run:
        shell("cut -f 2 {input.odd} > {output.odd}"  ),
        shell("cut -f 2 {input.even} > {output.even}")

