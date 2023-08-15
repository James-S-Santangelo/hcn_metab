######################
#### READ MAPPING ####
######################

rule star_genome_build:
    input:
        ref = config["ref"]
    output:
        directory(f"{STAR_DIR}/star_genome_build")
    log: f"{LOG_DIR}/star/star_genome_build.log"
    conda: "../envs/deg_analysis.yaml"
    threads: 8
    shell:
        """
        mkdir {output}
        STAR --runMode genomeGenerate \
            --genomeDir {output} \
            --genomeFastaFiles {input.ref} \
            --runThreadN {threads} &> {log}
        """

rule star_align_reads:
    input:
        star_build = rules.star_genome_build.output
        R1 = glob.glob(f"{config['raw_reads']}/*.{wildcards.sample}_R1.fastq.gz") 
        R2 = glob.glob(f"{config['raw_reads']}/*.{wildcards.sample}_R2.fastq.gz") 
    output:
        f"{STAR_DIR}/alignments/{{sample}}/{{sample}}_Aligned.sortedByCoord.out.bam"
    log: f"{LOG_DIR}/star/{{sample}}_star_align.log"
    conda: "../envs/deg_analysis.yaml"
    threads: 6
    params:
        out = f"{STAR_DIR}/alignments/{{sample}}/{{sample}}_"
    shell:
        """
        STAR --readFilesIn {input.R1} {input.R2} \
            --outSAMtype BAM SortedByCoordinate \
            --genomeDir {input.star_build} \
            --runThreadN {threads} \
            --readFilesCommand zcat \
            --alignIntronMax 100000 \
            --outFilterMatchNminOverLread 0.4 \
            --outFilterScoreMinOverLread 0.4 \
            --outReadsUnmapped Fastx \
            --outFileNamePrefix {params.out} &> {log} 
        """
