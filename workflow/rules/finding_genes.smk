rule star_genome_build_olsen:
    input:
        ref = config["olsen"]
    output:
        directory(f"{STAR_DIR}/star_genome_build_olsen")
    log: f"{LOG_DIR}/star_olsen/star_genome_build_olsen.log"
    conda: "../envs/finding_genes.yaml"
    threads: 8
    shell:
        """
        mkdir {output}
        STAR --runMode genomeGenerate \
            --genomeDir {output} \
            --genomeFastaFiles {input.ref} \
            --runThreadN {threads} &> {log}
        """

rule star_align_reads_toOlsen:
    input:
        star_build = rules.star_genome_build_olsen.output,
        R1 = lambda w: glob.glob(f"{config['raw_reads']}/*.{w.sample}_R1.fastq.gz"),
        R2 = lambda w: glob.glob(f"{config['raw_reads']}/*.{w.sample}_R2.fastq.gz") 
    output:
        f"{STAR_DIR}/alignments/olsen/{{sample}}/olsen_{{sample}}_Aligned.sortedByCoord.out.bam"
    log: f"{LOG_DIR}/star_olsen/{{sample}}_star_align.log"
    conda: "../envs/finding_genes.yaml"
    threads: 6
    params:
        out = f"{STAR_DIR}/alignments/olsen/{{sample}}/olsen_{{sample}}_"
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

rule minimap_align_hcn_genes:
    input:
        genes = config["hcn_genes"],
        ref = lambda w: config["olsen"] if w.ref == "olsen" else config["ref"]
    output:
        bam = f"{MINIMAP_DIR}/{{ref}}/hcn_genes_to{{ref}}.bam",
        bai = f"{MINIMAP_DIR}/{{ref}}/hcn_genes_to{{ref}}.bam.bai"
    conda: "../envs/finding_genes.yaml"
    log: f"{LOG_DIR}/minimap/{{ref}}_minimap_hcn_genes.log"
    threads: 8 
    shell:
        """
        ( minimap2 -ax asm5 {input.ref} {input.genes} \
            -t {threads} | samtools sort -@ {threads} - -o {output.bam} &&\
            samtools index {output.bam} ) 2> {log}
        """

rule subset_bam:
    input:
        rules.star_align_reads.output.bam
    output:
        bam = f"{READ_ANALYSIS_DIR}/{{sample}}_{{gene}}.bam",
        bai = f"{READ_ANALYSIS_DIR}/{{sample}}_{{gene}}.bam.bai"
    conda: "../envs/finding_genes.yaml"
    params:
        region = lambda w: "Chr02_Occ:9002690-9004325" if w.gene == "ugt" else ("Chr02_Occ:9107540-9109570" if w.gene == "cyp73" else ("Chr02_Occ:9222778-9224505" if w.gene == "cyp79" else "Chr04_Pall:28511241-28515143"))
    shell:
        """
        samtools view -hb {input} {params.region} > {output.bam} && samtools index {output.bam}
        """

rule get_num_alignments_and_mq:
    input:
        bam = rules.subset_bam.output
    output:
        f"{READ_ANALYSIS_DIR}/alignments_mq_{{sample}}_{{gene}}.txt"
    conda: "../envs/finding_genes.yaml"
    script:
        "../scripts/python/get_num_alignments_and_mq.py"

rule finding_genes_done:
    input:
        expand(rules.star_align_reads_toOlsen.output, sample=LEAF_SAMPLES),
        expand(rules.minimap_align_hcn_genes.output, ref=["olsen", "utm"]),
        expand(rules.get_num_alignments_and_mq.output, sample=SAMPLES, gene=["ugt", "cyp79", "cyp73", "li"])
    output:
        f"{MINIMAP_DIR}/minimap.done"
    shell:
        """
        touch {output}
        """
