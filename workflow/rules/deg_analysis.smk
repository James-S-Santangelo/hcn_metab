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
        star_build = rules.star_genome_build.output,
        R1 = lambda w: glob.glob(f"{config['raw_reads']}/*.{w.sample}_R1.fastq.gz"),
        R2 = lambda w: glob.glob(f"{config['raw_reads']}/*.{w.sample}_R2.fastq.gz") 
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

#####################
#### COUNT TABLE ####
#####################

rule generate_saf_file:
    input:
        gff = config["gff"]
    output:
        f"{DEG_DIR}/feature_counts/genes.saf"
    run:
        with open(output[0], "w") as fout:
            with open(input[0], "r") as fin:
                fout.write("GeneID\tChr\tStart\tEnd\tStrand\n")
                lines = fin.readlines()
                for l in lines:
                    if not l.startswith("#"):
                        sl = l.strip().split("\t")
                        feat = sl[2]
                        if feat == "gene":
                            chr = sl[0]
                            start = sl[3]
                            end = sl[4]
                            strand = sl[6]
                            id = re.search("ACLI19_g\d+", sl[8])[0]
                            fout.write(f"{id}\t{chr}\t{start}\t{end}\t{strand}\n")

rule feature_counts:
    input:
        bams = expand(rules.star_align_reads.output, sample=SAMPLES),
        saf = rules.generate_saf_file.output
    output:
        f"{DEG_DIR}/feature_counts/feature_counts.txt"
    log: f"{LOG_DIR}/feature_counts/feature_counts.log"
    conda: "../envs/deg_analysis.yaml"
    threads: 6
    shell:
        """
        featureCounts -T {threads} \
            -a {input.saf} \
            -t gene \
            -F SAF \
            -o {output} \
            -p \
            {input.bams} 2> {log}
        """
        
