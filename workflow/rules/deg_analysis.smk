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
