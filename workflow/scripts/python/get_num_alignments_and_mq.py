# Script to generate text file with alignments for each read
# Processes one BAM file at a time

import re
import pysam

bam = snakemake.input[0]
bam_file = pysam.AlignmentFile(bam, 'rb')

sample = snakemake.wildcards['sample']
gene = snakemake.wildcards['gene']

# pattern = 'utg\\d{6}l' if ref == 'stanford_full' else 'Scaffold_\\d+:\\d+-\\d+\\([+|-]\\)'
with open(snakemake.output[0], 'w') as fout:
    for align in bam_file:
        if not align.is_supplementary and align.is_mapped:
            qname = align.query_name
            mq = align.mapping_quality
            read = '1' if align.is_forward else '2'
            alignments = [align.reference_name]
            if align.has_tag('XA'):
                second_align = re.findall(pattern, align.get_tag('XA'))
                alignments = alignments + second_align
            fout.write(f"{sample}\t{gene}\t{qname}\t{read}\t{mq}\t{len(alignments)}\t{';'.join(alignments)}\n")
