import os

include:
    os.path.join('configs', config['config'] + '.py')

workdir: OUT_DIR

KMER_LENGTHS = range(20, 51, 5)

rule all:
    input:
        GEM_INDEX,
        expand('{prefix}_{kmer}.mappability', prefix=GEM_INDEX_PREFIX, kmer=KMER_LENGTHS),
        expand('bigWigs/{prefix}_{kmer}.bw', prefix=GEM_INDEX_PREFIX, kmer=KMER_LENGTHS),
        'summary_counts.tsv'


rule create_gem_index:
    input: GENOME_FA
    output: GEM_INDEX
    params:
        prefix = GEM_INDEX_PREFIX
    threads: 23
    shell:
        r'''gem-indexer -T {threads} -c dna -i {input} -o {params.prefix}'''


rule calc_gem_mappability:
    input: GEM_INDEX
    output: '{prefix}_{kmer}.mappability'
    params:
        kmer='{kmer}',
        prefix='{prefix}_{kmer}'
    threads: 23
    shell:
        r'''gem-mappability -T {threads} -I {input} -l {params.kmer} -o {params.prefix}'''


rule gem_2_wig:
    input:
        gem  = GEM_INDEX,
        mappability = '{prefix}_{kmer}.mappability'
    output: 'wigs/{prefix}_{kmer}.wig'
    params:
        prefix = 'wigs/{prefix}_{kmer}'
    shell:
        r'''gem-2-wig -I {input.gem} -i {input.mappability} -o {params.prefix}'''


rule wig_2_bw:
    input: 'wigs/{prefix}_{kmer}.wig'
    output: 'bigWigs/{prefix}_{kmer}.bw'
    shell:
        r'''wigToBigWig {input} {CHROM_SIZES} {output}'''


rule calculate_uniq_mappable:
    input: '{prefix}_{kmer}.mappability'
    output: '{prefix}_{kmer}.uniq.counts'
    shell:
        r'''gawk '{{c+=gsub(s,s)}}END{{print c}}' s='!' {input} > {output}'''


rule create_summary:
    input: expand('{prefix}_{kmer}.uniq.counts', prefix=GEM_INDEX_PREFIX, kmer=KMER_LENGTHS)
    output: 'summary_counts.tsv'
    run:
        summary = ''
        for inp in input:
            counts = open(inp).read().strip()
            summary += '{}\t{}\n'.format(inp.replace('.uniq.counts', ''), counts)
        with open(str(output), 'w') as f:
            f.write(summary)
