import os
from collections import defaultdict
import errno
import glob
import numpy as np
import pandas as pd
from os.path import join

include:
   config['config_path']

workdir: OUT_DIR
GEM_BINARY_DIR = '/home/cmb-06/as/skchoudh/software_frozen/GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin'
THREADS = 16
KMER_LENGTHS = range(25, 32, 1)
def mkdir_p(path):
  """Python version mkdir -p

  Parameters
  ----------

  path : str
  """
  if path:
    try:
      os.makedirs(path) 
    except OSError as exc:  # Python >2.5
      if exc.errno == errno.EEXIST and os.path.isdir(path):
        pass
      else:
        raise



mkdir_p(os.path.join(OUT_DIR, 'slurm-logs'))

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
    threads: THREADS
    shell:
        r'''{GEM_BINARY_DIR}/gem-indexer -T {threads} -c dna -i {input} -o {params.prefix}'''


"""
--mapping-mode={'fast'|'sensitive'|'customed'}
      Specify the mapping alignment approach used by the tool which ultimately 
      leads to a different accuracy/performance trade-off 
      [default=fast]

    -e, --alignment-max-error={FLOAT|INTEGER} 
      Maximum divergency rate allowed between the input sequence and the reported
      matches (i.e. maximum number of mismatches or indels allowed). All 
      global-matches above this error rate will be discarded 
      [default=0.12, 12%]

    --alignment-global-min-identity={FLOAT|INTEGER}
      Minimum global-alignment identity required (i.e. total  number of matching 
      bases). All global-matches below this identity will be discared 
      [default=80%]

    --alignment-global-min-score={FLOAT|INTEGER} 
      Minimum global-alignment score required (i.e. gap-affine score). All 
      global-matches with score below this threshold will be discarded 
      [default=0.20, 0.20*read_length*match_score]

    --alignment-local={'if-unmapped'|'never'} 
      Select whether the mapping algorithm should search for local alignments in
      case no global alignment is found, or never resort to local alignment search
      [default=if-unmapped]

    --alignment-local-min-identity={FLOAT|INTEGER} 
      Minimum local-alignment identity required (i.e. total number of matching 
      bases). All local-matches below this identity will be discarded
      [default=40]

    --alignment-local-min-score={FLOAT|INTEGER}  [default=20]
      Minimum global-alignment score required (i.e. gap-affine score). All 
      global-matches with score below this threshold will be discarded 
      [default=20]
"""

rule calc_gem_mappability:
    input: GEM_INDEX
    output: '{prefix}_{kmer}.mappability'
    params:
        kmer='{kmer}',
        prefix='{prefix}_{kmer}',
        min_match=18,
        max_mismatch=2,
    threads: THREADS
    shell:
        r'''{GEM_BINARY_DIR}/gem-mappability -T {threads} -m {params.max_mismatch} -I {input} -l {params.kmer} -o {params.prefix}'''


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
