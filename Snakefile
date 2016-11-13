import os
from glob import glob

if 'bdiv_env' in config['envs']:
    bdiv_env = config['envs']['bdiv_env']
if 'qiime_env' in config['envs']:
    qiime_env = config['envs']['qiime_env']
if 'codiv_env' in config['envs']:
    codiv_env = config['envs']['codiv_env']

pOTU_chunks = config['codiv_params']['pOTU_chunks']

codiv_scripts_dir = config['envs']['codiv_scripts_dir']

localrules: annotate_bdiv_tree, annotate_cotu_trees

rule all:
    input:
        'data/bdiv_summary/master_tree_annotated.pdf',
        expand('data/codiv/{width}/{test}_{field}/uncorrected_sig_nodes.txt',
                width=config['codiv_params']['pOTU_widths'],
                test=config['codiv_params']['codiv_test'],
                field=config['codiv_params']['collapse_field']),
        expand('data/codiv/{width}/{test}_{field}/{correction}/annotate_cotu_trees.done',
                width=config['codiv_params']['pOTU_widths'],
                test=config['codiv_params']['codiv_test'],
                field=config['codiv_params']['collapse_field'],
                correction=config['codiv_params']['annotate_trees'])

rule beta_div:
    input:
        'data/bdiv_summary/master_tree_annotated.pdf'


rule codiversification:
    input:
        expand('data/codiv/{width}/{test}_{field}/uncorrected_sig_nodes.txt',
                width=config['codiv_params']['pOTU_widths'],
                test=config['codiv_params']['codiv_test'],
                field=config['codiv_params']['collapse_field']),
        expand('data/codiv/{width}/{test}_{field}/{correction}/annotate_cotu_trees.done',
                width=config['codiv_params']['pOTU_widths'],
                test=config['codiv_params']['codiv_test'],
                field=config['codiv_params']['collapse_field'],
                correction=config['codiv_params']['annotate_trees'])


### Rules for generating param files

rule make_pOTU_picking_params:
    output:
        'data/params/{width}_pOTU_params.txt'
    params:
        method = config['codiv_params']['otu_picking_method']
    run:
        with open(output[0], 'w') as out:
            out.write('pick_otus:similarity {:3.2f}\n'.format(float(wildcards.width)/100))
            out.write('pick_otus:otu_picking_method {}\n'.format(params.method))


rule make_bdiv_params:
    output:
        'data/params/bdiv_params.txt'
    run:
        with open(output[0], 'w') as out:
            out.write('beta_diversity:metrics %s\n' %
                      ','.join(config['bdiv_params']['bdiv_metrics']))
            out.write('make_emperor:color_by %s\n' %
                      ','.join(config['bdiv_params']['color_by']))
            out.write('multiple_rarefactions_even_depth:num_reps %s\n' %
                      config['bdiv_params']['jk_permutations'])


### Prep rules used for both beta div and cospeciation

rule inflate_deblur:
    input:
        config['input_files']['deblurred_biom']
    output:
        'data/starting_files/inflated_deblurred.fna'
    run:
        shell("""
              set +u; {qiime_env}; set -u

              inflate_deblur.py {input} > {output}
              """)

rule remove_chimeras:
    input:
        'data/starting_files/inflated_deblurred.fna'
    output:
        'data/starting_files/inflated_deblurred.nochimera.fna'
    params:
        chimera_dir = 'data/starting_files/chimera_checked',
        ref_file = config['input_files']['chimera_ref'],
        non_chimeras_retention = 'intersection'
    log:
        chimera = 'logs/setup/remove_chimeras-identify_chimeric_seqs.log',
        filter = 'logs/setup/remove_chimeras-filter_fasta.log'
    benchmark:
        'benchmarks/setup/remove_chimeras.json'
    run:
        shell("""
              set +u; {qiime_env}; set -u

              identify_chimeric_seqs.py -i {input} \
              -m usearch61 \
              -o {params.chimera_dir} \
              -r {params.ref_file} \
              --non_chimeras_retention {params.non_chimeras_retention} 1> {log.chimera} 2>&1

              filter_fasta.py -f {input} \
              -o {output} \
              -s {params.chimera_dir}/chimeras.txt -n 1> {log.filter} 2>&1
              """)


rule pOTU_clustering:
    input:
        fasta = 'data/starting_files/inflated_deblurred.nochimera.fna',
        params = 'data/params/{width}_pOTU_params.txt'
    output:
        biom = 'data/pOTUs/{width}/otu_table.biom',
        otu_map = 'data/pOTUs/{width}/otu_map.txt',
        tree = 'data/pOTUs/{width}/rep_set.tre'
    params:
        output_dir = 'data/pOTUs/{width}'
    threads:
        8
    log:
        'logs/setup/pOTU_clustering-{width}.log'
    benchmark:
        'benchmarks/setup/pOTU_clustering-{width}.json'
    run:

        shell("""
              set +u; {qiime_env}; set -u

              pick_de_novo_otus.py -i {input.fasta} \
              -o {params.output_dir} \
              -p {input.params} \
              -aO {threads} -f 1> {log} 2>&1

              ln -s `readlink -f {params.output_dir}/*picked_otus/*_otus.txt` {output.otu_map}
              """)


### Beta-div sensitivity rules

rule jk_beta_diversity:
    input:
        tree = 'data/pOTUs/{width}/rep_set.tre',
        biom = 'data/pOTUs/{width}/otu_table.biom',
        md_map = config['input_files']['sample_map'],
        beta_params = 'data/params/bdiv_params.txt'
    output:
        ['data/pOTUs/{width}/jk_beta_div/%s/rare_upgma_consensus.tre' %
         x for x in config['bdiv_params']['bdiv_metrics']]
    params:
        depth = config['bdiv_params']['jk_depth'],
        output_dir = 'data/pOTUs/{width}/jk_beta_div'
    threads:
        8
    log:
        'logs/betadiv/jk_beta_diversity-{width}.log'
    benchmark:
        'benchmarks/betadiv/jk_beta_diversity-{width}.json'
    run:
        shell("""
              set +u; {qiime_env}; set -u

              jackknifed_beta_diversity.py -i {input.biom} \
              -o {params.output_dir} \
              -e {params.depth} \
              -m {input.md_map} \
              -t {input.tree} \
              -p {input.beta_params} \
              -aO {threads} -f 1> {log} 2>&1
              """)


rule compare_beta_diversity:
    input:
        upgma = 'data/pOTUs/{width}/jk_beta_div/{metric}/rare_upgma_consensus.tre',
        host_tree = config['input_files']['bdiv_host_tree_fp']
    output:
        support = 'data/bdiv_summary/{metric}_{width}/jackknife_support.txt',
        tree = 'data/bdiv_summary/{metric}_{width}/master_tree.tre',
    params:
        jack_dir = 'data/pOTUs/{width}/jk_beta_div/{metric}/rare_upgma',
        output_dir = 'data/bdiv_summary/{metric}_{width}'
    log:
        'logs/betadiv/compare_beta_diversity-{width}-{metric}.log'
    benchmark:
        'benchmarks/betadiv/compare_beta_diversity-{width}-{metric}.json'
    run:
        shell("""
              set +u; {qiime_env}; set -u

              tree_compare.py -m {input.host_tree} \
              -s {params.jack_dir} \
              -o {params.output_dir} 1> {log} 2>&1
              """)


rule annotate_bdiv_tree:
    input:
        expand('data/bdiv_summary/{metric}_{width}/master_tree.tre',
               metric=config['bdiv_params']['bdiv_metrics'],
               width=config['bdiv_params']['bdiv_OTU_widths'])
    output:
        'data/bdiv_summary/master_tree_annotated.pdf'
    params:
        input_dir = 'data/bdiv_summary',
        beta_metrics = ','.join(config['bdiv_params']['bdiv_metrics']),
        otu_widths = ','.join(str(x) for x in config['bdiv_params']['bdiv_OTU_widths'])
    log:
        'logs/betadiv/annotate_bdiv_tree.log'
    benchmark:
        'benchmarks/betadiv/annotate_bdiv_tree.json'
    run:
        shell("""
              set +u; {bdiv_env}; set -u

              annotate_bdiv_tree.py \
              --beta_metrics {params.beta_metrics} \
              --otu_widths {params.otu_widths} \
              --input_dir {params.input_dir} \
              --output_fp {output} \
              --tree_fp {input[0]} 1> {log} 2>&1
              """)


### Codiversification rules

rule collapse_pOTU_tables:
    input:
        biom = 'data/pOTUs/{width}/otu_table.biom',
        md_map = config['input_files']['sample_map'],
    output:
        biom = 'data/codiv/{width}/otu_table.%s.biom' %
                config['codiv_params']['collapse_field'],
        collapse_map = 'data/codiv/{width}/map.%s.txt' %
                config['codiv_params']['collapse_field']
    params:
        collapse_field = config['codiv_params']['collapse_field']
    log:
        'logs/codiv/collapse_pOTU_tables-{width}.log'
    benchmark:
        'benchmarks/codiv/collapse_pOTU_tables-{width}.json'
    run:
        shell("""
              set +u; {qiime_env}; set -u

              collapse_samples.py -b {input.biom} \
              -m {input.md_map} \
              --output_biom_fp {output.biom} \
              --output_mapping_fp {output.collapse_map} \
              --collapse_fields {params.collapse_field} 1> {log} 2>&1
              """)


rule rarify_pOTU_tables:
    input:
        biom = 'data/codiv/{width}/otu_table.%s.biom' %
                config['codiv_params']['collapse_field']
    output:
        biom = 'data/codiv/{width}/otu_table.%s.%s.biom' %
                (config['codiv_params']['collapse_field'],
                 config['codiv_params']['pOTU_rarify'])
    params:
        pOTU_rarify = config['codiv_params']['pOTU_rarify']
    log:
        'logs/codiv/rarify_pOTU_tables-{width}.log'
    benchmark:
        'benchmarks/codiv/rarify_pOTU_tables-{width}.json'
    run:
        shell("""
              set +u; {qiime_env}; set -u

              # rarify for even sampling across species
              single_rarefaction.py -i {input.biom} \
              -o {output.biom} \
              -d {params.pOTU_rarify} 1> {log} 2>&1
              """)


rule subset_pOTU_tables:
    input:
        biom = 'data/codiv/{width}/otu_table.%s.%s.biom' %
                (config['codiv_params']['collapse_field'],
                 config['codiv_params']['pOTU_rarify'])
    output:
        biom = 'data/codiv/{width}/otu_table.%s.%s.s%s.biom' %
                (config['codiv_params']['collapse_field'],
                 config['codiv_params']['pOTU_rarify'],
                 config['codiv_params']['min_obvs'])
    params:
        min_obvs = config['codiv_params']['min_obvs']
    log:
        'logs/codiv/subset_pOTU_tables-{width}.log'
    benchmark:
        'benchmarks/codiv/subset_pOTU_tables-{width}.json'
    run:
        shell("""
              set +u; {qiime_env}; set -u

              # subset filtered otu table by abundance
              filter_otus_from_otu_table.py \
              -i {input.biom} \
              -s {params.min_obvs} \
              -o {output.biom} 1> {log} 2>&1
              """)


rule split_pOTU_tables:
    input:
        biom = 'data/codiv/{width}/otu_table.%s.%s.s%s.biom' %
                (config['codiv_params']['collapse_field'],
                 config['codiv_params']['pOTU_rarify'],
                 config['codiv_params']['min_obvs'])
    output:
        ['data/codiv/{width}/temp_biom/chunk%s.biom' %
         x for x in list(range(1, pOTU_chunks + 1))]
    params:
        output_dir = 'data/codiv/{width}/temp_biom'
    log:
        'logs/codiv/split_pOTU_tables-{width}.log'
    benchmark:
        'benchmarks/codiv/split_pOTU_tables-{width}.json'
    run:
        shell("""
              set +u; {qiime_env}; set -u

              split_biom.py -i {input.biom} \
              -o {params.output_dir} \
              -n {pOTU_chunks} 1> {log} 2>&1
              """)


rule subcluster_pOTUs:
    input:
        biom = 'data/codiv/{width}/temp_biom/chunk{chunk}.biom',
        otu_map = 'data/pOTUs/{width}/otu_map.txt',
        fasta = 'data/starting_files/inflated_deblurred.nochimera.fna',
        params = config['input_files']['cOTU_params_fp']
    output:
        'data/codiv/{width}/subclustered_otus/chunk{chunk}.done'
    log:
        'logs/codiv/subcluster_pOTUs-{width}-{chunk}.log'
    benchmark:
        'benchmarks/codiv/subcluster_pOTUs-{width}-{chunk}.json'
    run:
        shell("""
              set +u; {qiime_env}; set -u

              python {codiv_scripts_dir}/otu_subcluster.py \
              -i {input.otu_map} \
              -o data/codiv/{wildcards.width}/subclustered_otus \
              -f {input.fasta} \
              -p {input.params} \
              --force \
              -b {input.biom} 1> {log} 2>&1

              touch {output}
              """)


rule test_cospeciation:
    input:
        lambda wildcards: expand('data/codiv/{width}/subclustered_otus/chunk{chunk}.done',
                                 chunk=list(range(1, pOTU_chunks + 1)),
                                 width=wildcards.width),
        biom = 'data/codiv/{width}/otu_table.%s.%s.s%s.biom' %
                (config['codiv_params']['collapse_field'],
                 config['codiv_params']['pOTU_rarify'],
                 config['codiv_params']['min_obvs']),
        host_tree = config['input_files']['codiv_host_tree_fp'],
        md_map = config['input_files']['sample_map'],
    output:
        'data/codiv/{width}/%s_%s/uncorrected_sig_nodes.txt' %
                            (config['codiv_params']['codiv_test'],
                             config['codiv_params']['collapse_field']),
        'data/codiv/{width}/%s_%s/bh_fdr_sig_nodes.txt' %
                            (config['codiv_params']['codiv_test'],
                             config['codiv_params']['collapse_field']),
        'data/codiv/{width}/%s_%s/bonferroni_sig_nodes.txt' %
                            (config['codiv_params']['codiv_test'],
                             config['codiv_params']['collapse_field']),
        'data/codiv/{width}/%s_%s/FDR_sig_nodes.txt' %
                            (config['codiv_params']['codiv_test'],
                             config['codiv_params']['collapse_field']),
        'data/codiv/{width}/%s_%s/cospeciation_results_summary.txt' %
                            (config['codiv_params']['codiv_test'],
                             config['codiv_params']['collapse_field'])
    params:
        in_dir = 'data/codiv/{width}/subclustered_otus',
        out_dir = 'data/codiv/{width}/%s_%s' %
                              (config['codiv_params']['codiv_test'],
                               config['codiv_params']['collapse_field']),
        codiv_test = config['codiv_params']['codiv_test'],
        collapse_field = config['codiv_params']['collapse_field'],
        permutations = config['codiv_params']['permutations']
    log:
        'logs/codiv/test_cospeciation-{width}.log'
    benchmark:
        'benchmarks/codiv/test_cospeciation-{width}.json'
    run:
        shell("""
              set +u; {qiime_env}; set -u

              python {codiv_scripts_dir}/test_cospeciation.py \
              -i {params.in_dir} \
              -p {input.biom} \
              --host_tree_fp {input.host_tree} \
              -o {params.out_dir} \
              -T {params.codiv_test} \
              -m {input.md_map} \
              --collapse_fields {params.collapse_field} \
              --permutations {params.permutations} \
              --force 1> {log} 2>&1
              """)

rule annotate_cotu_trees:
    input:
        table = 'data/codiv/{width}/%s_%s/{correction}.txt' %
                            (config['codiv_params']['codiv_test'],
                             config['codiv_params']['collapse_field']),
        host_tree = config['input_files']['codiv_host_tree_fp']
    output:
        'data/codiv/{width}/%s_%s/{correction}/annotate_cotu_trees.done' %
                            (config['codiv_params']['codiv_test'],
                             config['codiv_params']['collapse_field'])
    params:
        collapse_field = config['codiv_params']['collapse_field'],
        subcluster_dir = 'data/codiv/{width}/subclustered_otus',
        output_dir = 'data/codiv/{width}/%s_%s/{correction}' %
                            (config['codiv_params']['codiv_test'],
                             config['codiv_params']['collapse_field'])
    log:
        'logs/codiv/annotate_cotu_trees-{width}-{correction}.log'
    benchmark:
        'benchmarks/codiv/annotate_cotu_trees-{width}-{correction}.json'
    run:
        shell("""
              set +u; {bdiv_env}; set -u

              annotate_cotu_trees.py
              --pies \
              --labels \
              --collapse_field {params.collapse_field} \
              --results_table_fp {input.table} \
              --subcluster_dir {params.subcluster_dir} \
              --host_tree_fp {input.host_tree} \
              --taxonomy \
              --output_dir {params.output_dir}
              --force 

              touch {output}
              """)