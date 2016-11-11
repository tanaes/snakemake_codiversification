

scripts_dir = config['envs']['scripts_dir']

rule all:
    input:
        'data/bdiv_summary/master_tree_annotated.pdf'


### Rules for generating param files

rule make_pOTU_picking_params:
    output:
        'data/params/{width}_pOTU_params.txt'
    params:
        method = config['codiv_params']['otu_picking_method']
    run:
        with open(output, 'w') as out:
            out.write('pick_otus:similarity {wildcards.width}\n'
                      'pick_otus:method {params.method}')


rule make_bdiv_params:
    output:
        'data/params/bdiv_params.txt'
    run:
        with open(output, 'w') as out:
            out.write('beta_diversity:metrics %s' %
                      ','.join(config['bdiv_params']['bdiv_metrics']))
            out.write('make_emperor:color_by %s' %
                      ','.join(config['bdiv_params']['color_by']))
            out.write('multiple_rarefactions_even_depth:num_reps %s' %
                      config['bdiv_params']['jk_permutations'])


### Prep rules used for both beta div and cospeciation

rule remove_chimeras:
    input:
        config['input_files']['original_fasta']
    output:
        'data/starting_files/reduce_fasta_fasta.nochimera.fna'
    params:
        chimera_dir = 'data/starting_files/chimera_checked',
        ref_file = config['input_files']['chimera_ref'],
        non_chimears_retention = 'intersection'
    run:
        shell("""
              set +u; {QIIME_ENV}; set -u

              identify_chimeric_seqs.py -i {input} \
              -m usearch61 \
              -o {params.chimera_dir} \
              -r {params.ref_file} \
              --non_chimeras_retention {params.non_chimeras_retention}

              filter_fasta.py -f {input} \
              -o {output} \
              -s {params.chimera_dir}/chimeras.txt -n
              """)


rule pOTU_clustering:
    input:
        fasta = 'data/starting_files/reduce_fasta_fasta.nochimera.fna',
        params = 'data/params/{width}_pOTU_params.txt'
    output:
        biom = 'data/pOTUs/{width}/otu_table.biom',
        otu_map = 'data/pOTUs/{width}/otu_map.txt',
        rep_set = 'data/pOTUs/{width}/rep_set.fna',
        tree = 'data/pOTUs/{width}/rep_set.tre'
    params:
        output_dir = 'data/pOTUs/{width}'
    threads:
        8
    run:
        shell("""
              set +u; {QIIME_ENV}; set -u

              pick_de_novo_otus.py -i {input.fasta} \
              -o {params.output_dir} \
              -p {input.params} \
              -aO {threads}
              """)


### Beta-div sensitivity rules

rule jk_beta_diversity:
    input:
        tree = 'data/pOTUs/{width}/rep_set.tre',
        biom = 'data/pOTUs/{width}/otu_table.biom',
        md_map = config['input_files']['sample_map'],
        beta_params = 'data/params/bdiv_params.txt'
    output:
        'data/pOTUs/{width}/jk_beta_div/{metric}/rare_upgma_consensus.tre'
    params:
        depth = config['bdiv_params']['jk_depth'],
        output_dir = 'data/pOTUs/{width}/jk_beta_div'
    threads:
        8
    run:
        shell("""
              set +u; {QIIME_ENV}; set -u

              jackknifed_beta_diversity.py -i {input.biom} \
              -o {params.output_dir} \
              -e {params.depth} \
              -m {input.md_map} \
              -t {input.tree} \
              -p {input.beta_params} \
              -aO {threads}
              """)


rule compare_beta_diversity:
    input:
        upgma = 'data/pOTUs/{width}/jk_beta_div/{metric}/rare_upgma_consensus.tre',
        host_tree = config['input_files']['host_tree_fp']
    output:
        support = 'data/bdiv_summary/{metric}_{width}/jackknife_support.txt',
        tree = 'data/bdiv_summary/{metric}_{width}/master_tree.tre',
    params:
        jack_dir = 'data/pOTUs/{width}/jk_beta_div_{metric}/rare_upgma',
        output_dir = 'data/bdiv_summary/{metric_width}'
    run:
        shell("""
              set +u; {QIIME_ENV}; set -u

              tree_compare.py -m {input.host_tree} \
              -s {params.jack_dir} \
              -o {params.output_dir}
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
    run:
        shell("""
              set +u; {BDIV_ENV}; set -u

              annotate_bdiv_tree.py \
              --beta_metrics {params.beta_metrics} \
              --otu_widths {params.beta_metrics} \
              --input_dir ${params.input_dir} \
              --output_fp {output} \
              --tree_fp ${input[0]}
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
    run:
        shell("""
              set +u; {QIIME_ENV}; set -u

              collapse_samples.py -b {input.biom} \
              -m {input.md_map} \
              --output_biom_fp {output.biom} \
              --output_mapping_fp {output.collapse_map} \
              --collapse_fields {params.collapse_field}
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
    run:
        shell("""
              set +u; {QIIME_ENV}; set -u

              # rarify for even sampling across species
              single_rarefaction.py -i {input.biom} \
              -o {output.biom} \
              -d {params.pOTU_rarify}
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
    run:
        shell("""
              set +u; {QIIME_ENV}; set -u

              # subset filtered otu table by abundance
              filter_otus_from_otu_table.py \
              -i {input.biom} \
              -s {params.min_obvs} \
              -o {output.biom}
              """)


rule split_pOTU_tables:
    input:
        biom = 'data/codiv/{width}/otu_table.%s.%s.s%s.biom' %
                (config['codiv_params']['collapse_field'],
                 config['codiv_params']['pOTU_rarify'],
                 config['codiv_params']['min_obvs'])
    output:
        temp(lambda wildcards: expand('data/codiv/{width}/temp_biom/chunk{chunk}.biom',
                                      chunk=range(1, config['codiv_params']['pOTU_chunks'] + 1),
                                      width=wildcards.width))
    params:
        chunks = config['codiv_params']['pOTU_chunks']
    run:
        shell("""
              set +u; {QIIME_ENV}; set -u

              split_biom.py -i {input.biom} \
              -o data/codiv/temp_biom \
              -n {params.chunks}
              """)


rule subcluster_pOTUs:
    input:
        biom = 'data/codiv/{width}/temp_biom/chunk{chunk}.biom',
        otu_map = 'data/pOTUs/{width}/otu_map.txt',
        fasta = 'data/starting_files/reduce_fasta_fasta.nochimera.fna',
        params = config['input_files']['cOTU_params_fp']
    output:
        'data/codiv/{width}/subclustered_otus/chunk{chunk}.done'
    run:
        shell("""
              set +u; {QIIME_ENV}; set -u

              python ${scripts_dir}/otu_subcluster.py \
              -i {input.otu_map} \
              -o data/codiv/{wildcards.width}/subclustered_otus \
              -f {input.fasta} \
              -p {input.params} \
              --force \
              -b {input.biom}

              touch {output}
              """)


rule test_cospeciation:
    input:
        lambda wildcards: expand('data/codiv/{width}/subclustered_otus/chunk{chunk}.done',
                                 chunk=range(1, config['codiv_params']['pOTU_chunks'] + 1),
                                 width=wildcards.width),
        biom = 'data/codiv/{width}/otu_table.%s.%s.s%s.biom' %
                (config['codiv_params']['collapse_field'],
                 config['codiv_params']['pOTU_rarify'],
                 config['codiv_params']['min_obvs']),
        host_tree = config['input_files']['host_tree_fp'],
        md_map = config['input_files']['sample_map'],
    output:
        'data/codiv/{width}/%s_%s/uncorrected_sig_nodes.txt' %
        (config['codiv_params']['codiv_test'],
         config['codiv_params']['collapse_field'])
    params:
        in_dir = 'data/codiv/{width}/subclustered_otus',
        out_dir = 'data/codiv/{width}/%s_%s' %
                  (config['codiv_params']['codiv_test'],
                   config['codiv_params']['collapse_field']),
        codiv_test = config['codiv_params']['codiv_test'],
        collapse_field = config['codiv_params']['collapse_field'],
        perms = config['codiv_params']['permutations']
    run:
        shell("""
              set +u; {QIIME_ENV}; set -u

              python ${scripts_dir}/test_cospeciation.py \
              -i {params.in_dir} \
              -p {input.biom} \
              --host_tree_fp {input.host_tree} \
              -o {params.out_dir} \
              -T {params.codiv_test} \
              -m {input.md_map} \
              --collapse_fields {params.collapse_field} \
              --permutations {params.permutations}
              --force
              """)