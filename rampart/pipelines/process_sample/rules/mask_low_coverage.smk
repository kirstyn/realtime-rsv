rule map_to_cns:
    input:
        cns = config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/medaka.clean.fasta",
        reads = config["output_path"]+"/binned_{sample}/{analysis_stem}.primer_trimmed.fastq"
    output:
        config["output_path"] + "/binned_{sample}/medaka/{analysis_stem}/consensus.mapped.paf"
    shell:
        "minimap2 -x map-ont {input.cns} {input.reads} > {output}"

rule mask_low_coverage_regions:
    input:
        cns = config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/medaka.clean.fasta",
        paf = rules.map_to_cns.output
    params:
        path_to_script = workflow.current_basedir
    output:
        config["output_path"] +"/binned_{sample}/consensus_sequences/{analysis_stem}.masked.fasta"
    shell:
        """
        python {params.path_to_script}/mask_low_coverage.py \
        --cns {input.cns} \
        --paf {input.paf} \
        --min_coverage 30 \
        --masked_cns {output}
        """

rule cat_stems:
    input:
        expand(config["output_path"] +"/binned_{{sample}}/consensus_sequences/{analysis_stem}.masked.fasta",analysis_stem=config["analysis_stem"])
    output:
        config["output_path"] + "/consensus_sequences/{sample}.split.fasta"
    shell:
        "cat {input} > {output}"

rule aln_to_merge:
    input:
        config["output_path"] + "/consensus_sequences/{sample}.split.fasta"
    output:
        config["output_path"] + "/consensus_sequences/{sample}.split.aln.fasta"
    shell:
        "mafft {input} > {output}"

rule merge:
    input:
        rules.aln_to_merge.output
    params:
        path_to_script = workflow.current_basedir
    output:
        config["output_path"] + "/consensus_sequences/{sample}.fasta"
    shell:
        """
        python {params.path_to_script}/merge.py \
        -i {input} \
        -o {output} 
        """