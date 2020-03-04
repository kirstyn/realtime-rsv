
rule map_to_cns:
    input:
        cns = config["output_path"] + "/binned_{sample}/{analysis_stem}.consensus.fasta",
        reads = config["output_path"]+"/binned_{sample}/{analysis_stem}.primer_trimmed.fastq"
    output:
        paf = config["output_path"] + "/binned_{sample}/medaka/{analysis_stem}/consensus.mapped.paf"
    shell:
        "minimap2 -x map-ont {input.cns} {input.reads} > {output}"

rule mask_low_coverage_regions:
    input:
        cns = config["output_path"] + "/binned_{sample}/{analysis_stem}.consensus.fasta",
        paf = config["output_path"] + "/binned_{sample}/medaka/{analysis_stem}/consensus.mapped.paf"
    params:
        path_to_script = workflow.current_basedir,
        min_reads = config["min_reads"]
    output:
        config["output_path"] +"/binned_{sample}/medaka/{analysis_stem}/consensus_masked.fasta"
    shell:
        """
        python {params.path_to_script}/mask_low_coverage.py \
        --cns {input.cns} \
        --paf {input.paf} \
        --min_coverage 30 \
        --masked_cns {output}
        """

rule gather_files:
    input:
        expand(config["output_path"] +"/binned_{{sample}}/medaka/{analysis_stem}/consensus_masked.fasta", analysis_stem=config["analysis_stem"])
    output:
        config["output_path"] + "/consensus_sequences/{sample}.fasta"
    shell:
        "cat {input} > {output}"

