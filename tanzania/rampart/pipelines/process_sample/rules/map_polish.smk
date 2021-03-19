rule files:
    params:
        ref=config["output_path"] + "/binned_{sample}/{analysis_stem}.fasta",
        reads=config["output_path"]+"/binned_{sample}/{analysis_stem}.primer_trimmed.fastq"

rule minimap2_racon0:
    input:
        reads=rules.files.params.reads,
        ref=rules.files.params.ref
    output:
        config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/mapped.paf"
    shell:
        "minimap2 -x map-ont {input.ref} {input.reads} > {output}"

rule racon1:
    input:
        reads=rules.files.params.reads,
        fasta=rules.files.params.ref,
        paf= rules.minimap2_racon0.output
    output:
        config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/racon1.fasta"
    shell:
        "racon -m 8 -x -6 -g -8  --no-trimming -t 1 {input.reads} {input.paf} {input.fasta} > {output}"

rule mafft1:
    input:
       fasta = rules.racon1.output,
       ref = rules.files.params.ref
    params:
        temp_file = config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/temp.racon1.fasta"
    output:
        config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/racon1.aln.fasta"
    shell:
        "cat {input.ref} {input.fasta} > {params.temp_file} && "
        "mafft {params.temp_file} > {output} && "
        "rm {params.temp_file}"

rule clean1:
    input:
        aln = rules.mafft1.output,
        cns = rules.racon1.output
    params:
        path_to_script = workflow.current_basedir,
        seq_name = "{analysis_stem}"
    output:
        config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/racon1.clean.fasta"
    shell:
        "python {params.path_to_script}/clean.py "
        "--alignment_with_ref {input.aln} "
        "--name {params.seq_name} "
        "--output_seq {output} "
        "--polish_round 1"

rule medaka:
    input:
        basecalls=rules.files.params.reads,
        draft= rules.clean1.output
    params:
        outdir=config["output_path"] + "/binned_{sample}/medaka/{analysis_stem}"
    output:
        config["output_path"] + "/binned_{sample}/medaka/{analysis_stem}/consensus.fasta"
    threads:
        2
    shell:
        "medaka_consensus -i {input.basecalls} -d {input.draft} -o {params.outdir} -f -t 2 || touch {output}"

rule join_contigs_if_split:
    input:
        config["output_path"] + "/binned_{sample}/medaka/{analysis_stem}/consensus.fasta"
    output:
        config["output_path"] + "/binned_{sample}/medaka/{analysis_stem}/consensus_combined.fasta"
    run:
        joined_seq = ""
        start, end = 0,0
        last_end = 0

        with open(output[0],"w") as fw:

            for record in SeqIO.parse(input[0],"fasta"):

                name = record.description.split(":")[0]
                start,end = [float(i) for i in record.description.split(":")[1].split("-")]

                if not joined_seq:
                    joined_seq +="N"*int(start)
                    last_end = end
                else:
                    gap = int(start-last_end)
                    joined_seq +="N"*int(gap)
                    last_end = end
                joined_seq += record.seq
                
            fw.write(f">{name}\n{joined_seq}\n") 


rule mafft5:
    input:
       fasta = config["output_path"] + "/binned_{sample}/medaka/{analysis_stem}/consensus_combined.fasta",
       ref = rules.files.params.ref
    params:
        temp_file = config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/temp.medaka.fasta"
    output:
        config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/medaka.aln.fasta"
    shell:
        "cat {input.ref} {input.fasta} > {params.temp_file} && "
        "mafft {params.temp_file} > {output} && "
        "rm {params.temp_file}"

rule clean5:
    input:
        aln = rules.mafft5.output,
        cns = config["output_path"] + "/binned_{sample}/medaka/{analysis_stem}/consensus_combined.fasta"
    params:
        path_to_script = workflow.current_basedir,
        seq_name = "{analysis_stem}"
    output:
        config["output_path"] + "/binned_{sample}/{analysis_stem}.consensus.fasta"
    shell:
        "python {params.path_to_script}/clean.py "
        "--alignment_with_ref {input.aln} "
        "--name {params.seq_name} "
        "--output_seq {output} "
        "--polish_round medaka"

