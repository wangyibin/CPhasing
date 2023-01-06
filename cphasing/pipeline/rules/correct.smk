
rule correct:
    input: 
        fasta = config["draft_fasta"],
        pairs = config["draft_pairs"]
    output:
        f"{config["output"]}/correct/corrected.fasta",
        f"{config["output"]}/correct/corrected.bed",
        f"{config["output"]}/correct/corrected.pairs.gz",

    threads: ncpus
    params:
        resolutions = config["correct"]["resolutions"]
    shell:
        """
        cphasing correct {input.fasta} {input.pairs} -o {output[0]} -ob {output[1]}\
                         -op {output[0]} -t {threads} -r {params.resolutions}
        """
