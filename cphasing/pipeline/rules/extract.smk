
rule extract:
    input:
        pairs="workdir/{sample}.pairs.gz",
        contigsizes=f"workdir/{{ref_prefix}}.contigsizes"
    output:
        temp("workdir/{sample}.edges")
    threads: 4
    log: "logs/{sample}.log"
    shell:
        "cphasing extract -t {threads} {input.pairs} "
            "{input.contigsizes} {output} 2>{log}"
    

    
        
