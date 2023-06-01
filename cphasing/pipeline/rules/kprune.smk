
rule alleles:
    input:
        fasta = FASTA
    output:
        alleletable = f'workdir/{FASTA_PREFIX}.allele.table'
    log: f'logs/{FASTA_PREFIX}.allele.log'
    shell:
        """
        cphasing alleles -f {input.fasta} -o {output.alleletable} 2> {log}
        """

rule kprune:
    input:
        alleletable = f'{FASTA_PREFIX}.allele.table',
        cool = 'workdir/{SAMPLE_PREFIX}.wholecool'
    output:
        prunetable = 'workdir/prune.contig.table'
    log: f'logs/{FASTA_PREFIX}.prune.log'
    shell:
        """
        cphasing kprune {input.alleletable} {input.cool} -o {output.prunetable} 2> {log}
        """
