def input_sortBam(wildcards):
	if PRIMER:
		return os.path.join(DATAFOLDER["mapping"], "{sample}", "{sample}.primerclipped.bam")
	if PCR_DEDUP:
		return os.path.join(DATAFOLDER["dedup"], "{sample}", "{sample}.dedupped.bam")
	return os.path.join(DATAFOLDER["mapping"], "{sample}", "{sample}.bam")

rule sortBam:
    input:
        #input_sortBam
        os.path.join(DATAFOLDER["mapping"], "{sample}", "{sample}.dedupped.bam")
    output:
        temp(os.path.join(DATAFOLDER["mapping"], "{sample}", "{sample}.sort.bam"))
    log:
        os.path.join(DATAFOLDER["logs"], "mapping", "{sample}.sort.log")
    conda:
        "../envs/samtools.yaml"
    threads:
        10
    shell:
        r"""
            samtools sort -@ {threads} -o {output} {input} &> {log}
        """
