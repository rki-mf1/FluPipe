
rule minimap2bestRefHits:
    input:
        r1 = (os.path.join(DATAFOLDER["filtered"], "{sample}", "{sample}.R1.fastq.gz")),
        r2 = (os.path.join(DATAFOLDER["filtered"], "{sample}", "{sample}.R2.fastq.gz")),
        ref = (os.path.join(DATAFOLDER["mapping"],"5Refs", "{sample}", "{segment}.{sample}_best_refs.fasta"))
    output:
        temp(os.path.join(DATAFOLDER["mapping"],"Mapping_2","{sample}", "{segment}.{sample}.minimap.bam"))
    log:
        temp(os.path.join(DATAFOLDER["mapping"],"Logs", "mapping_2", "{segment}.{sample}.log"))
    conda:
        "../envs/minimap2.yaml"
    threads:
        50
    shell:
        r"""
            (   time \
                #-Y 	In SAM output, use soft clipping for supplementary alignments.
                #--sam-hit-only  	In SAM, dont output unmapped reads. 
                #I Load at most NUM target bases into RAM for indexing [4G].
                minimap2 -t {threads} -I10g -ax sr -Y --seed 42 --sam-hit-only {input.ref} {input.r1} {input.r2} |\
                samtools view -Sb -@ {threads} -o {output}
            ) &> {log}
        """