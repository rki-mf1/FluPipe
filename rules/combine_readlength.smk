
rule combine_readlength:
    input:
        PE1 = expand(os.path.join(DATAFOLDER["trimmed"], "{sample}", "read_length_{sample}.R1.fastq.txt"), sample=SAMPLES),
        PE2 = expand(os.path.join(DATAFOLDER["trimmed"], "{sample}", "read_length_{sample}.R2.fastq.txt"), sample=SAMPLES)
    output:
        os.path.join(DATAFOLDER["trimmed"], "merged_length.tsv")
    shell:
        """
            paste {input.PE1} {input.PE2} > {output}
        """
