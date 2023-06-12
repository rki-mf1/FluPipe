rule getReadLength:
    input:
        PE1 = (os.path.join(DATAFOLDER["trimmed"], "{sample}", "{sample}.R1.fastq.gz")),
        PE2 = (os.path.join(DATAFOLDER["trimmed"], "{sample}", "{sample}.R2.fastq.gz"))
    output:
        PE1_tmp = temp(os.path.join(DATAFOLDER["trimmed"], "{sample}", "read_length_{sample}.R1.fastq.txt")),
        PE2_tmp = temp(os.path.join(DATAFOLDER["trimmed"], "{sample}", "read_length_{sample}.R2.fastq.txt"))
    shell:
        r"""
            # 1. get length of each read per file and ouput read_length file with sample name on top
            zcat {input.PE1} | awk 'NR%4==2 {{print length}}' | sed -e '1s/^/{wildcards.sample}.R1\n/' > {output.PE1_tmp}
            zcat {input.PE2} | awk 'NR%4==2 {{print length}}' | sed -e '1s/^/{wildcards.sample}.R2\n/' > {output.PE2_tmp}
        """