rule getVirusReads:
    input:
        PE1 = os.path.join(DATAFOLDER["classified"], "{sample}", "{sample}.R_1.fastq.gz"),
        PE2 = os.path.join(DATAFOLDER["classified"], "{sample}", "{sample}.R_2.fastq.gz")
    output:
        PE1 = (os.path.join(DATAFOLDER["filtered"], "{sample}", "{sample}.R1.fastq.gz")),
        PE2 = (os.path.join(DATAFOLDER["filtered"], "{sample}", "{sample}.R2.fastq.gz"))
    log:
        os.path.join(DATAFOLDER["logs"], "filtered", "{sample}.extract.log")
    shell:
        r"""
            set +o pipefail
            zgrep -A3 'kraken:taxid|11308\|kraken:taxid|11320\|kraken:taxid|11520' {input.PE1} | sed -e 's/^--$//' | sed '/^\s*$/d' | gzip 1> {output.PE1} 2>> {log}
            zgrep -A3 'kraken:taxid|11308\|kraken:taxid|11320\|kraken:taxid|11520' {input.PE2} | sed -e 's/^--$//' | sed '/^\s*$/d' | gzip 1> {output.PE2} 2>> {log}
        """
