rule dedupReads:
	input:
		bam = os.path.join(DATAFOLDER["mapping"], "{sample}", "{sample}.bam")
	output:
		metric = temp(os.path.join(DATAFOLDER["mapping"], "{sample}", "{sample}.dedup.txt")),
		dedupped_bam = temp(os.path.join(DATAFOLDER["mapping"], "{sample}", "{sample}.dedupped.bam")),
		sortedbam = temp(os.path.join(DATAFOLDER["mapping"], "{sample}", "{sample}.picard.tmp.bam"))
	log:
		log1=os.path.join(DATAFOLDER["logs"], "dedup", "{sample}.dedup.log"),
		log2=os.path.join(DATAFOLDER["logs"], "dedup", "{sample}.dedup2.log")
	conda:
		"../envs/picard.yaml"
	threads: 1
	shell:
		r"""
			(picard SortSam -INPUT {input.bam} -OUTPUT {output.sortedbam} -SORT_ORDER coordinate &> {log.log1}) &&
			(picard MarkDuplicates -INPUT {output.sortedbam} -OUTPUT {output.dedupped_bam} --QUIET true -METRICS_FILE {output.metric} -REMOVE_DUPLICATES true &> {log.log2})
    	"""
