if REFERENCE:
    rule prepareReference:
        input:
            REFERENCE
        output:
            os.path.join(DATAFOLDER["mapping"], "{sample}", os.path.basename(REFERENCE))
        conda:
            "../envs/python.yaml"
        params:
            script = workflow.source_path("../scripts/dnasafe_fasta.py")
        script:
            "{params.script}"
else:
    rule prepareReference:
        input:
            os.path.join(DATAFOLDER["mapping"], "{sample}", "temp_final_refs_{sample}.fasta")
        output:
            os.path.join(DATAFOLDER["mapping"], "{sample}", "final_refs_{sample}.fasta")
        conda:
            "../envs/python.yaml"
        params:
            script = workflow.source_path("../scripts/dnasafe_fasta.py")
        script:
            "{params.script}"