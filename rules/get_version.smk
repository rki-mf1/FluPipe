#singularity: "docker://rkibioinf/general:3.6.0--53569a8"
# cd ${params.sd} command problematic in container

rule getVersion:
    output:
        os.path.join(PROJFOLDER, "pipeline.version")
    params:
        sd = srcdir(".")
    log:
        os.path.join(DATAFOLDER["logs"], "version", "git.repository.version.log")
    shell:
        r"""
            cd {params.sd};
            git describe 1> {output} 2> {log} || echo 'unknown_version' 1> {output} 2> {log}
        """