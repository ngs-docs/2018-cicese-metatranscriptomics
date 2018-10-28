from os.path import join

rule rcorrector_pe:
    """
    Run Rcorrector
    """
    input:
       r1= lambda wildcards: join(DATA_DIR, '{}_{}_1.trim.fq.gz'.format(wildcards.sample,wildcards.unit))
       r2= lambda wildcards: join(DATA_DIR, '{}_{}_2.trim.fq.gz'.format(wildcards.sample,wildcards.unit))
    output:
        r1=join(TRIM_DIR, "{sample}_{unit}_1.rcorr.fq.gz"),
        r2=join(TRIM_DIR, "{sample}_{unit}_2.rcorr.fq.gz"),
        r1_unpaired=join(TRIM_DIR, "{sample}_{unit}_1.se.rcorr.fq.gz"),
        r2_unpaired=join(TRIM_DIR, "{sample}_{unit}_2.se.rcorr.fq.gz"),
    message:
        """--- PE Rcorrector"""
    threads:4
    params:
        extra = ''
    log:
       join(LOGS_DIR, 'rcorrector/{sample}_{unit}_pe.log')
    conda: "rcorrector-env.yaml"
    script: "rcorrector.py"