configfile: "config.yaml"

SAMPLES, = glob_wildcards('reads/{sample}_R1.fastq.gz')

rule all:
    #input: expand('alignments/{sample}.bam.bai', sample=SAMPLES)
    input: "variant-calls/raw.vcf.gz"

rule bwa_map:
    input:
        r1="reads/{sample}_R1.fastq.gz",
        r2="reads/{sample}_R2.fastq.gz"
    output:
        "alignments/{sample}.bam"
    params:
        rg="@RG\tID:{sample}\tSM:{sample}"
    threads: 8
    resources:
        mem=lambda wildcards, attempt: 12*attempt,
        hours=lambda wildcards, attempt: 48*attempt
    shell:
        """
        module load bwa biodata samtools/1.9

        bwa mem -t {threads} -R '{params.rg}' {config[bwa_ref]} \
            {input.r1} {input.r2} | samtools view -bh - | samtools sort - \
            | samtools rmdup - {output}
        samtools index {output}
        """

# TODO make variable-sized regions based on bam coverage
rule generate_regions:
    output: "variant-calls/regions.list"
    shell:
        """
        module load biodata freebayes/1.2

        fasta_generate_regions.py {config[faidx_ref]} 1000000 \
            > {output}
        mkdir -p variant-calls/regions
        """

rule make_regions_files:
    input: "variant-calls/regions.list"
    output: dynamic("variant-calls/regions/{region}.reg")
    run:
        for i, line in enumerate(open(str(input))):
            with open("variant-calls/regions/{}.reg".format(i), 'w') as out:
                out.write(line)

rule freebayes:
    input:
        bams=expand("alignments/{sample}.bam", sample=SAMPLES),
        region="variant-calls/regions/{region}.reg"
    output: "variant-calls/regions/{region}.vcf.gz"
    resources:
        mem=lambda wildcards, attempt: 12*attempt,
        hours=lambda wildcards, attempt: 24*attempt
    shell:
        """
        module load freebayes/1.2 biodata

        ref={config[fasta_ref]}
        region=$(cat {input.region})
        freebayes --region $region -f $ref {input.bams} | gzip > {output}
        """

rule freebayes_join:
    input: dynamic("variant-calls/regions/{region}.vcf.gz")
    output: "variant-calls/raw.vcf.gz"
    resources:
        hours=80
    shell:
        """
        module load freebayes/1.2

        zcat variant-calls/regions/*.vcf.gz | vcffirstheader | vcfstreamsort \
            | vcfuniq | gzip > {output}
        """
