configfile: "config.yaml"

SAMPLES, = glob_wildcards('reads/{sample}_R1.fastq.gz')
REGIONS, = glob_wildcards('variant-calls/regions/{region}.reg')

rule all:
    #input: expand('alignments/{sample}.bam.bai', sample=SAMPLES)
    input: "variant-calls/variants.Q20.sorted.vcf.gz"

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
        hours=lambda wildcards, attempt: 76*attempt
    shell:
        """
        module load bwa biodata samtools/1.9

        bwa mem -t {threads} -R '{params.rg}' {config[bwa_ref]} \
            {input.r1} {input.r2} | samtools view -bh - | samtools sort - \
            | samtools rmdup - {output}
        samtools index {output}
        """

rule freebayes:
    input:
        bams=expand("alignments/{sample}.bam", sample=SAMPLES),
        region="variant-calls/regions/{region}.reg"
    output: "variant-calls/regions/{region}.vcf.gz"
    resources:
        mem=lambda wildcards, attempt: 16*attempt,
        hours=lambda wildcards, attempt: 72*attempt
    shell:
        """
        module load freebayes/1.2 biodata

        ref={config[fasta_ref]}
        region=$(cat {input.region})
        freebayes --region $region -f $ref {input.bams} | gzip > {output}
        """

rule freebayes_join:
    input: expand("variant-calls/regions/{region}.vcf.gz", region=REGIONS)
    output: "variant-calls/raw.vcf.gz"
    resources:
        hours=80
    shell:
        """
        module load freebayes/1.2

        zcat variant-calls/regions/*.vcf.gz | vcffirstheader | vcfstreamsort \
            | vcfuniq | gzip > {output}
        """

rule quality_filter_and_sort:
    input: 'variant-calls/raw.vcf.gz'
    output: 'variant-calls/variants.Q20.sorted.vcf.gz'
    resources: mem=16, hours=10
    shell:
        """
        module load bcftools/1.8

        bcftools filter -i 'QUAL>=20' -Ou {input} | \
            bcftools sort -Oz -o {output} -m 13G
        """
