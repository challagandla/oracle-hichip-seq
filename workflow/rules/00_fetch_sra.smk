# Stage 00 — Fetch FASTQs from SRA.
#
# Only runs for samples whose `srr` column is non-empty; leave it blank and drop
# your own FASTQs into fastq_dir to skip this entirely.
#
# A biological replicate sequenced as several technical replicates has its runs
# listed comma-separated, and they are concatenated here into one library. That
# is the correct unit: technical replicates share a library, so keeping them apart
# would inflate n and understate the between-donor variance the differential test
# has to beat. Biological replicates stay separate.


def _srrs_for(sample):
    value = SAMPLES.loc[sample].get("srr", "")
    if pd.isna(value) or not str(value).strip():
        return []
    return [s.strip() for s in str(value).split(",") if s.strip()]


SRA_SAMPLES = [s for s in SAMPLE_IDS if _srrs_for(s)]


rule fetch_sra:
    """Download and merge the runs of one biological replicate.

    R1/R2 line counts are compared before the files are published. fasterq-dump
    will happily emit a truncated pair on a network failure, and a truncated pair
    does not fail loudly downstream — bwa reads what it is given and pairtools
    silently loses the tail. Better to fail here.
    """
    output:
        r1 = FASTQ_DIR / "{sample}_R1.fastq.gz",
        r2 = FASTQ_DIR / "{sample}_R2.fastq.gz",
    wildcard_constraints:
        sample = "|".join(re.escape(s) for s in SRA_SAMPLES) if SRA_SAMPLES else "$^",
    params:
        srrs = lambda wc: " ".join(_srrs_for(wc.sample)),
        tmp = lambda wc: str(FASTQ_DIR.parent / f"sra_tmp/{wc.sample}"),
    threads: 6
    conda: "../envs/sra.yaml"
    log:
        RESULTS / "logs/fetch_sra/{sample}.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.r1}) {params.tmp} $(dirname {log})

        r1s=""; r2s=""
        for srr in {params.srrs}; do
            # A prior interrupted fasterq-dump may leave only one mate. Never
            # accept that half-download as a reusable cache entry.
            if [ ! -s {params.tmp}/${{srr}}_1.fastq ] || [ ! -s {params.tmp}/${{srr}}_2.fastq ]; then
                rm -f {params.tmp}/${{srr}}_1.fastq {params.tmp}/${{srr}}_2.fastq
                fasterq-dump --split-files --threads {threads} \
                    --temp {params.tmp} -O {params.tmp} "$srr" >> {log} 2>&1
            fi
            if [ ! -s {params.tmp}/${{srr}}_2.fastq ]; then
                echo "$srr is not paired-end; HiChIP requires paired reads" >> {log}
                exit 1
            fi
            r1s="$r1s {params.tmp}/${{srr}}_1.fastq"
            r2s="$r2s {params.tmp}/${{srr}}_2.fastq"
        done

        cat $r1s | pigz -p {threads} -c > {output.r1}.part
        cat $r2s | pigz -p {threads} -c > {output.r2}.part

        pigz -t {output.r1}.part
        pigz -t {output.r2}.part
        lines1=$(pigz -dc {output.r1}.part | wc -l)
        lines2=$(pigz -dc {output.r2}.part | wc -l)
        if [ $((lines1 % 4)) -ne 0 ] || [ $((lines2 % 4)) -ne 0 ]; then
            echo "FASTQ line count is not divisible by four: R1=$lines1 R2=$lines2" >> {log}
            rm -f {output.r1}.part {output.r2}.part
            exit 1
        fi
        n1=$(( lines1 / 4 ))
        n2=$(( lines2 / 4 ))
        if [ "$n1" -ne "$n2" ]; then
            echo "R1 has $n1 reads but R2 has $n2: truncated download" >> {log}
            rm -f {output.r1}.part {output.r2}.part
            exit 1
        fi
        echo "{wildcards.sample}: $n1 pairs from {params.srrs}" >> {log}

        mv {output.r1}.part {output.r1}
        mv {output.r2}.part {output.r2}
        rm -rf {params.tmp}
        """
