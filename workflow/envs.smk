# Data-independent installer workflow for every unique HiChIP rule environment.
# setup.sh uses this file only to create and smoke-test packages; it never reads
# biological data, reference files, or pipeline results.

rule all:
    input:
        ".snakemake/setup-env-checks/align.ok",
        ".snakemake/setup-env-checks/cooler.ok",
        ".snakemake/setup-env-checks/coolerpy.ok",
        ".snakemake/setup-env-checks/cooltools.ok",
        ".snakemake/setup-env-checks/coreutils.ok",
        ".snakemake/setup-env-checks/fithichip.ok",
        ".snakemake/setup-env-checks/figures.ok",
        ".snakemake/setup-env-checks/hicrep.ok",
        ".snakemake/setup-env-checks/macs3.ok",
        ".snakemake/setup-env-checks/multiqc.ok",
        ".snakemake/setup-env-checks/mustache.ok",
        ".snakemake/setup-env-checks/oracle-export.ok",
        ".snakemake/setup-env-checks/pandas.ok",
        ".snakemake/setup-env-checks/pydeseq2.ok",
        ".snakemake/setup-env-checks/pygenometracks.ok",
        ".snakemake/setup-env-checks/pyranges.ok",
        ".snakemake/setup-env-checks/qc.ok",
        ".snakemake/setup-env-checks/sra.ok",
        ".snakemake/setup-env-checks/stripenn.ok",

rule check_align_env:
    output:
        touch(".snakemake/setup-env-checks/align.ok")
    conda:
        "envs/align.yaml"
    shell:
        """
        for tool in bwa bwa-mem2 samtools bgzip pairtools pairix; do command -v "$tool" >/dev/null; done
        pairtools --version >/dev/null
        touch {output}
        """

rule check_cooler_env:
    output:
        touch(".snakemake/setup-env-checks/cooler.ok")
    conda:
        "envs/cooler.yaml"
    shell:
        """
        for tool in cooler pairix; do command -v "$tool" >/dev/null; done
        cooler --version >/dev/null
        touch {output}
        """

rule check_coolerpy_env:
    output:
        touch(".snakemake/setup-env-checks/coolerpy.ok")
    conda:
        "envs/coolerpy.yaml"
    shell:
        """
        python -c 'import cooler, matplotlib, numpy, pandas, pyBigWig'
        touch {output}
        """

rule check_cooltools_env:
    output:
        touch(".snakemake/setup-env-checks/cooltools.ok")
    conda:
        "envs/cooltools.yaml"
    shell:
        """
        command -v cooltools >/dev/null
        python -c 'import cooler, cooltools, pandas'
        touch {output}
        """

rule check_coreutils_env:
    output:
        touch(".snakemake/setup-env-checks/coreutils.ok")
    conda:
        "envs/coreutils.yaml"
    shell:
        """
        for tool in awk cat gzip sort curl wget; do command -v "$tool" >/dev/null; done
        touch {output}
        """

rule check_fithichip_env:
    output:
        touch(".snakemake/setup-env-checks/fithichip.ok")
    conda:
        "envs/fithichip.yaml"
    shell:
        """
        for tool in Rscript bedtools samtools cooler curl bedToBigBed; do command -v "$tool" >/dev/null; done
        Rscript -e 'pkgs <- c("optparse","ggplot2","data.table","fdrtool","GenomicRanges","edgeR"); stopifnot(all(vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)))'
        touch {output}
        """

rule check_figures_env:
    output:
        touch(".snakemake/setup-env-checks/figures.ok")
    conda:
        "envs/figures.yaml"
    shell:
        """
        python -c 'import matplotlib, numpy, pandas'
        touch {output}
        """

rule check_hicrep_env:
    output:
        touch(".snakemake/setup-env-checks/hicrep.ok")
    conda:
        "envs/hicrep.yaml"
    shell:
        """
        python -c 'from hicrep import hicrepSCC; from hicrep.utils import readMcool'
        touch {output}
        """

rule check_macs3_env:
    output:
        touch(".snakemake/setup-env-checks/macs3.ok")
    conda:
        "envs/macs3.yaml"
    shell:
        """
        for tool in macs3 bedtools samtools; do command -v "$tool" >/dev/null; done
        macs3 --version >/dev/null
        touch {output}
        """

rule check_multiqc_env:
    output:
        touch(".snakemake/setup-env-checks/multiqc.ok")
    conda:
        "envs/multiqc.yaml"
    shell:
        """
        command -v multiqc >/dev/null
        multiqc --version >/dev/null
        touch {output}
        """

rule check_mustache_env:
    output:
        touch(".snakemake/setup-env-checks/mustache.ok")
    conda:
        "envs/mustache.yaml"
    shell:
        """
        command -v mustache >/dev/null
        touch {output}
        """

rule check_oracle_export_env:
    output:
        touch(".snakemake/setup-env-checks/oracle-export.ok")
    conda:
        "envs/oracle_export.yaml"
    shell:
        """
        python -c 'import cooler, h5py, numpy, pandas, torch; from torch_geometric.data import HeteroData'
        touch {output}
        """

rule check_pandas_env:
    output:
        touch(".snakemake/setup-env-checks/pandas.ok")
    conda:
        "envs/pandas.yaml"
    shell:
        """
        python -c 'import numpy, pandas'
        touch {output}
        """

rule check_pydeseq2_env:
    output:
        touch(".snakemake/setup-env-checks/pydeseq2.ok")
    conda:
        "envs/pydeseq2.yaml"
    shell:
        """
        python -c 'import matplotlib, numpy, pandas, scipy; from pydeseq2.dds import DeseqDataSet'
        touch {output}
        """

rule check_pygenometracks_env:
    output:
        touch(".snakemake/setup-env-checks/pygenometracks.ok")
    conda:
        "envs/pygenometracks.yaml"
    shell:
        """
        command -v pyGenomeTracks >/dev/null
        python -c 'import pandas'
        touch {output}
        """

rule check_pyranges_env:
    output:
        touch(".snakemake/setup-env-checks/pyranges.ok")
    conda:
        "envs/pyranges.yaml"
    shell:
        """
        python -c 'import pandas, pyranges'
        touch {output}
        """

rule check_qc_env:
    output:
        touch(".snakemake/setup-env-checks/qc.ok")
    conda:
        "envs/qc.yaml"
    shell:
        """
        for tool in fastqc fastp; do command -v "$tool" >/dev/null; done
        fastqc --version >/dev/null
        fastp --version >/dev/null 2>&1
        touch {output}
        """

rule check_sra_env:
    output:
        touch(".snakemake/setup-env-checks/sra.ok")
    conda:
        "envs/sra.yaml"
    shell:
        """
        for tool in fasterq-dump pigz; do command -v "$tool" >/dev/null; done
        fasterq-dump --version >/dev/null
        touch {output}
        """

rule check_stripenn_env:
    output:
        touch(".snakemake/setup-env-checks/stripenn.ok")
    conda:
        "envs/stripenn.yaml"
    shell:
        """
        command -v stripenn >/dev/null
        python -c 'import cooler, numpy, pandas, stripenn'
        touch {output}
        """
