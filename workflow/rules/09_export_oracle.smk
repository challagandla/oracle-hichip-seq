# Stage 09 — Export to ORACLE Chromatin Operating System (COS) format
# Per sample, build a multi-resolution graph at 5 kb / 25 kb / 100 kb / 1 Mb
# with per-bin prototype node features, loop edges, and global tokens. Saved as
# PyTorch Geometric HeteroData (.pt) and mirrored to HDF5.

rule annotate_loops:
    """Annotate loop anchors with peak overlap and nearest genes."""
    input:
        loops = RESULTS / f"loops/{{sample}}/{{sample}}.interactions_FitHiC_{FITHICHIP_Q_LABEL}.bed",
        peaks = RESULTS / "peaks/{sample}_peaks.bed",
        gtf   = GENOME["gtf"]
    output:
        bedpe = RESULTS / "oracle_cos/intermediates/{sample}.annotated_loops.bedpe"
    threads: 2
    conda: "../envs/pyranges.yaml"
    log:
        RESULTS / "logs/annotate_loops/{sample}.log"
    script:
        "../scripts/bedpe_annotate.py"


rule export_oracle_cos:
    """
    Produce the canonical ORACLE input for a sample:
        - .pt    PyTorch Geometric HeteroData with node + edge attributes
        - .h5    HDF5 mirror for ablation/debugging
        - manifest JSON describing channels and limitations
    """
    input:
        mcool = RESULTS / "cool/{sample}.mcool",
        loops_annot = RESULTS / "oracle_cos/intermediates/{sample}.annotated_loops.bedpe",
        peaks = RESULTS / "peaks/{sample}_peaks.bed",
        insul = RESULTS / "qc/insulation/{sample}.insulation.tsv",
        eigs  = RESULTS / "qc/compartments/{sample}.cis.eigs.tsv",
        loop_qc = RESULTS / "qc/loop_qc/{sample}.json"
    output:
        pt = RESULTS / "oracle_cos/{sample}.pt",
        h5 = RESULTS / "oracle_cos/{sample}.h5",
        manifest = RESULTS / "oracle_cos/{sample}.manifest.json"
    params:
        bin_sizes_bp = ORACLE_BIN_SIZES_BP,
        chromsizes = GENOME["chromsizes"],
        # Per-assembly blacklist from genome.yaml (single source of truth).
        blacklist = GENOME.get("blacklist", ""),
        microbiome_tsv = config["oracle_export"].get("microbiome_metadata_tsv", ""),
        drop_chroms = config["oracle_export"]["drop_chromosomes"],
        emit_bigwigs = config["oracle_export"]["emit_bigwigs"]
    threads: 8
    conda: "../envs/oracle_export.yaml"
    log:
        RESULTS / "logs/export_oracle/{sample}.log"
    script:
        "../scripts/export_oracle_cos.py"
