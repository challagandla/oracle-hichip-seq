#!/usr/bin/env bash
# Download reference assets used by config/genome.yaml.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RESOURCES_DIR="${ROOT_DIR}/resources"

usage() {
    cat <<'EOF'
Usage:
  bash prepare_references.sh hg38 [mm10 ...]

Supported assemblies:
  hg38  Human GRCh38 primary assembly, GENCODE v46, ENCODE blacklist
  mm10  Mouse mm10/GRCm38, GENCODE vM25, ENCODE blacklist
EOF
}

download() {
    local url="$1"
    local dest="$2"
    local checksum="$3"
    local algorithm="$4"
    local actual part
    mkdir -p "$(dirname "${dest}")"
    if [[ -s "${dest}" ]]; then
        actual="$(${algorithm}sum "${dest}" | awk '{print $1}')"
        if [[ "$actual" == "$checksum" ]]; then
            echo "verified: ${dest}"
            return
        fi
        echo "invalid checksum; replacing: ${dest}" >&2
        rm -f "${dest}"
    fi
    echo "download: ${url}"
    part="${dest}.part"
    rm -f "$part"
    if command -v curl >/dev/null 2>&1; then
        curl -L --fail --retry 3 -o "$part" "${url}"
    elif command -v wget >/dev/null 2>&1; then
        wget -O "$part" "${url}"
    else
        echo "ERROR: curl or wget is required." >&2
        exit 1
    fi
    actual="$(${algorithm}sum "$part" | awk '{print $1}')"
    [[ "$actual" == "$checksum" ]] || {
        rm -f "$part"
        echo "ERROR: checksum failed for ${url}" >&2
        exit 1
    }
    [[ "$dest" != *.gz ]] || gzip -t "$part"
    mv "$part" "$dest"
}

decompress_fasta() {
    local gz="$1"
    local fasta="$2"
    if [[ -s "${fasta}" ]]; then
        echo "exists: ${fasta}"
        return
    fi
    echo "decompress: ${gz}"
    gzip -dc "${gz}" > "${fasta}.part"
    [[ -s "${fasta}.part" ]] || { rm -f "${fasta}.part"; exit 1; }
    mv "${fasta}.part" "${fasta}"
}

index_reference() {
    local fasta="$1"
    local bwa_prefix="$2"
    local bwamem2_prefix="$3"

    # Say which indexer is missing rather than skipping in silence. A skipped index
    # does not fail here — it fails hours later, inside the first alignment job,
    # as a missing-input error that points at the index rather than at this script.
    if ! command -v samtools >/dev/null 2>&1; then
        echo "ERROR: samtools not on PATH; cannot faidx ${fasta}." >&2
        exit 1
    fi
    if ! command -v bwa >/dev/null 2>&1 && ! command -v bwa-mem2 >/dev/null 2>&1; then
        echo "ERROR: neither bwa nor bwa-mem2 on PATH; no aligner index can be built." >&2
        echo "       Install one and rerun, or the alignment stage will fail." >&2
        exit 1
    fi

    if [[ ! -s "${fasta}.fai" ]]; then
        samtools faidx "${fasta}"
    fi

    mkdir -p "$(dirname "${bwa_prefix}")" "$(dirname "${bwamem2_prefix}")"
    if command -v bwa >/dev/null 2>&1 && [[ ! -s "${bwa_prefix}.bwt" ]]; then
        bwa index -p "${bwa_prefix}" "${fasta}"
    else
        echo "note: skipping bwa index (bwa not on PATH or index present)."
    fi
    if command -v bwa-mem2 >/dev/null 2>&1 && [[ ! -s "${bwamem2_prefix}.0123" ]]; then
        bwa-mem2 index -p "${bwamem2_prefix}" "${fasta}"
    else
        echo "note: skipping bwa-mem2 index (bwa-mem2 not on PATH or index present)."
    fi
}

write_chromsizes() {
    local fasta="$1"
    local chromsizes="$2"
    if [[ -s "${chromsizes}" ]]; then
        echo "exists: ${chromsizes}"
        return
    fi
    if [[ ! -s "${fasta}.fai" ]]; then
        echo "ERROR: missing ${fasta}.fai; install samtools and rerun." >&2
        exit 1
    fi
    cut -f1,2 "${fasta}.fai" > "${chromsizes}"
}

write_digest() {
    local chromsizes="$1"
    local fasta="$2"
    local digest="$3"
    local tmp="${digest%.gz}"
    if [[ -s "${digest}" ]]; then
        echo "exists: ${digest}"
        return
    fi
    if command -v cooler >/dev/null 2>&1; then
        cooler digest -o "${tmp}" "${chromsizes}" "${fasta}" MboI
        gzip -f "${tmp}"
    else
        echo "skip digest: install cooler and rerun to create ${digest}" >&2
    fi
}

download_hg38() {
    local dir="${RESOURCES_DIR}/hg38"
    local fasta="${dir}/GRCh38.primary_assembly.genome.fa"
    download "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz" "${fasta}.gz" "a445fcadf36bcbf1cdd7839ff2c8bf95" md5
    download "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.primary_assembly.annotation.gtf.gz" "${dir}/gencode.v46.primary_assembly.annotation.gtf.gz" "b4dd7c18c24c28d083a9418cd001dcfe" md5
    download "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/61a04d2c5e49341d76735d485c61f0d1177d08a8/lists/hg38-blacklist.v2.bed.gz" "${dir}/hg38-blacklist.v2.bed.gz" "c92e763af17271446194991e71917ac220593a5a3d40a06667be24178ef08cf2" sha256
    decompress_fasta "${fasta}.gz" "${fasta}"
    index_reference "${fasta}" "${dir}/bwa_index/GRCh38.primary_assembly.genome.fa" "${dir}/bwamem2_index/GRCh38.primary_assembly.genome.fa"
    write_chromsizes "${fasta}" "${dir}/hg38.chrom.sizes"
    write_digest "${dir}/hg38.chrom.sizes" "${fasta}" "${dir}/MboI.digest.hg38.bed.gz"
}

download_mm10() {
    local dir="${RESOURCES_DIR}/mm10"
    local fasta="${dir}/mm10.fa"
    download "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz" "${fasta}.gz" "db005b65828db31735f384e4c5787be5" md5
    download "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz" "${dir}/gencode.vM25.annotation.gtf.gz" "0c38fc4ccbc731a2708fc91e7f1c2efd" md5
    download "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/61a04d2c5e49341d76735d485c61f0d1177d08a8/lists/mm10-blacklist.v2.bed.gz" "${dir}/mm10-blacklist.v2.bed.gz" "febafb843c6df492f9a9fc418f8796762ee899d9864330fb509ae2d38ddc0b46" sha256
    decompress_fasta "${fasta}.gz" "${fasta}"
    index_reference "${fasta}" "${dir}/bwa_index/mm10.fa" "${dir}/bwamem2_index/mm10.fa"
    write_chromsizes "${fasta}" "${dir}/mm10.chrom.sizes"
    write_digest "${dir}/mm10.chrom.sizes" "${fasta}" "${dir}/MboI.digest.mm10.bed.gz"
}

if [[ "$#" -eq 0 ]]; then
    usage
    exit 1
fi

for assembly in "$@"; do
    case "${assembly}" in
        hg38) download_hg38 ;;
        mm10) download_mm10 ;;
        -h|--help) usage; exit 0 ;;
        *)
            echo "ERROR: unsupported assembly '${assembly}'." >&2
            usage >&2
            exit 1
            ;;
    esac
done

cat <<EOF

Reference download complete.

Paths are written under:
  ${RESOURCES_DIR}

These paths match config/genome.yaml for the supported assemblies.
EOF
