#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    echo "Usage: $0 SOURCE_DIR DEST_DIR [EXPECTED_BAMS] [--dry-run]" >&2
    exit 2
}

[[ $# -ge 2 && $# -le 4 ]] || usage
SRC=$1
DEST=$2
EXPECTED_BAMS=70
DRY_RUN=0
shift 2
for arg in "$@"; do
    case "$arg" in
        --dry-run) DRY_RUN=1 ;;
        *) [[ "$arg" =~ ^[0-9]+$ ]] || usage; EXPECTED_BAMS=$arg ;;
    esac
done

[[ -d "$SRC" ]] || { echo "ERROR: source directory not found: $SRC" >&2; exit 1; }
[[ "$DEST" == /* ]] || { echo "ERROR: destination must be an absolute path" >&2; exit 1; }
command -v rsync >/dev/null || { echo "ERROR: rsync is required" >&2; exit 1; }
command -v sha256sum >/dev/null || { echo "ERROR: sha256sum is required" >&2; exit 1; }

TMP=$(mktemp -d "${TMPDIR:-/tmp}/dedup-transfer.XXXXXX")
trap 'rm -rf "$TMP"' EXIT
LIST="$TMP/files.list"
SRC_HASH="$TMP/source.sha256"
DST_HASH="$TMP/destination.sha256"

mapfile -d '' BAMS < <(find "$SRC" -type f -name '*.dedup.bam' -print0 | sort -z)
[[ ${#BAMS[@]} -eq $EXPECTED_BAMS ]] || {
    echo "ERROR: expected $EXPECTED_BAMS BAMs; found ${#BAMS[@]}" >&2
    exit 1
}

: > "$LIST"
for bam in "${BAMS[@]}"; do
    rel=${bam#"$SRC"/}
    printf '%s\n' "$rel" >> "$LIST"
    if [[ -s "$bam.bai" ]]; then
        printf '%s\n' "$rel.bai" >> "$LIST"
    elif [[ -s "${bam%.bam}.bai" ]]; then
        printf '%s\n' "${rel%.bam}.bai" >> "$LIST"
    else
        echo "ERROR: no BAI index found for $rel" >&2
        exit 1
    fi
done
sort -u -o "$LIST" "$LIST"

hash_list() {
    local base=$1 output=$2
    (
        cd "$base"
        while IFS= read -r rel; do
            [[ -s "$rel" ]] || { echo "ERROR: missing or empty file: $base/$rel" >&2; exit 1; }
            sha256sum -- "$rel"
        done < "$LIST"
    ) > "$output"
}

hash_list "$SRC" "$SRC_HASH"
printf 'BAMs: %s; files including indexes: %s\n' "${#BAMS[@]}" "$(wc -l < "$LIST")"

if (( DRY_RUN )); then
    mkdir -p "$DEST"
    rsync -a --dry-run --partial --files-from="$LIST" "$SRC/" "$DEST/"
    echo "Dry run complete; no file content copied."
    exit 0
fi

mkdir -p "$DEST"
rsync -a --partial --progress --files-from="$LIST" "$SRC/" "$DEST/"
hash_list "$DEST" "$DST_HASH"
diff -u "$SRC_HASH" "$DST_HASH"
cp -- "$SRC_HASH" "$DEST/transfer_manifest.sha256"
printf 'PASS: %s BAMs, %s total files; SHA-256 verified.\n' "${#BAMS[@]}" "$(wc -l < "$LIST")"
