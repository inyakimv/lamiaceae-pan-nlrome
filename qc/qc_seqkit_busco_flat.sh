#!/bin/bash
# ==========================================================
# QC combinado: SeqKit + BUSCO
# Versión adaptada para carpeta PLANA (todos los .fna.gz
# directamente en el directorio actual, sin subcarpetas)
#
# Requiere (instalar una sola vez):
#   conda install -c bioconda seqkit
#   conda create -n busco_env -c bioconda -c conda-forge busco=5.7.1
#   conda activate busco_env
#   busco --download eudicots_odb10
#   conda deactivate
#
# Uso:
#   Ponlo en la misma carpeta donde están los *_genomic.fna.gz
#   bash qc_seqkit_busco_flat.sh
# ==========================================================

set -uo pipefail

LINEAGE="eudicots_odb10"      # cambia si alguna especie no es eudicotiledonea
THREADS=32                     # ajustado para servidor compartido de 256 nucleos
BUSCO_ENV="busco_env"

GENOMES_DIR="00_genomes"
RESULTS_DIR="01_qc"
mkdir -p "$RESULTS_DIR"

SUMMARY="$RESULTS_DIR/qc_summary.csv"
echo "assembly_file,num_seqs,sum_len_bp,N50,GC_percent,busco_complete_pct,busco_single_pct,busco_duplicated_pct,busco_fragmented_pct,busco_missing_pct" > "$SUMMARY"

CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"

shopt -s nullglob
FILES=("$GENOMES_DIR"/*_genomic.fna.gz)

if [ ${#FILES[@]} -eq 0 ]; then
    echo "ERROR: no encuentro archivos *_genomic.fna.gz en esta carpeta."
    exit 1
fi

echo "Se encontraron ${#FILES[@]} genomas."

for FNA_GZ in "${FILES[@]}"; do
    FILENAME=$(basename "$FNA_GZ")
    BASENAME="${FILENAME%.fna.gz}"   # ej: GCA_051938595.1_PrHapA_genomic
    echo "=============================================="
    echo ">> Procesando $BASENAME"
    echo "=============================================="

    # Descomprimir solo si no existe ya el .fna (se deja junto al .gz en 00_genomes)
    FNA="$GENOMES_DIR/${BASENAME}.fna"
    if [ ! -f "$FNA" ]; then
        gunzip -k "$FNA_GZ"
    fi

    # ---------------------------
    # 1. SeqKit
    # ---------------------------
    echo "  -> SeqKit stats..."
    STATS=$(seqkit stats -a -T "$FNA" | tail -n1)
    NUM_SEQS=$(echo "$STATS" | awk -F'\t' '{print $4}')
    SUM_LEN=$(echo "$STATS"  | awk -F'\t' '{print $5}')
    N50=$(echo "$STATS"      | awk -F'\t' '{print $12}')
    GC=$(echo "$STATS"       | awk -F'\t' '{print $NF}')

    # ---------------------------
    # 2. BUSCO
    # ---------------------------
    echo "  -> BUSCO (puede tardar)..."
    conda activate "$BUSCO_ENV"

    BUSCO_OUT="busco_${BASENAME}"
    rm -rf "$RESULTS_DIR/$BUSCO_OUT"

    busco -i "$FNA" \
          -o "$BUSCO_OUT" \
          -l "$LINEAGE" \
          -m genome \
          -c "$THREADS" \
          --out_path "$RESULTS_DIR" \
          -f > "$RESULTS_DIR/${BASENAME}_busco_log.txt" 2>&1

    conda deactivate

    SUMMARY_FILE=$(find "$RESULTS_DIR/$BUSCO_OUT" -name "short_summary*.txt" 2>/dev/null | head -n1)

    if [ -n "$SUMMARY_FILE" ]; then
        C_PCT=$(grep -oP 'C:\K[0-9.]+(?=%)' "$SUMMARY_FILE")
        S_PCT=$(grep -oP 'S:\K[0-9.]+(?=%)' "$SUMMARY_FILE")
        D_PCT=$(grep -oP 'D:\K[0-9.]+(?=%)' "$SUMMARY_FILE")
        F_PCT=$(grep -oP 'F:\K[0-9.]+(?=%)' "$SUMMARY_FILE")
        M_PCT=$(grep -oP 'M:\K[0-9.]+(?=%)' "$SUMMARY_FILE")
    else
        echo "  [!] BUSCO no genero resumen para $BASENAME. Revisa $RESULTS_DIR/${BASENAME}_busco_log.txt"
        C_PCT="NA"; S_PCT="NA"; D_PCT="NA"; F_PCT="NA"; M_PCT="NA"
    fi

    echo "$BASENAME,$NUM_SEQS,$SUM_LEN,$N50,$GC,$C_PCT,$S_PCT,$D_PCT,$F_PCT,$M_PCT" >> "$SUMMARY"
    echo "  -> Listo: $BASENAME (BUSCO completo: ${C_PCT}%)"
    echo ""
done

echo "=============================================="
echo "QC terminado. Resumen en: $SUMMARY"
echo "Logs y resultados detallados de BUSCO en: $RESULTS_DIR/"
echo "=============================================="
column -s',' -t "$SUMMARY"

