#!/usr/bin/env bash
# Verifica el protocolo documentado, por CLI y con receptor persistente.
set -uo pipefail

PY="$HOME/rescoring_env/bin/python"
LOCAL="/mnt/c/Users/Fredy/masive-als"
S="$LOCAL/analysis/rescoring_local/mmgbsa_openff_gb.py"
POSE="$LOCAL/gpu_dock/resultados_libreria/results_SOD1/CHEMBL4584906_out.pdbqt"
REC="$LOCAL/gpu_dock/SOD1.pdbqt"
FIJO="$LOCAL/analysis/rescoring_local/receptores_fijos/SOD1_fijo.pdb"
SM='N=C1N/C(=N/NC(=O)c2ccc(CN3C(=O)c4cccc5cccc3c45)cc2)c2ccccc21'

mkdir -p "$(dirname "$FIJO")"

echo "=== 1) preparar el receptor UNA vez (persistente, no /tmp) ==="
"$PY" "$S" --preparar-receptor "$REC" --pose "$POSE" --out "$FIJO"
md5sum "$FIJO"

echo
echo "=== 2) rescoring con receptor fijo, CUDA (2 corridas) ==="
for i in 1 2; do
  echo -n "  corrida $i: "
  start=$(date +%s.%N)
  MMGBSA_FINAL_PLATFORM=CUDA "$PY" "$S" --pose "$POSE" --receptor "$FIJO" \
      --fijo --smiles "$SM" --out json 2>/dev/null
  end=$(date +%s.%N)
  echo "    tiempo: $(echo "$end - $start" | bc) s"
done

echo
echo "=== 3) la misma corrida en Reference (para comparar) ==="
echo -n "  Reference: "
MMGBSA_FINAL_PLATFORM=Reference "$PY" "$S" --pose "$POSE" --receptor "$FIJO" \
    --fijo --smiles "$SM" --out json 2>/dev/null
