#!/usr/bin/env python3
"""Valida resultados de cribado MASIVE-ALS sin modificar los datos de entrada."""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


def numeric(value: str) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def validate_csv(path: Path) -> dict:
    report = {
        "path": str(path),
        "exists": path.is_file(),
        "rows": 0,
        "valid_rows": 0,
        "invalid_rows": 0,
        "columns": [],
        "targets": {},
        "duplicate_keys": 0,
        "issues": [],
    }
    if not path.is_file():
        report["issues"].append("archivo no encontrado")
        return report

    with path.open(newline="", encoding="utf-8-sig", errors="replace") as handle:
        reader = csv.DictReader(handle)
        report["columns"] = reader.fieldnames or []
        required = {"ligand", "target"}
        energy_column = next(
            (name for name in ("energy", "best_energy", "vina_affinity", "vina_affinity_min") if name in report["columns"]),
            None,
        )
        missing = sorted(required - set(report["columns"]))
        if missing:
            report["issues"].append("faltan columnas: " + ", ".join(missing))
        if energy_column is None:
            report["issues"].append("no se encontró una columna de energía reconocida")

        keys = Counter()
        for row in reader:
            report["rows"] += 1
            ligand = (row.get("ligand") or "").strip()
            target = (row.get("target") or "").strip()
            key = (ligand, target, (row.get("seed") or row.get("best_seed") or "").strip())
            keys[key] += 1
            valid = bool(ligand and target)
            if energy_column is not None:
                valid = valid and numeric((row.get(energy_column) or "").strip())
            if valid:
                report["valid_rows"] += 1
                report["targets"][target] = report["targets"].get(target, 0) + 1
            else:
                report["invalid_rows"] += 1

        report["duplicate_keys"] = sum(count - 1 for count in keys.values() if count > 1)
        if report["duplicate_keys"]:
            report["issues"].append("hay claves duplicadas")
        if report["invalid_rows"]:
            report["issues"].append("hay filas sin ligando, diana o energía numérica")
    report["status"] = "válido" if report["rows"] and not report["issues"] else "requiere revisión"
    return report


def inspect_output_directory(path: Path) -> dict:
    csvs = sorted(path.rglob("*.csv")) if path.is_dir() else []
    logs = sorted(path.rglob("*.log")) if path.is_dir() else []
    poses = sorted(path.rglob("*_out.pdbqt")) if path.is_dir() else []
    return {
        "path": str(path),
        "exists": path.is_dir(),
        "csv_files": [str(item) for item in csvs],
        "log_count": len(logs),
        "pose_count": len(poses),
        "csv_reports": [validate_csv(item) for item in csvs],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", type=Path, action="append", default=[])
    parser.add_argument("--output-dir", type=Path, action="append", default=[])
    parser.add_argument("--report", type=Path, default=Path("analysis/validacion_cribado.json"))
    args = parser.parse_args()

    report = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "scientific_rule": "No se aceptan resultados sin energía numérica ni se mezclan plataformas/protocolos sin trazabilidad.",
        "csv_reports": [validate_csv(path) for path in args.csv],
        "output_directories": [inspect_output_directory(path) for path in args.output_dir],
    }
    report["status"] = "válido" if all(
        item.get("status") == "válido" for item in report["csv_reports"]
    ) and not any(
        item.get("exists") and item.get("pose_count", 0) == 0 and item.get("csv_files")
        for item in report["output_directories"]
    ) else "pendiente o requiere revisión"
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
