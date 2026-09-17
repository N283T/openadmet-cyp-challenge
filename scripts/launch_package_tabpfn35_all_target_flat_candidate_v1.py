"""Launch and register the all-target TabPFN 3.5 flat candidate."""

import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REG = ROOT / "artifacts/runs/regression"
SHARED = ROOT / "artifacts/runs/shared"
CARUANA = REG / "20260912T074545177210Z-0d561e17"
inputs = {
    "reference.csv": ROOT / "data/raw/cyp-challenge-TEST-BLINDED.csv",
    "aliases.parquet": SHARED / "20260909T132734119233Z-a7063ccd/outputs/aliases.parquet",
    "base-candidate.csv": ROOT
    / "artifacts/candidates/regression/20260914T024735834387Z-f87e9914/predictions.csv",
    "jeremy-frozen-test.parquet": REG
    / "20260912T074409483212Z-f3024998/outputs/test-predictions.parquet",
    "jeremy-scores-test.parquet": REG
    / "20260912T074557646349Z-da857069/outputs/test-predictions.parquet",
    "jeremy-finetuned-test.parquet": REG
    / "20260912T074624192288Z-bd58b762/outputs/test-predictions.parquet",
    "chemeleon-pca-test.parquet": REG
    / "20260912T074558733332Z-98ab8ceb/outputs/test-predictions.parquet",
    "chemeleon-raw-test.parquet": REG
    / "20260917T124509252191Z-72581a49/outputs/test-predictions.parquet",
    "t35-fast-pub030-test.parquet": REG
    / "20260917T123836807365Z-b8a43a60/outputs/test-predictions.parquet",
    "t35-fast-jeremy-raw-test.parquet": REG
    / "20260917T123919466013Z-afe1b6fd/outputs/test-predictions.parquet",
    "t35-fast-anchor-test.parquet": REG
    / "20260917T124004117613Z-279d6c0f/outputs/test-predictions.parquet",
    "t35-fast-jeremy-scores-test.parquet": REG
    / "20260917T124139840429Z-396892be/outputs/test-predictions.parquet",
    "t35-full-pub030-test.parquet": REG
    / "20260917T124150857848Z-fafa7ed1/outputs/test-predictions.parquet",
    "t35-full-jeremy-raw-test.parquet": REG
    / "20260917T124251141007Z-ef60e265/outputs/test-predictions.parquet",
    "t35-full-anchor-test.parquet": REG
    / "20260917T124342046803Z-4f691c7d/outputs/test-predictions.parquet",
    "t35-full-jeremy-scores-test.parquet": REG
    / "20260917T124432805624Z-b4f422fc/outputs/test-predictions.parquet",
    "t35-fast-chemeleon-test.parquet": REG
    / "20260917T125149042397Z-c7343ca2/outputs/test-predictions.parquet",
    "owned-2c9-test.parquet": REG
    / "20260917T124709539239Z-15a821ad/outputs/test-predictions.parquet",
    "d6-atomic-test.parquet": REG
    / "20260917T115944893661Z-9791da86/outputs/atomic-test-predictions.parquet",
    "d6-fits.csv": REG / "20260917T123605451416Z-d9ca2bee/outputs/fits.csv",
    "d6-bank-config.json": ROOT / "configs/regression/flat-bank-2d6-tabpfn35-v3.json",
    "a4-fits.csv": REG / "20260917T123614914030Z-ca4f4ba1/outputs/fits.csv",
    "a4-bank-config.json": ROOT / "configs/regression/flat-bank-3a4-tabpfn35-v3.json",
    "caruana-compounds.parquet": CARUANA / "inputs/000-compounds.parquet",
    "caruana-test-mapping.parquet": CARUANA / "inputs/002-test_mapping.parquet",
    "caruana-class-features.parquet": CARUANA / "inputs/004-external-embeddings.parquet",
    "a4-class-quantile-065.txt": CARUANA
    / "outputs/CYP3A4_pIC50_direct_inhibition/class_quantile_0.65.txt",
    "a4-morgan-quantile-065.txt": CARUANA
    / "outputs/CYP3A4_pIC50_direct_inhibition/morgan_quantile_0.65.txt",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Package all-target candidate from target-specific flat TabPFN 3.5 atomic banks",
    "--config",
    str(ROOT / "configs/regression/tabpfn35-all-target-flat-candidate-v1.json"),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
result = subprocess.run(
    command + ["--", "python", "scripts/package_tabpfn35_all_target_flat_candidate.py"],
    cwd=ROOT,
    text=True,
    capture_output=True,
)
print(result.stdout, flush=True)
print(result.stderr, flush=True)
result.check_returncode()
run_ids = re.findall(r"run_id=(\S+)", result.stdout)
assert len(run_ids) == 1
run_id = run_ids[0]
subprocess.run(
    [
        "pixi",
        "run",
        "cyp",
        "candidate",
        "--track",
        "regression",
        "--file",
        str(REG / run_id / "outputs/candidate-placed.csv"),
        "--reference",
        str(ROOT / "data/raw/cyp-challenge-TEST-BLINDED.csv"),
        "--run",
        run_id,
        "--notes",
        (
            "All four endpoints updated by independently evaluated flat atomic banks; "
            "target distributions aligned to the prior favorable placement candidate "
            "while preserving the new ranks."
        ),
    ],
    cwd=ROOT,
    check=True,
)
