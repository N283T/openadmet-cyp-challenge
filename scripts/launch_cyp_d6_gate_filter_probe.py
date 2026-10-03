"""Launch the offline evaluation of a gate-filtered D6 call vector."""

import argparse
import subprocess
from pathlib import Path

DATA = "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
SPLITS = "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs"
POOL_RUN = "artifacts/runs/classification/20261002T215345688167Z-7920b49a/outputs"
TESTLIKE = "artifacts/tracked/cyp-d6-testlike-pool-20261002T215345688167Z-7920b49a"
SEP15_BLIND = "artifacts/runs/classification/20260924T045344348206Z-7e2db311/inputs"
NEIGHBOURS = "artifacts/runs/classification/20261002T223515148996Z-ccb05cf3/outputs"
SUBMISSIONS = "artifacts/submissions/classification"
PAYLOADS = {
    "payload_sep13.parquet": "20260913T113807374506Z-77f4825a",
    "payload_sep14.parquet": "20260914T020759952207Z-37dd7bd8",
    "payload_sep15.parquet": "20260915T130004513006Z-75b25070",
    "payload_sep18.parquet": "20260918T094451131332Z-f23dd61b",
    "payload_sep23.parquet": "20260923T180233070503Z-552d9af4",
    "payload_oct01.parquet": "20261001T141041310215Z-c7faa79f",
    "payload_oct02.parquet": "20261002T092348406577Z-4ea606c1",
    "payload_oct03.parquet": "20261002T223229553466Z-a08301d2",
}
INPUTS = {
    "testlike-pool-oof.parquet": f"{TESTLIKE}/oof_scores.parquet",
    "test_scores.parquet": f"{POOL_RUN}/test_scores.parquet",
    "sep15_blind_scores.parquet": f"{SEP15_BLIND}/blind-scores.parquet",
    "assignments.parquet": f"{SPLITS}/assignments.parquet",
    "test_mapping.parquet": f"{DATA}/test_mapping.parquet",
    "blind_neighbour_rates.parquet": f"{NEIGHBOURS}/blind_neighbour_rates.parquet",
    **{alias: f"{SUBMISSIONS}/{run}/predictions.parquet" for alias, run in PAYLOADS.items()},
}
PURPOSE = (
    "Evaluate a gate-filtered D6 call vector offline against the deployed rank-mean recipe at "
    "an unchanged call count: four p_active eligibility thresholds plus the product and "
    "square-root product forms, scored on a rate grid extended to 50% under an "
    "official-labels-only reweighted metric and a composition-matched mixture metric on both "
    "fixed split families, with a paired group bootstrap, the eight-transition historical "
    "concordance check and a blind swap inventory; read-only, no candidate or upload"
)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/classification/d6-gate-filter-probe-v1.json")
    args = parser.parse_args()
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        PURPOSE,
        "--config",
        str(root / args.config),
        "--input-mode",
        "hardlink",
    ]
    for alias, relative in INPUTS.items():
        path = root / relative
        if not path.is_file():
            raise FileNotFoundError(path)
        command.extend(["--input", f"{alias}={path}"])
    subprocess.run(
        command + ["--", "python", "scripts/evaluate_cyp_d6_gate_filter_probe.py"],
        cwd=root,
        check=True,
    )


if __name__ == "__main__":
    main()
