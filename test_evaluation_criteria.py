"""Execute official DE functions on fixtures and quantify archived validation omissions.
This is a metric-layer audit, not an end-to-end run of the production scoring service.
"""
from pathlib import Path
import importlib, json, logging, math, sys, types, platform
import scipy
import numpy as np
import pandas as pd
import polars as pl
from scipy import stats, sparse
from audit_evaluation_criteria import ROOT, OUT, PIN, sha
from v6_methods import direction_metrics, de
from audit_transfer_validation import pds


def official_modules():
 sources = json.loads((OUT / "sources.json").read_text())
 assert sources["official_commit"] == PIN
 for name, record in sources["sources"].items():
  assert sha(OUT / "official" / name) == record["sha256"], name
 # Avoid importing the package's unrelated CLI/GPU stack. These namespace packages
 # load the three unchanged official modules, including their real relative imports.
 package = "_criteria_official"
 for name, folder in [(package, "src/cell_eval2"), (package+".metrics", "src/cell_eval2/metrics")]:
  module = types.ModuleType(name)
  module.__path__ = [str(OUT / "official" / folder)]
  sys.modules[name] = module
 return importlib.import_module(package+".metrics.de"), importlib.import_module(package+".metrics.direction")


def fixture(n=14, real_calls=None, pred_calls=None):
 genes = np.array(["TARGET"] + [f"G{i:03}" for i in range(n)])
 def side(calls):
  q = np.ones(n+1)
  q[0] = .001
  q[1:1+calls] = .001
  return {"eligible": np.ones(n+1, bool), "q": q, "p": q/2, "lfc": np.ones(n+1)}
 return genes, side(n if real_calls is None else real_calls), side(n if pred_calls is None else pred_calls)


def frame(genes, values):
 keep = values["eligible"]
 return pl.DataFrame({"target": ["TARGET"] * int(keep.sum()), "feature": genes[keep],
  "log2_fold_change": values["lfc"][keep], "p_adj": values["q"][keep], "p_value": values["p"][keep]})


def evaluate(dm, direction, genes, truth, pred):
 args = {"de_real": frame(genes, truth), "de_pred": frame(genes, pred)}
 return {"fidelity": direction.de_direction_fidelity_yield_raw(**args)["TARGET"],
  "reach": direction.de_direction_reach(corrected=False, **args)["TARGET"],
  "jaccard": dm.de_sig_jaccard(**args)["TARGET"],
  "nmae": dm.de_lfc_nmae(**args).get("TARGET", math.nan)}


def number(value):
 if math.isnan(float(value)): return None
 if math.isinf(float(value)): return str(value)
 return float(value)


def main():
 logging.getLogger("_criteria_official").setLevel(logging.ERROR)
 dm, direction = official_modules()
 cases = []
 def record(label, genes, truth, pred, mismatches=()):
  with np.errstate(divide="ignore", invalid="ignore"):
   old = direction_metrics(pred, truth, genes, 0)
  exact = evaluate(dm, direction, genes, truth, pred)
  for key in exact:
   if key in mismatches:
    assert not np.isclose(old[key], exact[key], equal_nan=True), (label, key)
   else:
    np.testing.assert_allclose(old[key], exact[key], equal_nan=True, err_msg=label+":"+key)
  cases.append({"case": label, "local": {k: number(old[k]) for k in exact},
   "official": {k: number(v) for k,v in exact.items()}, "expected_mismatches": list(mismatches)})
  return exact
 for n in [9, 10]:
  genes, truth, pred = fixture(n)
  pred["lfc"][:] = 0
  result = record(f"NMAE gate {n} after own-gene exclusion", genes, truth, pred, ("nmae",) if n == 9 else ())
  assert math.isnan(result["nmae"]) if n == 9 else result["nmae"] == 1
 genes, truth, pred = fixture(14, real_calls=4, pred_calls=14)
 result = record("FID rewards correct signs beyond reference-significant genes", genes, truth, pred, ("nmae",))
 assert result["fidelity"] == 1 and result["jaccard"] == 4/14
 genes, truth, pred = fixture(14, real_calls=14, pred_calls=4)
 result = record("Reach includes predicted-nonsignificant reference responses", genes, truth, pred)
 assert result["reach"] == 1 and result["fidelity"] == 4/14
 genes, truth, pred = fixture(10)
 pred["lfc"][1] = -2
 pred["q"][1] = .00001
 result = record("Reach can recover to 90 percent purity after a leading miss", genes, truth, pred)
 assert result["reach"] == 1
 genes, truth, pred = fixture(14)
 truth["lfc"][:] = 0
 result = record("Significant real genes with no adjudicable direction give reach zero", genes, truth, pred, ("reach", "nmae"))
 assert result["reach"] == 0 and math.isnan(result["nmae"])
 # Own-gene-only changes cannot improve the four DE metrics (own stays in BH earlier).
 genes, truth, pred = fixture(14)
 before = evaluate(dm, direction, genes, truth, pred)
 pred["lfc"][0] = -100
 after = record("Own-target fold change is excluded from all four DE metrics", genes, truth, pred)
 assert before == after
 rng = np.random.default_rng(20261009)
 control = rng.normal(size=100)
 measured = rng.normal(size=(12,100))
 constant = np.repeat(rng.normal(size=(1,100)),12,axis=0)
 np.testing.assert_allclose(pds(constant, measured, control), .5)
 # Analytic CPU-Wilcoxon comparison: scanpy's rank z-score has no continuity term.
 # This does NOT identify which DE backend the hosted service currently uses.
 boundary = None
 for n in range(4, 41):
  for shift in range(1,n):
   a, b = np.arange(n), np.arange(n)+shift
   pc = stats.mannwhitneyu(a,b,method="asymptotic",use_continuity=True).pvalue
   pn = stats.mannwhitneyu(a,b,method="asymptotic",use_continuity=False).pvalue
   if pn < .05 <= pc:
    # Equal depths make this an exact two-gene count fixture, including BH.
    x, c = sparse.csr_matrix(np.c_[a,100-a]), sparse.csr_matrix(np.c_[b,100-b])
    local = de(x,c,np.array(["G", "BACKGROUND"]))
    np.testing.assert_allclose(local["q"],pc)
    boundary = {"cells_per_group": n, "shift": shift, "local_continuity_q": float(pc),
     "no_continuity_q": float(pn), "scope": "Analytic scanpy-style comparator; production backend not run"}
    break
  if boundary is not None: break
 assert boundary is not None
 v11 = ROOT / "experiments/v11_consensus_coverage/analysis"
 inputs = {}
 changes = []
 for stage in ["development", "confirmation"]:
  path = v11 / f"{stage}_scores.csv"
  inputs[str(path.relative_to(ROOT))] = sha(path)
  saved = pd.read_csv(path)
  for keys, group in saved.groupby(["candidate", "reference_mode"]):
   eligible = group.reference_calls >= 10
   retained = group[eligible & group.nmae.notna()]
   # Saved counts are finite epsilon-LFC outputs: ref_calls equals this gate size.
   filtered = group.assign(nmae=group.nmae.where(eligible))
   def aggregate(df):
    return df.groupby(["family", "context"]).nmae.mean().groupby("family").mean().mean()
   old, new = aggregate(group), aggregate(filtered)
   changes.append({"stage": stage, "candidate": keys[0], "reference_mode": keys[1],
    "rows": len(group), "previously_scored_small_gates": int(((group.reference_calls < 10) & group.nmae.notna()).sum()),
    "zero_reference_calls": int((group.reference_calls == 0).sum()),
    "official_gate_eligible_rows": len(retained),
    "eligible_families": "|".join(sorted(retained.family.unique())),
    "eligible_contexts": "|".join(sorted(retained.context.unique())),
    "old_nmae_equal_family_context": old, "gate10_nmae_equal_family_context": new,
    "difference": new-old})
 pd.DataFrame(changes).to_csv(OUT/"archived_nmae_gate_audit.csv",index=False)
 status_path = v11 / "published_api_status.json"
 status = json.loads(status_path.read_text())
 keys = ["pds", "mse", "nmae", "fid", "reach", "jac"]
 scaled = {k: status["score_"+k] for k in keys}
 np.testing.assert_allclose(np.mean(list(scaled.values())),status["score_avg"],atol=1e-14)
 contributions = {k: v/6 for k,v in scaled.items()}
 checks = {}
 model = json.loads((v11/"reference_model.json").read_text())
 for name, digest in model["dependencies"].items():
  assert sha(ROOT/name) == digest, name
 checks["frozen_v11_dependencies"] = len(model["dependencies"])
 manifest = json.loads((v11/"publication_manifest.json").read_text())
 for name, digest in manifest["generation_snapshot_sha256"].items():
  assert sha(v11/name) == digest, name
 checks["frozen_generation_snapshots"] = len(manifest["generation_snapshot_sha256"])
 research = ROOT/"analysis/post_v11_research/future_experiment_plan.json"
 inputs[str(research.relative_to(ROOT))] = sha(research)
 result = {"passed": True, "script_sha256": sha(__file__), "official_commit": PIN,
  "official_module_scope": "Unchanged official de.py, metrics/de.py, metrics/direction.py executed; package CLI initialization bypassed, no hosted or full counts pipeline parity claim",
  "runtime": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "polars": pl.__version__, "scipy": scipy.__version__},
  "local_dependencies": {name: sha(ROOT/name) for name in ["v6_methods.py", "audit_transfer_validation.py", "audit_evaluation_criteria.py"]},
  "fixtures": cases, "constant_prediction_pds": .5, "continuity_boundary": boundary,
  "archive_input_sha256": inputs, "publication_input_sha256": sha(status_path),
  "published_overall": status["score_avg"], "scaled_components": scaled,
  "overall_contributions": contributions, "pds_fraction_of_net_overall": contributions["pds"]/status["score_avg"],
  "frozen_provenance_verified": checks}
 (OUT/"test_results.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
 print(json.dumps({"fixtures": len(cases), "continuity_boundary": boundary, "frozen_checks": checks,
  "archived_v11": [row for row in changes if row["candidate"] == "consensus20_noise65"],
  "pds_share": result["pds_fraction_of_net_overall"]},indent=2))

if __name__ == "__main__":
 sys.stdout.reconfigure(encoding="utf-8")
 main()
