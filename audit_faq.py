"""Review the live public FAQ and local final-round readiness without changing models.
Use --show to print extracted FAQ text for manual review; full answers are not archived.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse, hashlib, json, re, sys
import numpy as np
import requests

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "analysis/faq_audit"
FAQ = "https://virtualcellchallenge.org/faq"
ASSET = "/_next/static/chunks/334-291236033da9ac0d.js"
PIN = "8ab6ab697c9efc18fa4ba8981a0b1912ee7683823c4932aa26f9e100f697deab"


def sha(path):
 return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fetch(url):
 response = requests.get(url, timeout=45)
 response.raise_for_status()
 return response


def answer_end(content, start):
 # Stop at the containing FAQ object's closing brace, ignoring delimiters in strings.
 stack = []
 quote = None
 escaped = False
 pairs = {"}": "{", ")": "(", "]": "["}
 for i in range(start, len(content)):
  char = content[i]
  if quote is not None:
   if escaped: escaped = False
   elif char == chr(92): escaped = True
   elif char == quote: quote = None
   continue
  if char in [chr(34), chr(39), chr(96)]: quote = char
  elif char in "{([": stack.append(char)
  elif char in "})]":
   if not stack:
    if char == "}": return i
    raise ValueError("Unexpected FAQ object delimiter")
   assert stack.pop() == pairs[char]
 raise ValueError("Unterminated FAQ answer")


def extract(content):
 matches = list(re.finditer(r'question:("(?:\\.|[^"\\])*"),answer:', content))
 entries = []
 skip = {"div", "p", "a", "span", "br", "li", "ul", "b", "i", "strong", "_blank", "noopener noreferrer", "code"}
 for i, match in enumerate(matches):
  block = content[match.end():answer_end(content, match.end())]
  strings = []
  for literal in re.findall(r'"(?:\\.|[^"\\])*"', block):
   try: value = json.loads(literal)
   except json.JSONDecodeError: continue
   if value not in skip and not any(t in value for t in ["font-", "bg-", "text-", "padding:", "className"]):
    strings.append(value)
  entries.append({"question": json.loads(match.group(1)), "answer_text": " ".join(strings)})
 return entries


def main(show=False):
 html = fetch(FAQ)
 if ASSET not in html.text:
  raise RuntimeError("FAQ asset changed; inspect current page dependencies before updating the pin")
 asset = fetch("https://virtualcellchallenge.org" + ASSET)
 assert hashlib.sha256(asset.content).hexdigest() == PIN, "FAQ changed; review before updating evidence"
 entries = extract(asset.text)
 assert len(entries) == 50, len(entries)
 assert "use strict" not in entries[-1]["answer_text"]
 def answer(part):
  found = [e for e in entries if part.lower() in e["question"].lower()]
  assert len(found) == 1, part
  return found[0]["answer_text"]
 checks = {
  "release_october_22_2026": "October 22, 2026" in answer("What dates"),
  "final_submission_november_5_2026": "November 5, 2026" in answer("What dates"),
  "validation_and_final_target_panels_differ": "different panels" in answer("How much data"),
  "final_contexts_D_E_F": "D, E and F" in answer("How many cell lines"),
  "identities_not_published": "do not publish which cell lines" in answer("How many cell lines"),
  "no_separate_challenge_training_set": "no separate Challenge-specific training set" in answer("What data sets"),
  "control_guide_structure": "46 non-targeting guides" in answer("How much data"),
  "all_contexts_share_within_round_panel": "share the same perturbation panel" in answer("How much data"),
  "context_label_swaps_not_diagnosed_by_score": "nothing in the score" in answer("Do I submit one file"),
  "external_data_allowed_with_rights": "appropriate permissions" in answer("Can I use other data"),
  "matched_crispri_flex_assay": "CRISPRi and 10x Flex" in answer("How was the benchmark"),
  "two_daily_results": "two results per day" in answer("Can I try multiple approaches"),
  "one_prize_considered_entry": "Only one submission" in answer("Can I try multiple approaches"),
  "explicit_stored_zeros_count_toward_limit": "zero you saved explicitly still counts" in answer("rejected as too dense"),
  "validation_test_scores_not_directly_comparable": "different perturbation panels" in answer("compare my validation score"),
  "faq_high_end_wording_is_unbounded": "unbounded above" in answer("Why can a score go above"),
 }
 assert all(checks.values()), checks
 if show:
  for entry in entries:
   print(entry["question"], "\n", entry["answer_text"], "\n")
 spec = json.loads((ROOT/"manifest.json").read_text())
 v11 = ROOT/"experiments/v11_consensus_coverage/analysis"
 model = json.loads((v11/"reference_model.json").read_text())
 import pandas as pd
 genes = pd.read_csv(ROOT/"gene_names.csv").gene_name.to_numpy(dtype=str)
 targets = pd.read_csv(ROOT/"pert_counts.csv").target_gene.to_numpy(dtype=str)
 with np.load(v11/"reference_model_arrays.npz", allow_pickle=False) as arrays:
  np.testing.assert_array_equal(arrays["genes"], genes)
  np.testing.assert_array_equal(arrays["targets"], targets)
 with np.load(ROOT/"experiments/v7_multisource_adaptive/analysis/source_cache.npz", allow_pickle=False) as arrays:
  np.testing.assert_array_equal(arrays["targets"], targets)
 contexts = [r["context"] for r in model["context_original_inputs"]]
 assert contexts == spec["contexts"] == ["A", "B", "C"]
 inputs = {}
 for name, digest in model["dependencies"].items():
  assert sha(ROOT/name) == digest, name
  inputs[name] = digest
 publication = json.loads((v11/"publication_manifest.json").read_text())
 for name, digest in publication["generation_snapshot_sha256"].items():
  assert sha(v11/name) == digest, name
 previous = json.loads((ROOT/"analysis/evaluation_criteria_audit/audit_manifest.json").read_text())
 for name, digest in previous["files"].items():
  assert sha(ROOT/name) == digest, name
 frozen_plan = "analysis/post_v11_research/future_experiment_plan.json"
 assert sha(ROOT/frozen_plan) == previous["frozen_plan_sha256"]
 inputs[frozen_plan] = sha(ROOT/frozen_plan)
 for name in ["manifest.json", "gene_names.csv", "pert_counts.csv", "fit_v11_submission.py", "build_v11_submission.py", "build_submission.py"]:
  inputs[name] = sha(ROOT/name)
 independent = json.loads((v11/"independent_validation.json").read_text())
 assert independent["passed"] and independent["all_packaged_values_indices_pointers_and_axes_match"]
 assert independent["packaged_nonzero_values_checked"] < 4_750_000_000
 reviewed = {}
 for name, needle in [("build_submission.py", "assert set(obj.obs['context'].astype(str)) == {context}"),
  ("build_v11_submission.py", "np.testing.assert_array_equal(z['targets'],targets)"),
  ("fit_v11_submission.py", "np.testing.assert_array_equal(data['targets'],targets)")]:
  lines = (ROOT/name).read_text(encoding="utf-8").splitlines()
  reviewed[name] = {"line": next(i+1 for i,line in enumerate(lines) if needle in line), "sha256": sha(ROOT/name)}
 # Positive integer counts at valid cell depths cannot have per-cell CPM below 1.
 count = np.array([1, 1, 7], dtype=np.int64)
 depth = np.array([20_000, 1_000_000, 999_999], dtype=np.int64)
 assert np.all(count*1_000_000/depth >= 1)
 record = {"retrieved_utc": datetime.now(timezone.utc).isoformat(), "faq_url": FAQ,
  "html_sha256": hashlib.sha256(html.content).hexdigest(), "content_asset_url": "https://virtualcellchallenge.org"+ASSET,
  "content_asset_sha256": PIN, "questions_reviewed": len(entries), "checks": checks,
  "current_partition": spec["partition"], "current_context_order": contexts,
  "current_target_count": len(targets), "current_target_and_source_axes_match": True,
  "stored_entries_from_archived_full_validation": independent["packaged_nonzero_values_checked"],
  "stored_entry_limit_fraction": independent["packaged_nonzero_values_checked"]/4_750_000_000,
  "strict_per_cell_below_1_CPM_pruning_is_noop_at_valid_raw_count_depths": True,
  "script_sha256": sha(__file__),
  "source_code_evidence": reviewed, "input_sha256": inputs,
  "frozen_dependencies_verified": len(model["dependencies"]),
  "frozen_generation_snapshots_verified": len(publication["generation_snapshot_sha256"]),
  "prior_evaluation_audit_and_research_plan_unchanged": True,
  "scope": "Read-only FAQ and axis/provenance audit. Previous full-payload validation read, not rerun. No final inputs available or final-target overlap inferred."}
 OUT.mkdir(parents=True, exist_ok=True)
 (OUT/"audit.json").write_text(json.dumps(record,indent=2)+"\n",encoding="utf-8")
 print(json.dumps({"faq_questions_reviewed":len(entries),"semantic_checks":len(checks),
  "all_checks_passed":all(checks.values()), "current_model_partition":spec["partition"],
  "frozen_dependencies_verified":len(model["dependencies"]),
  "frozen_generation_snapshots_verified":len(publication["generation_snapshot_sha256"])},indent=2))

if __name__ == "__main__":
 sys.stdout.reconfigure(encoding="utf-8")
 parser = argparse.ArgumentParser()
 parser.add_argument("--show",action="store_true")
 main(parser.parse_args().show)
