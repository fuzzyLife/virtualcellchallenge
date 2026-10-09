"""Fetch pinned official scoring evidence; audit archived V11 diagnostics without refitting."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, re, sys
from concurrent.futures import ThreadPoolExecutor
import requests

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "analysis/evaluation_criteria_audit"
PIN = "0ce25cf4c495e76e60dfca4ddc8f822507aa0e18"
FILES = ["LICENSE", "pyproject.toml", "docs/vcc2026_metrics/vcc2026-metrics.md",
 "docs/vcc2026_metrics/vcc2026-metrics-brief.md", "src/cell_eval2/de.py",
 "src/cell_eval2/competition.py", "src/cell_eval2/scoring.py",
 "src/cell_eval2/configs/vcc2026.yaml", "src/cell_eval2/de_compute.py",
 "src/cell_eval2/metrics/de.py", "src/cell_eval2/metrics/direction.py",
 "src/cell_eval2/metrics/discrimination.py", "src/cell_eval2/metrics/delta.py"]

def sha(path):
 return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def get(url):
 response = requests.get(url, timeout=60)
 response.raise_for_status()
 return response

def main():
 OUT.mkdir(parents=True, exist_ok=True)
 head = get("https://api.github.com/repos/ArcInstitute/cell-eval2/commits/main").json()
 if head["sha"] != PIN:
  raise RuntimeError("Official scorer changed; review the new revision before continuing")
 def fetch(name):
  url = f"https://raw.githubusercontent.com/ArcInstitute/cell-eval2/{PIN}/{name}"
  target = OUT / "official" / name
  target.parent.mkdir(parents=True, exist_ok=True)
  target.write_bytes(get(url).content)
  old = ROOT / "reference/v8_metric_review" / name
  return name, {"url": url, "sha256": sha(target),
    "same_as_v8_review": sha(old) == sha(target) if old.exists() else None}
 with ThreadPoolExecutor(max_workers=4) as pool:
  sources = dict(pool.map(fetch, FILES))
 url = "https://virtualcellchallenge.org/evaluation"
 page = get(url)
 assets = re.findall(r'<script[^>]+src="([^"]+)"', page.text)
 chunk = next(a for a in assets if "/evaluation/page-" in a)
 chunk_url = "https://virtualcellchallenge.org" + chunk
 body = get(chunk_url)
 # Keep a hash and selected short strings, not a republished copy of the website.
 rules = ["The overall score is a single unweighted average over the six metrics and the three cell contexts",
          "4,750,000,000 stored entries", "400 cells for every perturbation", "18,533 genes",
          "1,000,000 counts", "No non-targeting rows."]
 website = {"url": url, "html_sha256": hashlib.sha256(page.content).hexdigest(),
  "asset_url": chunk_url, "asset_sha256": hashlib.sha256(body.content).hexdigest(),
  "checked_strings": {r: r in body.text for r in rules}}
 record = {"retrieved_utc": datetime.now(timezone.utc).isoformat(),
  "official_commit": PIN, "official_commit_date": head["commit"]["committer"]["date"],
  "existing_v8_review_commit": json.loads((ROOT/"analysis/v8_official_metric_review.json").read_text())["commit"],
  "sources": sources, "website": website}
 (OUT / "sources.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
 print(json.dumps({"official_commit": PIN, "unchanged_prior_source_files": sum(v["same_as_v8_review"] is True for v in sources.values()), "website_checks": website["checked_strings"]}, indent=2))

if __name__ == "__main__":
 sys.stdout.reconfigure(encoding="utf-8")
 main()
