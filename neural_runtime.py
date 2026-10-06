"""Activate the optional task-local cache of the unchanged Torch package."""
import json
from pathlib import Path
import sys
import tempfile


def activate_local_runtime():
    config = Path(__file__).resolve().parent / 'analysis/gene300_neural/runtime_config.json'
    if not config.exists():
        return
    cache = Path(json.loads(config.read_text(encoding='utf-8'))['cache_path']).resolve()
    if not cache.is_relative_to(Path(tempfile.gettempdir()).resolve()):
        raise ValueError('Torch cache must stay inside the local temporary directory')
    if (cache / 'torch/__init__.py').is_file():
        sys.path.insert(0, str(cache))
