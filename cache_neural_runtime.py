"""Cache the official Torch runtime locally to avoid loading DLLs over SMB.
Copies unchanged installed files into a task-specific temporary directory.
"""
from pathlib import Path
import json
import shutil
import sys
import tempfile
from build_submission import ROOT, sha256, write_json


def main():
    packages = Path(sys.prefix) / 'Lib/site-packages'
    source = packages / 'torch'
    if not source.is_dir():
        raise FileNotFoundError('Install the pinned official Torch wheel first')
    destination = Path(tempfile.gettempdir()) / 'vcc_gene300_torch_210_cu126'
    destination.mkdir(parents=True, exist_ok=True)
    if shutil.disk_usage(destination).free < 8 * 1024**3 and not (destination / 'torch').exists():
        raise OSError('Local temporary disk needs at least 8 GiB free for the runtime cache')
    def ignore(directory, names):
        excluded = [n for n in names if n == '__pycache__' or n.endswith('.pyc')]
        if Path(directory) == source:
            excluded += [n for n in ('include', 'share') if n in names]
        return excluded
    for name in ['torch', 'functorch', 'torchgen', 'torch-2.10.0+cu126.dist-info']:
        print(f'Copying {name} into local task runtime cache', flush=True)
        shutil.copytree(packages / name, destination / name, dirs_exist_ok=True, ignore=ignore)
    hashes = {}
    for original in sorted((source / 'lib').glob('*.dll')):
        copied = destination / 'torch/lib' / original.name
        before, after = sha256(original), sha256(copied)
        if before != after:
            raise ValueError(f'Native binary checksum mismatch: {original.name}')
        hashes[original.name] = before
    out = ROOT / 'analysis/gene300_neural'
    out.mkdir(parents=True, exist_ok=True)
    write_json(out / 'runtime_config.json', {'cache_path': str(destination)})
    write_json(out / 'runtime_cache_audit.json', {'official_torch_version': '2.10.0+cu126',
        'all_native_binary_hashes_match': True, 'dll_sha256': hashes,
        'excluded': ['C++ include headers', 'share', 'Python bytecode'], 'binaries_modified': False})
    print('Local runtime ready; all copied native binaries are unchanged.', flush=True)


if __name__ == '__main__':
    main()
