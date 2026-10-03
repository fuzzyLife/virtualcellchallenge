"""Download and verify the original authors' genome-wide CRISPR pseudobulk data.

The smaller essential-gene reference has zero overlap with this competition panel.
This CC-BY-4.0 file covers genome-wide K562 interventions, without requiring the
65.8 GB single-cell file. Only the compact pseudobulk release is downloaded.
"""
from pathlib import Path
import argparse
import hashlib
import json
import time
import urllib.request

ROOT = Path(__file__).resolve().parent
ARTICLE_URL = 'https://api.figshare.com/v2/articles/20029387'
FILENAME = 'K562_gwps_raw_bulk_01.h5ad'
FILE_ID = 35774443
EXPECTED_MD5 = '4570b53c9d62ff6df281e622f0350060'


def md5_file(path):
    h = hashlib.md5()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=ROOT)
    root = parser.parse_args().root.resolve()
    dest = root/'reference'
    dest.mkdir(exist_ok=True)
    with urllib.request.urlopen(ARTICLE_URL,timeout=60) as response:
        metadata = json.load(response)
    (dest/'replogle_figshare_metadata.json').write_text(json.dumps(metadata,indent=2),encoding='utf-8')
    info = next(x for x in metadata['files'] if x['id']==FILE_ID and x['name']==FILENAME)
    path = dest/FILENAME
    expected = info.get('computed_md5') or info.get('supplied_md5')
    if expected != EXPECTED_MD5:
        raise RuntimeError('Pinned source checksum changed; review the public reference before rebuilding')
    if path.exists() and path.stat().st_size==info['size'] and md5_file(path)==expected:
        print(f'Reusing verified public reference: {path}',flush=True)
    else:
        partial = path.with_suffix('.h5ad.partial')
        print(f'Downloading {FILENAME}: {info["size"]:,} bytes',flush=True)
        for attempt in range(4):
            offset = partial.stat().st_size if partial.exists() else 0
            headers = {'Range':f'bytes={offset}-'} if offset else {}
            try:
                request = urllib.request.Request(info['download_url'],headers=headers)
                with urllib.request.urlopen(request,timeout=120) as response:
                    if offset and response.status!=206:
                        offset=0
                    mode='ab' if offset else 'wb'
                    with partial.open(mode) as f:
                        last = time.time()
                        while True:
                            block=response.read(4*1024*1024)
                            if not block:break
                            f.write(block)
                            if time.time()-last>=10:
                                print(f'Downloaded {f.tell():,}/{info["size"]:,} bytes',flush=True)
                                last=time.time()
                if partial.stat().st_size!=info['size']:
                    raise RuntimeError('Public reference download truncated')
                actual=md5_file(partial)
                if expected and actual!=expected:
                    raise RuntimeError('Public reference checksum mismatch')
                partial.replace(path)
                break
            except Exception:
                if attempt==3:raise
                time.sleep(2**attempt)
    provenance = {'source':'Replogle et al. 2022 original processed genome-wide Perturb-seq release',
                  'article_url':ARTICLE_URL,'article_doi':metadata.get('doi'),
                  'file_id':info['id'],'file_name':FILENAME,'download_url':info['download_url'],
                  'bytes':info['size'],'md5':expected,'license':metadata['license'],
                  'challenge_ground_truth_used':False}
    (dest/'provenance.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
    print('Public reference verified',flush=True)


if __name__=='__main__':
    main()
