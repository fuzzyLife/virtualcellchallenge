"""Official-environment lossless validation/decompression only; no inference."""
from pathlib import Path
import argparse,tarfile,json
import zstandard as zstd
from vcc.vccfile import validate_vcc,read_vcc_meta

def main():
 p=argparse.ArgumentParser();p.add_argument('package');p.add_argument('destination');args=p.parse_args();source=Path(args.package);dest=Path(args.destination);assert not dest.exists();validate_vcc(source);meta=read_vcc_meta(source)
 with tarfile.open(source,'r:') as archive:
  item=archive.extractfile('pred.h5ad.zst');assert item is not None
  with item,dest.open('wb') as out:zstd.ZstdDecompressor().copy_stream(item,out)
 print(json.dumps(meta),flush=True)
if __name__=='__main__':main()
