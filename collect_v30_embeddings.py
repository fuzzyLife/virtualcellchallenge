"""Fetch the pinned original SLIM STRING artifact; verify size, SHA256 and Git blob."""
from pathlib import Path
import urllib.request,hashlib,json
ROOT=Path(__file__).resolve().parent;D=ROOT/'reference/v30_structured_response';REV='5a7e9ade5d0a6b6331e6dbc81181450605047bcc';NAME='gene_string_embeddings.v0.3.h5';SHA='789416877b8701ef6f800106d26bf7bb97ea8e72744e6ab93e24933a717f247d';BLOB='8282fd198c2713f6d0a823c88211b6d76308a7df';SIZE=19502608

def main():
 D.mkdir(parents=True,exist_ok=True);p=D/NAME;url=f'https://raw.githubusercontent.com/RasmussenLab/SLIM/{REV}/data/{NAME}'
 if not p.exists():
  body=urllib.request.urlopen(url,timeout=120).read(SIZE+1);assert len(body)==SIZE;assert hashlib.sha256(body).hexdigest()==SHA;p.write_bytes(body)
 body=p.read_bytes();assert len(body)==SIZE and hashlib.sha256(body).hexdigest()==SHA;assert hashlib.sha1(f'blob {len(body)}\0'.encode()+body).hexdigest()==BLOB
 license_path=D/'LICENSE'
 if not license_path.exists():license_path.write_bytes(urllib.request.urlopen(f'https://raw.githubusercontent.com/RasmussenLab/SLIM/{REV}/LICENSE',timeout=60).read())
 receipt={'url':url,'bytes':SIZE,'revision':REV,'git_blob_sha1':BLOB,'sha256':SHA,'source':'SLIM author STRING embeddings; first64 coordinates match author default'}; (D/'embeddings.json').write_text(json.dumps(receipt,indent=2)+'\n');print('Original STRING artifact verified.')
if __name__=='__main__':main()
