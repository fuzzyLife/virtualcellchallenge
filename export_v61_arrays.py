"""Copy checksum-verified V61 state/response/count arrays into durable project storage."""
import shutil
import v61_pooled_state_response as v
m=v.m

def main():
 v.activate();files=[]
 for source in sorted(v.CACHE.rglob('*.npz')):
  relative=source.relative_to(v.CACHE);dest=v.A/'arrays'/relative;digest=m.sha(source);dest.parent.mkdir(parents=True,exist_ok=True)
  if dest.exists():assert m.sha(dest)==digest
  else:shutil.copyfile(source,dest);assert m.sha(dest)==digest
  files.append(dict(file=str(relative),bytes=dest.stat().st_size,sha256=digest))
 for source in sorted(v.parent.PARENT_CACHE.joinpath('raw').glob('*.npz')):
  relative=m.Path('reused_raw')/source.name;dest=v.A/'arrays'/relative;digest=m.sha(source);dest.parent.mkdir(parents=True,exist_ok=True)
  if dest.exists():assert m.sha(dest)==digest
  else:shutil.copyfile(source,dest);assert m.sha(dest)==digest
  files.append(dict(file=str(relative),bytes=dest.stat().st_size,sha256=digest))
 m.write(v.A/'array_exports.json',dict(files=files,total_bytes=sum(f['bytes'] for f in files),exporter_sha256=m.sha(__file__),scientific_plan_sha256=m.sha(v.A/'plan.json'),scope='V61 state bases, source mechanisms, public predictions and count statistics. Verified V59 raw CPM caches are included under reused_raw; production arrays and inputs are already inside the project.',D3_opened=False))
 print('V61 durable array export:',len(files),'files',sum(f['bytes'] for f in files),'bytes',flush=True)
if __name__=='__main__':main()
