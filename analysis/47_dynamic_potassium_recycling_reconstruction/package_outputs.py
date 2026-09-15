"""Deterministic lossless CSV compression for compact Git checkpoints."""
from common import HERE
import gzip
import hashlib
import json

def main():
    manifest={}
    for p in sorted((HERE/'output').glob('*/*.csv')):
        if p.name not in ('states.csv','trajectory.csv'):continue
        raw=p.read_bytes();compressed=gzip.compress(raw,compresslevel=9,mtime=0)
        gz=p.with_suffix(p.suffix+'.gz')
        if gz.exists():assert gz.read_bytes()==compressed
        else:gz.write_bytes(compressed)
        assert gzip.decompress(gz.read_bytes())==raw
        manifest[str(gz.relative_to(HERE))]=dict(uncompressed_bytes=len(raw),uncompressed_sha256=hashlib.sha256(raw).hexdigest(),compressed_sha256=hashlib.sha256(compressed).hexdigest())
    (HERE/'output/compression_manifest.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
    print('Lossless CSV archives:',len(manifest))

if __name__=='__main__':main()
