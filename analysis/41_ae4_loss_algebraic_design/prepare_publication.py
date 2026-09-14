"""Prepare compact UTF-8 Git tree entries; publication uses GitHub tools."""
from pathlib import Path
import hashlib,json

ROOT=Path(__file__).resolve().parents[2]
BATCHES=ROOT.parent/'task41_publication_batches'
AUDIT=ROOT/'results/41_ae4_loss_algebraic_design/publication_audit.json'


def main():
    paths=[]
    for directory in ('analysis/41_ae4_loss_algebraic_design','results/41_ae4_loss_algebraic_design'):
        paths.extend(p for p in (ROOT/directory).rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.png' and p!=AUDIT)
    paths.extend(ROOT/p for p in ('src/modern_full_model/nkcc1_capacity_limited.py',
        'src/modern_full_model/ae4_cacc_recruitment.py','src/modern_full_model/task41_selected.py',
        'tests/test_task41_nkcc1_capacity.py','tests/test_task41_ae4_cacc_recruitment.py','tests/test_task41_selected_driver.py'))
    paths=sorted(set(paths));manifest=[]
    for p in paths:
        data=p.read_bytes();data.decode('utf-8')
        manifest.append({'path':str(p.relative_to(ROOT)),'bytes':len(data),
            'sha256':hashlib.sha256(data).hexdigest(),
            'git_blob_sha1':hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()})
    audit={'repository':'esig626/ae4-salivary-transport-control',
        'branch':'codex/task-40-ae4-equal-cation-routing',
        'base_commit':'76f3144a4e662eb12d6da26c83ebf6f2ff47a59d',
        'base_tree':'ac130d6d97bf1d6b6b8b94996e469abb0d58d35a',
        'main_observed_before_publication':'fe51473d7cc201e2e04007aa62fbc1b45098732e',
        'local_parent_files_compared_with_remote_blobs':189,'changed_parent_files':[],
        'remote_write_method':'Connected GitHub integration; Git trees/commit/non-forced branch update',
        'scope':'New Task 41 files only; preserve every parent blob; no main merge',
        'artifact_policy':'UTF-8 source, tests, Markdown, JSON, CSV, logs and SVG. PNG is a local preview of the published SVG.',
        'files_excluding_this_audit':manifest}
    AUDIT.write_text(json.dumps(audit,indent=2,sort_keys=True)+'\n');paths.append(AUDIT)
    entries=[{'path':str(p.relative_to(ROOT)),'mode':'100644','type':'blob','content':p.read_text()} for p in paths]
    BATCHES.mkdir(exist_ok=True)
    groups=[];group=[];size=0
    for entry in entries:
        n=len(entry['content'].encode())
        if group and size+n>200000:groups.append(group);group=[];size=0
        group.append(entry);size+=n
    if group:groups.append(group)
    for i,group in enumerate(groups):
        (BATCHES/f'batch_{i:02d}.json').write_text(json.dumps(group,ensure_ascii=False))
    (BATCHES/'index.json').write_text(json.dumps({'batches':len(groups),'entries':len(entries),
        'bytes':sum(len(e['content'].encode()) for e in entries),'files':[e['path'] for e in entries]}))
    print(json.dumps({'batches':len(groups),'entries':len(entries),'directory':str(BATCHES)}))

if __name__=='__main__':main()
