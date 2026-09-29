"""Check the copied scientific inputs against their original Git blob hashes."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def verify():
    manifest = json.loads((ROOT / 'SOURCE_MANIFEST.json').read_text())
    for relative, expected in manifest['git_blobs'].items():
        data = (ROOT / relative).read_bytes()
        actual = hashlib.sha1(b'blob ' + str(len(data)).encode('ascii') + b'\0' + data).hexdigest()
        if actual != expected:
            raise RuntimeError('Missing or changed scientific input: ' + relative)
    print('Verified', len(manifest['git_blobs']), 'copied scientific files.')


if __name__ == '__main__':
    verify()
