"""Export the constructed principal parameter objects without integrating an ODE."""
from dataclasses import asdict, is_dataclass
from enum import Enum
from pathlib import Path
import argparse
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'model'))


def encode(value):
    if is_dataclass(value):
        return asdict(value)
    if isinstance(value, Enum):
        return value.value
    raise TypeError('Unsupported parameter object: ' + type(value).__name__)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    from verify_files import verify
    verify()
    from task52_model import make_parent, CARRIER
    parent = make_parent('CCH_IPR')
    freeze = json.loads((ROOT / 'parameters/parameter_and_case_freeze.json').read_text())
    record = {'whole_cell_parameters': parent.parameters,
              'Ae4_carrier_parameters': parent.ae4_parameters,
              'NBC_parameters': parent.nbc_parameters,
              'NHE1_carrier_amount_fmol': CARRIER,
              'case_matrix': freeze['case_matrix'],
              'central_onset_states': {g: freeze['projected_onsets']['central'][g]['state']
                                      for g in ('WT', 'KO')},
              'note': 'Original constructor evaluated; no trajectory integrated or parameter fitted.'}
    # Exclusive creation avoids overwriting a previous parameter export.
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(record, handle, default=encode, indent=2, allow_nan=False)
        handle.write('\n')


if __name__ == '__main__':
    main()
