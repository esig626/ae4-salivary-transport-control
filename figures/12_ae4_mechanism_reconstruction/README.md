# Phase 12 diagnostic figures

| Figure | Generated from | Interpretation |
| --- | --- | --- |
| `ae4_candidate_predictions.png` | `model_comparison.csv` | Endpoint and integrated held-out ratios for the primary C1--C8 rows; Round 2 has separate diagnostics |
| `ae4_closure_prediction_tradeoff.png` | `model_comparison.csv` | Primary-candidate WT closure versus held-out prediction, with overlapping C4--C8 defaults grouped for readability |
| `ae4_k_fraction_closure.png` | `round2_envelope.csv` | Exact fixed-row penalty from transported K fraction |
| `ae4_pka_activation_envelope.png` | `pka_activation_envelope.csv` | Positive PKA activation moves tested rest-pass predictions in the adverse direction |
| `ae4_rest_capacity_tradeoff.png` | `rest_capacity_envelope.csv` | Refined finite moving-rest capacity bands and primary chloride/pH gates |

These images are deterministic audit diagnostics, not primary experimental
figures or manuscript candidates. The controlling numerical values are the
machine-readable tables, not values read from image pixels.

The closure panel can be regenerated without repeating simulations:

```bash
PYTHONPATH=src python -c 'from ae4_mechanism_reconstruction.run_analysis import render_closure_prediction_tradeoff_from_frozen_results as render; render()'
```
