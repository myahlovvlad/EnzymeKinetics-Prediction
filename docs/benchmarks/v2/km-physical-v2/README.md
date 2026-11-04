# km: completed training benchmark

Selected by validation: **mlp**.

| Model | Validation MAE | Test MAE | Test RMSE | Test R² | Unseen-sequence MAE |
| --- | --- | --- | --- | --- | --- |
| mean | 1.0138 | 1.0136 | 1.2838 | -0.0017 | 1.0149446725845337 |
| ridge | 0.8371 | 0.7932 | 1.0284 | 0.3572 | 0.8252847194671631 |
| mbpls | 0.8371 | 0.7987 | 1.0358 | 0.3479 | 0.8275995115294088 |
| opls | 0.8597 | 0.8271 | 1.0708 | 0.3031 | 0.8787063132526279 |
| ec_mbpls | 0.8199 | 0.7499 | 0.9815 | 0.4145 | 0.8357688729300028 |
| mlp | 0.7460 | 0.6794 | 0.8961 | 0.5119 | 0.7154613733291626 |

All errors are in log10 target units. Test was not used for model selection.
The unseen subset excludes exact sequences present anywhere in trainval; homologues can still overlap.
Single seed; no claim of calibrated uncertainty or superiority to published models.

Rows: {'train': 29557, 'validation': 7499, 'test': 4118}. Wall time: 96.7 s.
