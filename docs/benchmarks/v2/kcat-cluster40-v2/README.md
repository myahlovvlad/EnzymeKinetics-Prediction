# kcat: completed training benchmark

Selected by validation: **mlp**.

| Model | Validation MAE | Test MAE | Test RMSE | Test R² | Unseen-sequence MAE |
| --- | --- | --- | --- | --- | --- |
| mean | 1.3889 | 1.3031 | 1.7018 | -0.0066 | 1.3031260967254639 |
| ridge | 1.2939 | 1.2309 | 1.6073 | 0.1021 | 1.230928659439087 |
| mbpls | 1.2927 | 1.2253 | 1.6054 | 0.1042 | 1.2252515192849556 |
| opls | 1.3060 | 1.2111 | 1.5825 | 0.1296 | 1.2111305438287965 |
| ec_mbpls | 1.4450 | 1.3615 | 1.8119 | -0.1411 | 1.3615129482886157 |
| mlp | 1.2533 | 1.2057 | 1.5888 | 0.1227 | 1.2056688070297241 |

All errors are in log10 target units. Test was not used for model selection.
The unseen subset excludes exact sequences present anywhere in trainval; homologues can still overlap.
Single seed; no claim of calibrated uncertainty or superiority to published models.

Rows: {'train': 14571, 'validation': 4182, 'test': 4398}. Wall time: 69.6 s.
