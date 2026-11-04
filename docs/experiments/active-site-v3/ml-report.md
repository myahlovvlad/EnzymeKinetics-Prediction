# kcat: completed training benchmark

Selected by validation: **mlp**.

| Model | Validation MAE | Test MAE | Test RMSE | Test R² | Unseen-sequence MAE |
| --- | --- | --- | --- | --- | --- |
| mean | 1.2533 | 1.2508 | 1.6291 | -0.0000 | 1.3506556749343872 |
| ridge | 1.1960 | 1.1297 | 1.4728 | 0.1827 | 1.1371694803237915 |
| mlp | 1.1190 | 0.9731 | 1.2906 | 0.3724 | 1.0610873699188232 |

All errors are in log10 target units. Test was not used for model selection.
The unseen subset excludes exact sequences present anywhere in trainval; homologues can still overlap.
Single seed; no claim of calibrated uncertainty or superiority to published models.

Rows: {'train': 16731, 'validation': 4104, 'test': 2316}. Wall time: 165.5 s.
