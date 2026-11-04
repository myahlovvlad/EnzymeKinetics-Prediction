# kcat: completed training benchmark

Selected by validation: **mlp**.

| Model | Validation MAE | Test MAE | Test RMSE | Test R² | Unseen-sequence MAE |
| --- | --- | --- | --- | --- | --- |
| mean | 1.2533 | 1.2508 | 1.6291 | -0.0000 | 1.3506556749343872 |
| ridge | 1.1966 | 1.1301 | 1.4744 | 0.1809 | 1.1361552476882935 |
| mbpls | 1.1937 | 1.1317 | 1.4810 | 0.1735 | 1.1421483564914638 |
| opls | 1.1930 | 1.1568 | 1.5192 | 0.1304 | 1.2213996430065868 |
| ec_mbpls | 1.2346 | 1.0559 | 1.4006 | 0.2608 | 1.2629112739446944 |
| mlp | 1.1119 | 0.9883 | 1.3154 | 0.3481 | 1.0951789617538452 |

All errors are in log10 target units. Test was not used for model selection.
The unseen subset excludes exact sequences present anywhere in trainval; homologues can still overlap.
Single seed; no claim of calibrated uncertainty or superiority to published models.

Rows: {'train': 16731, 'validation': 4104, 'test': 2316}. Wall time: 77.2 s.
