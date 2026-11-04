# kcat: completed training benchmark

Selected by validation: **opls**.

| Model | Validation MAE | Test MAE | Test RMSE | Test R² | Unseen-sequence MAE |
| --- | --- | --- | --- | --- | --- |
| mean | 1.3103 | 1.3092 | 1.7000 | -0.0239 | 1.3108346462249756 |
| ridge | 1.2676 | 1.2795 | 1.6551 | 0.0294 | 1.279571294784546 |
| mbpls | 1.2223 | 1.2439 | 1.6119 | 0.0795 | 1.2443032145951771 |
| opls | 1.2149 | 1.2451 | 1.6174 | 0.0731 | 1.2468767583393832 |
| ec_mbpls | 1.2223 | 1.2439 | 1.6119 | 0.0795 | 1.2443032145951771 |
| mlp | 1.2598 | 1.2302 | 1.5911 | 0.1030 | 1.2318156957626343 |

All errors are in log10 target units. Test was not used for model selection.
The unseen subset excludes exact sequences present anywhere in trainval; homologues can still overlap.
Single seed; no claim of calibrated uncertainty or superiority to published models.

Rows: {'train': 17924, 'validation': 1883, 'test': 3344}. Wall time: 78.6 s.
