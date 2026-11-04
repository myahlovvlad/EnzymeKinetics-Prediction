# H1-E4 — proposed matched-head development benchmark

Status: design only. No mass EnzGFM encoding, EnzGFM regression fitting or model ranking has been performed. Current runtime is ready; independent scientific validation is blocked by data lineage/units/licenses and split provenance. Existing official test remains unopened.

## Frozen comparisons and partitions

B0: existing handcrafted FeatureBuilder, including its common RDKit/Morgan, sequence, EC, conditions/missingness and titration blocks. B1: exactly B0 + frozen ESM-2 8M. B1-E: exactly B0 + frozen EnzGFM 1.5B. Geometry/physics remains unavailable; never synthesize it. Preserve existing B2/ChemBERTa nomenclature.

Use the exact four partitions in the new `baseline/{kcat,km}/partition.csv`, whose bytes match historical H1. Target-specific heads; no comparison across different units. kcat/Km are separate; Km remains dataset-native log10 concentration units pending provenance verification. No new label-dependent row filtering. The union has 271 unique sequences across 384 target-specific records; reuse embeddings across targets only when sequence normalization and complete encoder identity are identical.

Primary matched head: Ridge only, via the existing `training.run` model list (`models=['ridge']`). Preserve its existing alpha grid 1/10/100 and validation-MAE selection for all arms; do not choose among head families. Alpha may differ between arms under this common policy. Seeds 11/33/71 and split seed 42 remain fixed; Ridge can give identical seed results, which are not independent datasets. A fixed-alpha=10 secondary comparison would require an explicit existing-pipeline option in a later PR; it must not be falsely advertised as currently supported. Historical adaptive mean/Ridge/MLP scores are contextual, not the matched-head result.

Fit imputation/scaling/EC vocabulary only on the identical training rows, with the same preprocessing rules in every arm. Numeric scaler parameters are arm-specific because embeddings/features differ. Validation selects alpha only. Calibration is untouched during model selection; retain predeclared 90% cluster-maximum residual intervals and report row/cluster coverage and width, without adjusting quantiles against development evaluation. Current Km undercoverage is a known limitation; a better calibration experiment must be separately predeclared, with adequate calibration groups and split replication.

## Evaluation and data gates

Separate log10 MAE, RMSE, R², Spearman, empirical interval coverage/width, encoder extraction time, head fitting/inference time, process CPU-seconds, RAM/VRAM and throughput by length. Paired evaluation differences B1-E−B0 and B1-E−B1 use identical evaluation rows and the existing 500-draw cluster bootstrap; report intervals conditional on this split. Multiple seeds do not measure split uncertainty. Record EC/scaffold overlap. Supplied cluster ID disjointness does not establish an independently audited 40% sequence-identity threshold.

Before broader interpretation: source release/hash and license lineage, assay/reference/mutant provenance, kcat and Km transforms/units, sequence/pair duplicates, independent sequence similarity and pretraining-overlap audit. Unknown values remain unknown. Official test cannot guide encoder, pooling, head, calibration or cohort selection. This already-inspected development evaluation is not a new blind holdout.

## Integration boundary and immutable cache

Reuse FeatureBuilder and the existing training/predict pipeline. The current H1 coordinator only supports B0/B1 and limits rows/length to 256; B1-E or 1024/target is not an existing CLI mode. A separately reviewed adapter/coordinator extension is needed, preserving the current 18-file staged H1 proposal. Do not add another optimizer loop.

The successful one-record API check fed the cached 2880-vector to the unchanged FeatureBuilder; it did not train an EnzGFM head. A new feature dimension requires new train-fitted preprocessing/head artifacts, not reuse of an old ESM checkpoint. Proposed encoder contract: `encode(sequences) -> vectors + metadata`. Bridge the isolated Python 3.10 encoder to the main Python 3.13 pipeline with numeric arrays and checksums; no unknown pickle checkpoints. Retain frozen/eval/inference_mode, local-only files, explicit old-to-new config/class mapping, CPU Mamba override, no truncation, residue masks and versioned pooling. Cache identity covers shard/config/tokenizer/code SHA-256, sequence list/normalization, dtype/backend, runtime and pooling version. Shape/identity mismatch is an error, not fallback.

Current 271-vector payload is 3,121,920 bytes, excluding metadata. Up to 2048 distinct 2880-FP32 vectors would be 23,592,960 bytes. Cache does not improve scientific accuracy: it avoids repeated encoder work across heads/seeds/targets, and allows head-only workers to avoid loading ~6 GiB of weights. For S identical repeated encoding passes the ideal encoder compute saving is (S−1)/S; actual end-to-end savings and cache-read time are not measured.

## Cost and next authorization boundary

The only real EnzGFM observation is 19 residues / 21 tokens: forward 1.132 s, load+weight validation 9.070 s, total child lifetime 12.522 s, peak RSS 6339.5 MiB. It is not a throughput curve or full-model reproducibility test. No full forward was retried.

Linear token-count sensitivity gives ~46.8 minutes for 271 deduplicated sequences (52,123 tokens), ~66.7 minutes without cross-target deduplication. A hypothetical 2048 sequences at mean 194 tokens gives ~5.95 hours. These are uncertain estimates from one short sequence; attention, Python SSM loops, cache/allocator behavior and long-sequence memory may change them substantially. They exceed the original H1 30-minute planning budget and do not authorize execution.

Next minimal separately authorized experiment: at most one sequence each near 64/128/256 residues selected label-free from the existing development training data, same 2-thread/12-GiB-RSS/16-GiB-AS/300-s guards; no regression fitting or evaluation-driven decisions. Use results to freeze an encoding budget. If full-cohort cost is rejected, predeclare a label-free hash-selected paired subset from each existing partition, use the same subset for every arm, document cluster counts and limited calibration precision, and make no population-level claims. Do not choose the subset using errors.

Default recommendation: keep ESM-2 8M as the inexpensive development baseline; test EnzGFM as an ablation only after a bounded length profile and data audit. A smaller EnzGFM checkpoint needs its own availability/license/architecture/resource preflight; its filename alone does not establish a lower total cost. Small target-specific Ridge/MLP or bounded fusion projection heads are reasonable future trainable components after matched evidence; foundation fine-tuning is not justified here.
