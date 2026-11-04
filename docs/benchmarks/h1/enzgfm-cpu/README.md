# H1-E4: controlled CPU execution evidence

Engineering status: **READY_FOR_MATCHED_BENCHMARK**, development research only.
Independent scientific validation remains **BLOCKED_BY_DATA**. This records a
technical smoke test; EnzGFM has not been evaluated for kcat/Km accuracy.

The unchanged H1 B0/B1 coordinator was rerun on the same 192 sequences per target,
splits and seeds 11/33/71. Both partition CSVs and all four ensemble prediction
CSVs matched historical H1 byte for byte; reported metrics matched within 1e-12.
Fitted preprocessing matched in all 12 cases; serialized builders differ only in
new output/data/absent-geometry paths. All 144 historical output files passed
their original manifest checks. Official test access and outgoing network
connections were blocked; neither was attempted.

| Target / arm | log10 MAE | RMSE | R² | Spearman | 90% coverage / width |
|---|---:|---:|---:|---:|---:|
| kcat / B0 | 1.1039 | 1.3894 | -0.4601 | -0.1253 | 97.3% / 6.0385 |
| kcat / B1 | 0.9324 | 1.2169 | -0.1199 | -0.0410 | 100.0% / 7.3380 |
| Km / B0 | 0.9791 | 1.3010 | 0.0701 | 0.5604 | 69.0% / 3.2660 |
| Km / B1 | 1.0026 | 1.3312 | 0.0265 | 0.3452 | 73.8% / 3.6618 |

Baseline rerun: 48.69 seconds, peak process RSS 792.4 MiB, two CPU threads.
Baseline CPU-seconds were not measured separately. Existing fast periodic
geometry/trajectory/Glu43/clash tests: 18 passed. No MD/QM/MM was run; this does
not establish physical convergence. The historical adaptive head scores above
are not a matched-head comparison with EnzGFM.

## Actual EnzGFM 1.5B smoke

Code: [DeepBxM/EnzGFM](https://github.com/DeepBxM/EnzGFM/tree/0631099828e1b2aadbfe56751eec5f67a7bfdf0c),
commit 0631099828e1b2aadbfe56751eec5f67a7bfdf0c. Original model/config mathematics
were unchanged. Isolated Python 3.10.14, torch 2.1.2+cpu, transformers 4.40.1,
safetensors 0.4.3; no CUDA/Triton/Mamba extensions or trust_remote_code.

The locally supplied archive matched [Zenodo record 22042585](https://zenodo.org/records/22042585)
size/MD5; ZIP SHA-256:
`004bcaff1e142343107a47d7232383ffc72062488106bad5f771e8dc4da115eb`.
Only config/index and two safetensors shards were extracted, hashed and made
read-only. Their hashes remained unchanged after execution. No training_args.bin
or unknown pickle checkpoint was deserialized.

Exactly one attempt and one full forward were authorized and executed, exit 0.
A label-free SHA rule selected a 19-residue canonical sequence from development
training data. Batch 1, no padding, FP32, two CPU threads, frozen/eval/inference_mode,
use_mamba_kernels=False, eager attention, no cache/all-layer hidden states.
A meta constructor and safetensors CPU storage assignment avoided a second
initialized model copy. All 387 tensors and their finite FP32 values were checked;
386 encoder keys loaded with strict=True/assign=True after explicit model-prefix
removal and exclusion of exactly lm_head.weight. Missing/unexpected/shape mismatch: 0.

| Measurement | Observed |
|---|---:|
| Load and weight validation | 9.070 s |
| One forward and residue pooling | 1.132 s |
| Whole child process | 12.522 s |
| Sampled process CPU time, including imports/load/checks | 19.54 CPU-s |
| Peak process RSS | 6339.5 MiB / 6.19 GiB |
| Hidden / pooled shape | [1,21,2880] / [1,2880] |

Output was finite and residue pooling excluded special tokens. No outgoing
network/DNS occurred. The process used a tested external RSS watchdog (12 GiB),
300-second timeout and supported RLIMIT_AS (16 GiB address space, not RSS).
Small child self-tests confirmed timeout and RSS termination before the full run.
There was no automatic retry. Full-model repeatability and encoder padding
invariance were not measured; tiny repeatability was established separately.

The unchanged FeatureBuilder accepted the cached 2880-vector for one label-free
training record. This is an API check, not an integrated EnzGFM adapter or trained
head. New feature dimensions require new train-fitted preprocessing and heads.

## Scope, reproducibility and next step

[results.json](results.json) publishes aggregate measurements, checksums and
limitations. Local runtime, archive/shards, embeddings, raw rows, checkpoints and
full logs remain excluded from Git. This report is not a standalone full-model
reproduction kit: the single full forward was not repeated, and the guarded
runtime scripts/environment are retained locally. The B0/B1 CLI and its complete
protocol are documented in [the H1 protocol](../PROTOCOL.md).

The current cohort contains 271 unique sequences / 52,123 tokens across targets.
A linear sensitivity estimate from the single short forward is about 47 minutes
for cached, deduplicated encoding; it is not measured full-cohort throughput.
A hypothetical 2048 unique sequences at 194 tokens averages about six hours under
that approximation. Long-sequence cost and peak memory remain unmeasured. Current
H1 coordinator caps cohort size and sequence length at 256; 1024/target is not a
supported configuration. No bulk EnzGFM inference or regression fitting occurred.

[Matched-head proposal](MATCHED_PROTOCOL.md): B0/B1/B1-E, identical partitions,
Ridge-only family with the existing common validation alpha search, independent
calibration, paired cluster bootstrap and train-only feature fitting. Keep ESM-2
8M as the inexpensive baseline. Next: separately authorized 64/128/256 length
profile before approving a larger encoding budget or predeclared paired subset.

Data release/license/assay lineage, Km units, independent sequence similarity,
pretraining overlap and calibration remain unresolved. Different supplied cluster
IDs do not prove a 40% identity threshold. Negative kcat R² and Km undercoverage
remain visible. No claim that EnzGFM improves accuracy is supported here.
