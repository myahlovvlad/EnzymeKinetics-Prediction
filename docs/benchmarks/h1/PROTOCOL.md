# Benchmark-H1 — frozen encoders and small kinetics heads

Status: development pilot implemented; full independent benchmark NOT validated.
Protocol v1, 2026-11-04. Base main 181c45e, branch benchmark/h1 in a separate
worktree. C1 PR #2 and proposed C2/ADR-0002 are separate work. No MD/QM/MM,
foundation-model fitting, environment installation or weight download is needed.
No historical report is edited. The coordinator calls learning.training.run and
learning.predict.predict; it contains no alternative optimizer/training loop.

## Registered feature comparisons

| Arm | Protein | Ligand | Other inputs |
|---|---|---|---|
| B0 | current FeatureBuilder sequence composition/length | current Morgan radius 2, 256 bits + six RDKit descriptors/status | EC, T/pH/missingness, current titration features; geometry unavailable in this pilot |
| B1 | B0 + frozen ESM-2 8M, residue mean, 320 dimensions | same RDKit as B0 | identical B0 covariates |
| B2 (deferred) | same frozen protein encoder | frozen ChemBERTa replaces RDKit ligand block | same covariates; separate controlled RDKit+ChemBERTa ablation if registered before evaluation |
| B3 (deferred) | same as B2 | same as B2 | validated available target-independent physical sidecars and explicit missingness |

B0 is an H1 run of the existing handcrafted feature builder and heads, not a
rerun of a published v2 score or of the old Progres representation. Progres is
explicitly disabled rather than relabelled ESM-2. Conditions and titration are
shared control features, not new evidence for an enzyme mechanism. Geometry
and physical sidecars have zero coverage here; no missing physics is synthesized.

Only mean/Ridge/KineticsMLP heads fit. Ridge alpha uses the existing validation
search (1,10,100); MLP is the existing small 128/64 head, AdamW, <=20 epochs,
patience 5. Head and checkpoint selection use validation MAE only. Seeds
11/33/71 are fixed; split seed 42 is independent. Equal-weight prediction
ensemble across all three seeds is registered in advance, not chosen using
holdout scores. Deterministic head selections can have identical seed results;
three seeds are not three independent datasets.

## Two levels of evidence

**Current pilot:** only local kcat/km random_trainval.csv files are read. Their
upstream commit and assay-unit lineage are not verified; local CatPred folder
has no .git identity. Eligibility: length 1..256, one row per exact sequence
(first source row), order by SHA-256(sequence), first 192. No label affects this
choice. This deliberately biased short-protein pilot checks implementation,
not representative dataset performance. All arms use the exact same rows.

Provided sequence_40cluster groups split this cohort into fitting pool,
calibration and development evaluation; a further group split creates tuning
validation. Approximate proportions are 45/15/20/20 but cluster sizes determine
actual counts. Original row IDs and partitions are saved. Group labels are
used as supplied; neither their provenance nor true maximum cross-split identity
has been independently established. No claim of <=40% identity is made.

Official random_test.csv is never opened, hashed, repartitioned or used by the
pilot. training --development-only enforces this boundary and rejects
--split-scheme grouped-three-way. Fitting imputer/scaler/EC vocabulary sees train
only. Validation is for selection, calibration for residual intervals, development
evaluation for metrics after selection. Inference receives only an explicit
label/product-free feature-column whitelist. Train-derived aggregate targets,
reaction products and dataset target grouping columns are excluded.

**Full benchmark (blocked pending data/protocol verification):** pin and hash an
upstream release with official independent train/validation/test or create a
new documented external holdout without altering any official test. Previous
v2 reports admit the old test was already inspected. Do not market that test as
a new blind confirmation. Freeze head search, feature arms, preprocessing,
thresholds and calibration policy before one final official/external evaluation.
Do not rank results across CatPred, KinBench or EnzGFM publication tables with
different targets, transforms, units, filters or splits.

## Leakage and provenance gate

The pilot persists pairwise overlap counts for exact sequence, supplied 40%
cluster, canonical enzyme/substrate pair, EC and nonempty Murcko scaffold.
Sequence/cluster/pair overlap blocks the run; EC/scaffold overlap is measured,
not claimed disjoint. Acyclic/invalid molecules do not contribute a nonempty
scaffold and cannot be claimed scaffold-novel. Invalid SMILES retain the existing
explicit indicator and are reported, not silently filtered differently per arm.
Exact sequence uniqueness prevents duplicate pairs across the four partitions.
Repeated assays are not arbitrarily averaged or used to increase sample size.

Before the full benchmark additionally require:

- original assay/reference identifiers, exact and canonical-pair duplicates,
  conflicting repeat values, EC overlap by level and mutant/WT lineage;
- independently audited clustering/alignment identity and scaffold grouping;
- assay pH/temperature/method/organism/construct, missingness and reference overlap;
- kcat s^-1 and each Km concentration unit verified from the source, including
  any existing log transform; never silently rename dataset-native Km units to M;
- train/val/test source checksums, stable row IDs, release and license lineage;
- overlap against protein pretraining sequences and ligand pretraining corpus
  where membership is available. Missing corpus membership remains UNKNOWN.
  Sequence exposure in unsupervised pretraining differs from kinetic-label leakage.

Label-free full-test feature audits may later read an explicit column whitelist;
labels remain sealed during preprocessing/selection. This pilot does not read
that file at all. Supplied cluster labels alone do not prove no data leakage.

## Evaluation and uncertainty

kcat and Km have separate heads, reports and metrics. Errors are log10 MAE/RMSE,
R² and Spearman; constant-vector Spearman is undefined/null. Model spread across
seeds is reported descriptively; it is not an independent-replicate confidence
interval. 500 percentile cluster-bootstrap draws on development evaluation yield
conditional 95% CIs; undefined draws are excluded and counted. These intervals
omit dataset selection, training-set and pretraining uncertainty.

A predeclared 90% calibration radius uses each calibration cluster's maximum
absolute residual, rank ceil((G+1)*.9). If that rank exceeds G, radius is null
rather than an invented finite interval. Report empirical row and all-members
cluster coverage and width. Exchangeability/stationarity are unverified; short
clustered pilot coverage is not a guarantee of calibrated deployment uncertainty.
B1-minus-B0 MAE uses paired cluster resampling on the same evaluation rows.

Report head/inference wall time separately from ESM encoding. RSS is Linux
process-lifetime high-water mark, not incremental per-model memory or GPU memory.
The first B0 seed is retrained and loaded predictions must match exactly before
B1 runs. Official test selection and comparisons to publication scores are
outside this development evidence.

## Resources and reproduction

Pilot config: enzyme-kinetics/configs/benchmark-h1-pilot.json. CPU only, two
threads, batch 2 for ESM, <=192 sequences/target and <=256 residues, three seeds,
small heads. Planning budget: 8 GiB RAM / 30 minutes, not a measured hard cap.
Recorded RSS/runtime must be checked before scaling. No all-dataset encoding,
new GPU job or foundation fine-tuning is authorized by this pilot.

Use an existing environment (no automatic sync/download):

```sh
PYTHONPATH=enzyme-kinetics/src HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 \
  uv run --no-sync python -m enzymekinetics_prediction.learning.benchmark_h1 \
  --config enzyme-kinetics/configs/benchmark-h1-pilot.json \
  --data-dir /path/to/existing/enzyme-kinetics/data \
  --output enzyme-kinetics/outputs/benchmarks/h1-pilot-new
```

In the isolated worktree the actual command used the original repository's
.venv/bin/python with the same PYTHONPATH, without installing a second environment.
The frozen encoder must already exist in Hugging Face cache at the pinned
revision. Offline cache miss blocks B1; never fall back to random weights.
Output must not exist. Incomplete attempts are preserved and a new output
version is required. Models, raw source rows and embeddings stay ignored/local;
portable aggregate results and checksums go into this report directory.

## Future fusion and own weights

Extend the existing FeatureBuilder with explicit ligand-embedding blocks or
label-free sidecars keyed by canonical substrate and model revision. Reuse
train-only preprocessing, existing ConsensusMBPLS block normalization and small
Ridge/MLP heads. Keep substrate-component handling explicit; do not discard
water/co-substrates or choose a different salt policy per arm. B3 requires a
common frozen cohort or a matched-coverage sensitivity analysis for all arms,
not just filtering favorable B3 rows. Cached physical results need method,
settings, source and units validation; no new MD/QM/MM is part of H1.

Own weights worth considering first: compact target-specific fusion heads,
block projections and calibration artifacts, after verified splits/units and
stable ablations. Reuse frozen protein/ligand encoders. GNNs, LoRA/deep tuning
or an enzyme foundation model need a separate resource/data/validation decision;
this pilot cannot justify them. Separate proposed PRs: H1 development boundary
and B0/B1 evidence; verified external benchmark/data provenance; licensed frozen
ChemBERTa integration; matched-coverage physical ablation. C1/C2/Harness remain
separate and no commit/push/merge is performed here.

## Acceptance revision: cache identity and historical v1 evidence

The engineering candidate adds strict H1 embedding identity: immutable model
revision, weight/config/vocab/tokenizer_config/special_tokens_map hashes,
normalized sorted sequence-set hash/count, pooling/normalization policy,
algorithm ID plus exporter source SHA and software versions. Before B1 fit,
training.read_embeddings validates that identity, tensor payload checksum,
sequence membership, float32 finite values and dimension=320. Cache paths are
not identity. Missing/mismatched identity fails rather than producing missing
or random embeddings. No cache-hit inference speedup was claimed for pilot v1:
its exporter always refused existing outputs.

Legacy read_embeddings(path) and trusted preprocessing/checkpoint readers stay
compatible. Strict checks are opt-in and mandatory in the H1 coordinator. Old
v1 sidecars lack the new identity and are deliberately rejected as strict caches;
they remain historical inputs readable through the legacy API. We do not
retroactively attach a new algorithm identity to those files. Preserved v1
manifest/code hashes describe the code used for the published pilot; current
candidate hashes and acceptance checks are recorded separately in evidence.json.
No foundation fitting or new H1 performance experiment was run during acceptance.

## H1.1 — registered next experiment, not executed

### Gates before compute

1. Audit local CSV lineage against a pinned upstream release, documented ETL,
   source-row/reference identifiers and checksums. Verify kcat time units, Km
   concentration units, log base, temperature conventions and repeat assays.
   Check licenses for source data separately from repository code. Unknown unit
   or missing redistribution right is an explicit blocker, not a guessed value.
2. Reuse the v1 acceptance findings as development evidence. Its evaluation
   labels have now been inspected and cannot become a fresh blind confirmation.
   Keep official test sealed. Define a new external holdout only with verified
   lineage/rights and preregister evaluation before inspecting its labels.
3. Freeze cohort, feature arms, transformations, score/interval definitions,
   splits and budget in a new versioned config before model selection. The
   current H1 pilot guard must not be bypassed silently for a larger experiment.

### Controlled comparison and partitions

- Up to 1024 exact-sequence-unique records per target, max 512 residues, selected
  by a label-free documented rule. Report excluded/covered rows by length and
  group rather than calling this the complete CatPred benchmark. Do not increase
  the cohort until units/lineage gates permit interpretation of its targets.
- Independently compute/audit homology grouping on permitted trainval sequences
  with pinned tool/version, identity denominator, alignment coverage, sensitivity
  and exact inputs. Reject unknown/sentinel/inconsistent assignments. Measure
  cross-partition maximum similarity explicitly; cluster labels alone prove no
  40% identity bound. Document treatment of fragments, mutations and constructs.
- Three predetermined outer development folds with nonoverlapping evaluation
  groups; each remaining pool has separate train/validation/calibration groups.
  Save exact row IDs. Training folds overlap, so fold scores are not independent
  experimental replicates. Disjoint experimental cohorts would be required for
  such a claim. Seeds 11/33/71 repeat stochastic heads within each split, not the
  dataset itself. Never pool all seed/fold rows as iid observations.
- Primary matched-head comparison: same mean anchor and Ridge family, fixed
  alpha=10, identical common-block preprocessing and fixed splits. Adding the
  320 frozen coordinates is the sole primary feature change. Secondary Ridge
  search [1,10,100] uses equal validation budgets and is labelled adaptive.
  Any MLP comparison fixes architecture/epochs/stopping and reports parameter
  counts and training budget; no winner is picked from evaluation labels.
- Cache protein vectors once per frozen model/sequence-set/algorithm identity,
  reuse across splits, and verify hashes on reuse. All statistical transforms
  still fit per-fold train only. Assert B0/B1 input rows, assay covariates,
  common-block imputer/scaler and target transformations match per fold.
- Quantify exact-sequence/pair/reference overlap, sequence similarity, EC-level
  and substrate scaffold overlap separately. A homology split is not EC/scaffold
  OOD; separate EC/scaffold challenges need their own frozen group policy.

### Calibration acceptance and diagnostics

First reproduce the finite-sample quantile rule and define the prediction unit:
row coverage and all-members cluster coverage are different estimands. Choose
alpha=.1 beforehand. Fit/select the point predictor using train/validation only,
then freeze it. Use >=100 independent calibration groups per target/fold where
available; absence is a feasibility failure, not permission to reuse evaluation.
Size is a planning minimum, not a theorem guaranteeing calibration.

Report calibration/evaluation residual distributions, bias, group sizes,
length/EC/scaffold/condition coverage and interval widths. Report coverage with
cluster-aware uncertainty, not a bare percentage. Never increase the radius
until the already inspected v1 evaluation set reaches 90%. A single observed
coverage below 90% is not by itself a proof of implementation error; finite
samples, dependence and distribution mismatch require separate diagnosis.

Keep constant-radius group calibration as the reference. A possible predeclared
secondary candidate is residual-scale normalization learned using training OOF
residuals only; calibration labels determine only its final quantile. The scale
model/checkpoint cannot be selected on calibration or evaluation coverage.
CQR is a later candidate requiring its own head/objective contract, not an
unannounced patch to H1. No empirical calibration success is claimed before
validation on the intended prediction unit/population.

Primary foundations: [Angelopoulos & Bates](https://arxiv.org/abs/2107.07511)
and [Romano, Patterson & Candes](https://arxiv.org/abs/1905.03222). Their guarantees
must not be transferred to unknown dependent/OOD clusters without assumptions.
Only HTML abstracts were consulted here; no PDF analysis or new method fit.

### Budget, outputs and stopping rules

Proposed ceiling: CPU 2 threads, 8 GiB RSS, 2 GiB new local artifacts, 30 minutes,
<=1024 records/target, ESM batch 2, no model downloads/GPU/foundation tuning.
These are proposed stopping limits requiring measurement/enforcement in H1.1,
not a claim that current H1 enforces OS memory or wall-clock limits. Reuse the
existing training loop; add a new bounded config and tests in its own PR.

Stop on lineage/unit/license failure, identity/partition mismatch, inadequate
calibration groups or budget overrun. Preserve failed attempts. Artifacts:
resolved config, lineage/unit/license decisions, cluster-tool manifest, frozen
partition map, per-arm/seed metrics and predictions, paired cluster-bootstrap
contrasts, residual/coverage/width diagnostics, runtime/RSS and source hashes.
CIs must state whether they condition on a split/model or include refitting;
percentile resampling is not a causal test of the encoder.

## H1.2 — Small Fusion Model contract, not implemented

Prerequisites: H1.1 data/units/identity gates, a reviewed calibration protocol,
and a separate implementation authorization. No new training pipeline, external
executor or C2 dependency is introduced for this research experiment.

### FeatureBuilder boundary

- Frozen protein block: ESM-2 8M, exact model/tokenizer/sequence/algorithm identity;
  320-vector + availability mask; no gradients and no random-weight fallback.
- Ligand block: current 256-bit Morgan + six RDKit descriptors and validity flag;
  canonical substrate/mixture policy fixed for every arm. Conditions T/pH and
  their missingness remain explicit. Existing EC/sequence/titration controls are
  either shared by all arms or removed in a preregistered matched ablation.
- Optional physics: <=32 validated label-independent sidecar features, explicit
  keys, units, method/settings/source hashes and per-feature missingness. No
  on-the-fly MD/QM/MM; absent physics is unavailable, not an invented zero.
- Fit imputation/scaling and category vocabulary on train only. Missingness
  indicators distinguish numerical imputation from measured physical values.
  Freeze block slices/units/order in preprocessing metadata for predict/reload.

### Small model and fitting contract

Use existing models.py/training.py dispatch: protein 320->32, ligand 263->32,
conditions/masks ->8, optional physics+mask ->8; concatenate available projected
blocks and shared masks, then a small 64->32->1 log10 head. Cap the complete
model at 100k trainable parameters per target; record actual dimensions/count
and fail the cap rather than resizing from evaluation results. The code may
share a class but kcat/Km models and fitted weights are separate: no unreviewed
cross-target label sharing through a jointly trained trunk.

Use the existing train-only target standardization, MSE/AdamW, validation MAE
checkpoint selection, bounded epochs/patience and fixed seeds. Freeze a small
search budget before evaluation. Compare to B0/B1 on the same cohort, partitions,
source versions, target units and measured compute caps. Additional F0/F1 fusion
ablation keeps the same architecture/capacity, masking the protein slot for F0;
its results do not replace the historical B0/B1 rows. Physics comparisons require
identical eligibility or a preregistered matched-coverage sensitivity analysis.

Intervals are post-hoc, target-specific and calibrated after the predictor is
fixed, using untouched calibration groups. Report empirical row/group coverage,
interval width and uncertainty; a point head or seed spread is not a calibrated
uncertainty model. Failed/rejected calibration remains a publishable negative
result and prevents deployment claims.

### Future implementation scope and proof obligations

Extend only existing FeatureBuilder (features.py), model class (models.py),
training dispatch/metadata (training.py) and compatible inference (predict.py),
plus bounded experiment config and synthetic tests. Proposed exact paths must
be approved in that later PR; none of these fusion changes is implemented here.
Preserve old feature/checkpoint versions and add a new explicit model/feature
schema instead of changing how old pickles/checkpoints are interpreted.

Required tests: frozen gradients; label/product exclusion; block dimensions and
masks; train-only transforms; old/new reload agreement; identical matched splits;
optional-sidecar join/units/hash failure; trainable-parameter/compute limits;
selection without calibration/evaluation labels; interval quantile correctness;
no output overwrite and no hidden weight/network access. Rollback removes the
new optional class/dispatch and leaves all old runs/readers untouched.
