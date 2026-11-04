# Дневник вычислительных экспериментов

Время ниже: Europe/Moscow. Полные параметры и UTC — events.jsonl.

| Время | Эксперимент | Событие | Подробности |
| --- | --- | --- | --- |
| 2026-11-01 08:01:52 | registration | registered |  |
| 2026-11-01 08:02:22 | H1-static-sites | started |  |
| 2026-11-01 08:02:50 | H1-static-sites | result |  |
| 2026-11-01 08:02:50 | H1-static-sites | completed | 28.1 s |
| 2026-11-01 08:02:50 | H2-conditions-factorial | started |  |
| 2026-11-01 08:02:50 | H2-conditions-factorial | failed | PDB2PQR failed: see /home/vladmyahlov/Документы/GitHub/EnzymeKinetics-Prediction/enzyme-kinetics/outputs/experiments/active-site-v3/ph-5/protein.stderr.log |
| 2026-11-01 08:03:02 | H1-static-sites | interrupted | Preparation profiling found recursive copying of parent structure; fixed with detached Bio.PDB residue.copy before producing results |
| 2026-11-01 08:03:03 | H1-static-sites | started |  |
| 2026-11-01 08:03:03 | H1-static-sites | result |  |
| 2026-11-01 08:03:03 | H1-static-sites | completed | 0.2 s |
| 2026-11-01 08:03:03 | H2-conditions-factorial | started |  |
| 2026-11-01 08:03:04 | H2-conditions-factorial | failed | PDB2PQR failed: see /home/vladmyahlov/Документы/GitHub/EnzymeKinetics-Prediction/enzyme-kinetics/outputs/experiments/active-site-v3/ph-5/protein.stderr.log |
| 2026-11-01 08:03:40 | H1-static-sites | correction | Prior process had exited with a PDB2PQR hydrogen-connectivity failure before attempted stop. Earlier interrupted entry was a profiling annotation, not an actual cancellation. |
| 2026-11-01 08:03:40 | preparation | protocol-amendment | 8CHO includes partial deposited hydrogens; strip H/D and reassign at stated pH. Protein backbone C-N continuity checked in both chains; no gaps >1.8 A. |
| 2026-11-01 08:03:41 | H1-static-sites | started |  |
| 2026-11-01 08:03:41 | H1-static-sites | result |  |
| 2026-11-01 08:03:41 | H1-static-sites | completed | 0.2 s |
| 2026-11-01 08:03:41 | H2-conditions-factorial | started |  |
| 2026-11-01 08:05:14 | WT-dimer-preparation | started |  |
| 2026-11-01 08:05:40 | H2-conditions-factorial | result |  |
| 2026-11-01 08:05:40 | H2-conditions-factorial | completed | 118.8 s |
| 2026-11-01 08:05:40 | study | static-stages-complete |  |
| 2026-11-01 08:06:13 | WT-dimer-preparation | result |  |
| 2026-11-01 08:06:13 | WT-dimer-preparation | completed | 58.6 s |
| 2026-11-01 08:06:46 | H3-replicated-md | started |  |
| 2026-11-01 08:10:49 | H5-projected-proton-surface | started |  |
| 2026-11-01 08:10:49 | H5-projected-proton-surface | failed | Unproductive transferred analog pose: C4--Asp38 oxygen distance 6.979 A > 4.0 A; do not compute a reaction barrier from this pose |
| 2026-11-01 08:12:01 | H2b-temperature-dielectric | registered |  |
| 2026-11-01 08:12:01 | H2b-temperature-dielectric | started |  |
| 2026-11-01 08:12:25 | H2b-temperature-dielectric | result |  |
| 2026-11-01 08:12:25 | H2b-temperature-dielectric | completed | 24.0 s |
| 2026-11-01 08:15:04 | H5-mechanism-status | result |  |
| 2026-11-01 08:15:49 | H4-conditional-torsions | started |  |
| 2026-11-01 08:16:10 | performance | finding |  |
| 2026-11-01 08:16:10 | performance | implementation-improvement |  |
| 2026-11-01 08:17:26 | H4-conditional-torsions | result |  |
| 2026-11-01 08:17:26 | H4-conditional-torsions | completed | 96.3 s |
| 2026-11-01 08:20:59 | H2c-grid-refinement | registered |  |
| 2026-11-01 08:20:59 | H2c-grid-refinement | started |  |
| 2026-11-01 08:22:47 | H2c-grid-refinement | result |  |
| 2026-11-01 08:22:47 | H2c-grid-refinement | completed | 108.7 s |
| 2026-11-01 08:23:48 | H2d-focused-grid | registered |  |
| 2026-11-01 08:23:48 | H2d-focused-grid | started |  |
| 2026-11-01 08:24:31 | H3-replicated-md | result |  |
| 2026-11-01 08:24:31 | H3-replicated-md | completed | 1064.8 s |
| 2026-11-01 08:24:43 | H2d-focused-grid | result |  |
| 2026-11-01 08:24:43 | H2d-focused-grid | completed | 55.3 s |
| 2026-11-01 08:25:33 | performance-WT-throughput | started |  |
| 2026-11-01 08:25:37 | replica-101-analysis | result |  |
| 2026-11-01 08:25:39 | replica-202-analysis | result |  |
| 2026-11-01 08:25:41 | replica-303-analysis | result |  |
| 2026-11-01 08:25:50 | performance-WT-throughput | result |  |
| 2026-11-01 08:25:50 | performance-WT-throughput | completed | 17.9 s |
| 2026-11-01 08:27:17 | H2e-refined-factorial | registered | Refine all 27 conditions after focused grid passed at central probe |
| 2026-11-01 08:27:17 | H2e-refined-factorial | started |  |
| 2026-11-01 08:31:10 | H2e-refined-factorial | interrupted | Profiling found GIL-bound OpenDX parsing with per-line numpy.prod; stopped driver while native solvers were idle, retaining all DX/logs/results. Replaced numeric parsing with numpy.fromstring; resume validates recovered inputs. |
| 2026-11-01 08:31:10 | H2e-refined-factorial-resume | started |  |
| 2026-11-01 08:33:20 | H1-static-sites | reanalyze | Complete chi1/chi2 atom-name maps for binding-contact residues; original coordinate-based role annotations unchanged |
| 2026-11-01 08:33:20 | H1-static-sites | result |  |
| 2026-11-01 08:33:43 | H2e-refined-factorial | result |  |
| 2026-11-01 08:33:43 | H2e-refined-factorial-resume | completed | 153.2 s |
| 2026-11-01 08:36:58 | H5-alternative-product-template | registered |  |
| 2026-11-01 08:36:58 | H5-alternative-product-template | result |  |
| 2026-11-01 08:37:22 | H5-projected-proton-surface | started |  |
| 2026-11-01 08:38:04 | H5-projected-proton-surface | result |  |
| 2026-11-01 08:38:04 | H5-projected-proton-surface | completed | 42.0 s |
| 2026-11-01 08:39:59 | H5-projected-proton-surface | metadata-correction | Replaced stale hardcoded saturated-analog text by actual ASD product template; energies and gradients unchanged, prior result retained in history |
| 2026-11-01 08:39:59 | H5-projected-proton-surface | result |  |
| 2026-11-01 08:41:00 | H1-static-sites | reanalyze | Add directly measured backbone phi/psi and peptide omega, including cis-Pro39; no coordinate changes |
| 2026-11-01 08:41:00 | H1-static-sites | result |  |
| 2026-11-01 08:43:17 | environment | post-run-record |  |
