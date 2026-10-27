# Changelog

## v1.1.1 — medical attribution and publication-status correction

The eight original masks from P120/P125 were retained without contour changes following a physician decision communicated by the author.

* The physician is the decision-maker and the author the communication source. Eight retain actions, original mask bytes, absence of merging/aliases and source clinical counts are unchanged. This does not establish independent clinical identity for P120L03/L04 or provide a physician decision date, signature, decision-specific blinding or independence.
* `PHYSICIAN_DECISION_RELAYED_BY_AUTHOR` and explicit role/evidence fields replace the previous current attribution. Compatibility checks preserve all input-integrity gates, distinguish legacy records and do not treat a retention decision as new complete medical documentation.
* Repository status now distinguishes the already published v1.1.0 from the assembly and actual publication evidence of v1.1.1. Exact original top-level files remain under `historical_reference/v1.1.0_pre_correction_metadata`.
* Scientific numerical values, lists, representatives, masks, model settings/parameters, random seeds, completed validation and original checkpoint identities are unchanged. No feature extraction or model training is performed for this patch. Changed metadata/code hashes are separate from the preserved original scientific-run hashes.
* Authorisation to publish the dedicated package remains an author action. Physician work on the prior 57 QC items and17 final masks, AI-assisted work and all historical execution/failure evidence remain disclosed.

Actual v1.1.1 commit, Release asset hash and anonymous-download verification are recorded separately after publication in `PUBLICATION_VERIFICATION_v1.1.1.json`, without inserting a commit's own hash into itself or moving earlier tags.

## v1.1.0 — published corrected-cohort scientific release

Published fixed Release: https://github.com/Yangjunwei-GZ/TACE-2004217-R1/releases/tag/v1.1.0. Commit: `b9f39ac0f3a72b457b044bd2c5d4e611c694ec26`. ZIP SHA256: `412ddb43397690af0a79692dff2b9d74ca6ec031bbd5a5673b657f4e843601fa`. Prior anonymous verification was recorded at 2026-10-27T00:47:58+08:00. The preparation descriptions originally shipped in that archive are historical snapshots, not its actual current publication status.

This release completed a new 120-patient, 85-death, 210-ROI analysis with 21 model specifications, 10 five-fold point repetitions and 500 five-fold patient bootstrap tasks. All 210 feature rows were reused only after source/cache identity checks. It contains verified numerical references, portable matrix/official-CT routes, stored baseline/survival evidence support and checkpoint identity lists. It does not include original CT, private physician forms or checkpoint bodies. The earlier statement that the author medically decided to retain the eight original masks is corrected above; the adoption actions and scientific results remain the same.

## v1.0.0 — preserved historical release

The published tag, commit and archive are unchanged. Its 118-patient numerical results and checkpoint identities are not the current 120-patient results. Historical recovery metadata and the old cohort-specific Figure 1 helper remain isolated under `historical_reference`.

## Verified 120-patient results retained without recomputation




The verified current run includes 120 patients, 85 recorded deaths and 210 active ROI records. The primary paired comparison adds training-developed radiomics to clinical variables, maximum3D diameter and mean CT attenuation: delta C-index -0.011579818031430933,95% percentile CI [-0.08051608944641396, 0.03963575576919167]. Actual new validation contains 10 point repetitions and 500 patient-bootstrap repetitions, with five outer folds each. Additional predictive value was not established; a zero-crossing interval is not proof of no effect. Current secondary estimates, common evaluable bootstrap counts and statistical NA reasons are supplied in secondary_uncertainty/reference_results. They do not replace the primary paired comparison.



All 13 secondary endpoint intervals were estimable under the frozen rules; 13 paired intervals include zero. Twelve-month AUC uses 480/500 jointly evaluable replicates (20 NA), and 36-month grouped calibration errors use 489/500 (11 NA). The current per-model endpoint records retain 64 NA rows. Statistics reuse the actual newly fitted baselines and probabilities; no regression coefficient is refitted in this supplementary stage. The archived first primary attempt contains ten audit-gate failures, and a separate secondary CLI startup import failure is preserved. Success of the completed attempts is not an overall zero-exception claim.

