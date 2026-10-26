# Changelog

## v1.1.0 — authorised release preparation

The author has authorised publication of this dedicated package to Yangjunwei-GZ/TACE-2004217-R1. Code/input preparation is complete; current primary and secondary result attachment and independent numerical checks are complete; final sealed-archive checks remain separate. No actual new remote tag or anonymous-access success is claimed yet.

* The frozen input changes from the prior118-patient/202-record cohort to120 patients,85 deaths and210 active ROI records. Six of126 candidate patients are excluded. Twenty of230 source records are unused: eight from excluded patients and12 from retained patients.
* Eight original source masks for patients83/94 are adopted by author instruction, without new contour edits, merging or aliases. This is research-input adoption, not new independent physician boundary or clinical lesion-identity confirmation. Source clinical tumour counts are unchanged. P112's previously corrected union volume is not subtracted again.
* The original17 adopted masks are retained byte for byte. All210 feature rows are reused only after CT/mask/settings/cache identity verification. A new full model-development and validation run has completed for the changed cohort; no earlier model or checkpoint is passed off as new.
* Portable matrix and official-CT routes bind the new lists, representative records and input fingerprints. The raw route uses exact source files and packaged final masks, can re-extract all107 features, restores frozen patient order and stops on mismatch or failed extraction.
* New checkpoints store fitted baseline and monthly survival evidence. The supplementary module consumes those newly generated objects, preserves statistical NA counts and reasons, uses the fixed1–36-month IBS window, and matches point0 to500 patient-bootstrap intervals. Its actual current numerical output and NA denominators are included after completion.
* Code, source data and third-party licences remain separate. Private physician forms, original CT, full local archives and checkpoint bodies are excluded. Public checkpoint identities will come only from this new verified run.

## v1.0.0 — preserved historical release

The existing published tag, commit and archive are unchanged. Its118-patient numerical results and old checkpoint identities are not current v1.1.0 results. Historical recovery metadata and the old cohort-specific Figure1 helper are explicitly isolated under `historical_reference`.

## Verified v1.1.0 result summary

The verified current run includes 120 patients, 85 recorded deaths and 210 active ROI records. The primary paired comparison adds training-developed radiomics to clinical variables, maximum3D diameter and mean CT attenuation: delta C-index -0.011579818031430933,95% percentile CI [-0.08051608944641396, 0.03963575576919167]. Actual new validation contains 10 point repetitions and 500 patient-bootstrap repetitions, with five outer folds each. Additional predictive value was not established; a zero-crossing interval is not proof of no effect. Current secondary estimates, common evaluable bootstrap counts and statistical NA reasons are supplied in secondary_uncertainty/reference_results. They do not replace the primary paired comparison.

All 13 secondary endpoint intervals were estimable under the frozen rules; 13 paired intervals include zero. Twelve-month AUC uses 480/500 jointly evaluable replicates (20 NA), and 36-month grouped calibration errors use 489/500 (11 NA). The current per-model endpoint records retain 64 NA rows. Statistics reuse the actual newly fitted baselines and probabilities; no regression coefficient is refitted in this supplementary stage. The archived first primary attempt contains ten audit-gate failures, and a separate secondary CLI startup import failure is preserved. Success of the completed attempts is not an overall zero-exception claim.
