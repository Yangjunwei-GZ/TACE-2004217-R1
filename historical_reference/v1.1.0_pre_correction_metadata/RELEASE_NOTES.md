# TACE manuscript2004217 — v1.1.0 authorised release preparation

The author has authorised this dedicated code/data payload for the existing public repository `Yangjunwei-GZ/TACE-2004217-R1`. The v1.0.0 tag will not be moved or overwritten. This record describes preparation; actual commit, fixed Release URL, asset checksum and unauthenticated download verification must be recorded after the real publication operation.

Frozen input:120 patients,85 deaths,210 active ROI records,107 features per record, and the same17 adopted masks. Eight original masks for patients83/94 are obtained through the official WAW-TACE source and are adopted as research inputs by the author; they are not packaged as new contours or described as new physician medical confirmation. Clinical tumour counts remain source values and are separate from ROI-record counts.

Both frozen-matrix and official-CT reproduction routes are provided. For this analysis,210 feature rows were reused after source and settings identity checks. The changed cohort requires entirely new training, tuning and full-flow validation; current reference values are attached from that completed run. Old118-patient results are not substituted. The main comparison and secondary estimates are reported as actually calculated, with all evaluability limitations retained.

The new secondary module retains all planned repetitions and statistical NA reasons. AUC/Brier at12/24/36 months, fixed-window IBS1–36 months and grouped calibration uncertainty use the documented paired estimands. Medical provenance gaps, absent external validation and selection-stability limits are not resolved by code publication.

Use Windows x64,CPython3.10.11 and the supplied lock; R with survival is required for the complete PH-diagnostic route. Obtain original CT and source masks from the WAW-TACE data repository, https://zenodo.org/records/12741586. That data DOI is not a code DOI. No new DOI is requested.

The public package contains code, de-identified public-source derived inputs/outcomes, frozen lists and maps,17 final masks, reference results after verification, and checkpoint identities. It contains no original CT, private doctor forms, full local archive or original/current checkpoint bodies. Licence notices distinguish MIT analysis code, source/derived data terms and third-party software.

Before release sealing, execute final fresh-package acceptance and regenerate exact upload inventory and checksums. Current-run main and secondary reference verification is complete. After actual publication, perform anonymous page and asset-download checks; no success is claimed beforehand.

## Verified current numerical result

The verified current run includes 120 patients, 85 recorded deaths and 210 active ROI records. The primary paired comparison adds training-developed radiomics to clinical variables, maximum3D diameter and mean CT attenuation: delta C-index -0.011579818031430933,95% percentile CI [-0.08051608944641396, 0.03963575576919167]. Actual new validation contains 10 point repetitions and 500 patient-bootstrap repetitions, with five outer folds each. Additional predictive value was not established; a zero-crossing interval is not proof of no effect. Current secondary estimates, common evaluable bootstrap counts and statistical NA reasons are supplied in secondary_uncertainty/reference_results. They do not replace the primary paired comparison.

All 13 secondary endpoint intervals were estimable under the frozen rules; 13 paired intervals include zero. Twelve-month AUC uses 480/500 jointly evaluable replicates (20 NA), and 36-month grouped calibration errors use 489/500 (11 NA). The current per-model endpoint records retain 64 NA rows. Statistics reuse the actual newly fitted baselines and probabilities; no regression coefficient is refitted in this supplementary stage. The archived first primary attempt contains ten audit-gate failures, and a separate secondary CLI startup import failure is preserved. Success of the completed attempts is not an overall zero-exception claim.
