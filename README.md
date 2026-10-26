# TACE manuscript 2004217: new author-input reanalysis candidate

**Author-authorised v1.1.0 release preparation.** The author has explicitly authorised publication of this dedicated payload in [Yangjunwei-GZ/TACE-2004217-R1](https://github.com/Yangjunwei-GZ/TACE-2004217-R1), replacing the prior local-retention restriction for this publication workflow. Code, frozen inputs and reproduction-route preparation are complete. Current-run primary and secondary result attachment is complete; final archive and remote-publication verification are separately recorded. No new tag, commit, successful public download or anonymous-access verification is claimed at this preparation stage; actual publication evidence will be recorded separately after release. The earlier v1.0.0 tag and archive remain unchanged. This is not a submission-readiness certification.

## Scope and evidence

The frozen author-specified inputs contain 120 patients, 85 recorded deaths and 210 active ROI records from 126 candidates and 230 source records, with 107 Original radiomics features per ROI. Source clinical tumour count is preserved independently of ROI count. The prior 17 adopted masks are included unchanged. The eight original masks for patients 83/94 are retained from the official source; they are not packaged as new contour revisions. L03/L04 are retained as separate segmentation records, without inferred clinical lesion identity, merging or alias assignment. P112's previously corrected union volume and clinical count are unchanged.

The author adoption instruction determines research inputs; it does not supply a new physician boundary adjudication, independence, missing review dates or unknown blinding. Those evidence limitations remain explicit. Previously confirmed physician QC/contour work retains its physician attribution. Computational assistance includes code, saving/resampling where documented, feature extraction, statistics and audit. No extra ICC or repeat-reader experiment is claimed.

The current-run model protocol uses 21 specifications, ten five-fold point repetitions and 500 five-fold patient bootstrap tasks, all fit within training folds. These planned counts are now matched by actual new-run validation evidence. No successful-repeat-only substitution is permitted. Outcomes and features are source-public de-identified derivatives. This is prognostic research after TACE, not treatment-benefit prediction. Interpretation must use the completed new-run paired results; old 118-patient values cannot be substituted.

The current input freeze reused all 210 ROI feature rows from explicitly verified caches: the previous 202 rows plus eight newly admitted rows, with the same 107-feature definitions and source CT/mask/settings hashes. No new feature extraction was performed merely to construct the package. Because the patient cohort changed, the current scientific analysis starts a new full development and validation run; it does not reuse old coefficients, scores, bootstrap results or checkpoint bodies. The public raw-CT route remains available to freshly re-extract all 210 records and verify the rebuilt matrix.

## Contents

* `code/data/analysis_matrix.csv` and `code/reporting_inputs/02_extraction/arterial_lesion_features.csv`: exact frozen numerical cells and IDs; only CT/mask paths relocated.
* `frozen_inputs`: frozen patients, source ROI records, representatives, public filename/hash mapping, explicit input-adoption scope and numerical comparison protocol.
* `final_masks`: the same 17 adopted masks. Raw route applies these files by SHA256 before extraction.
* `current_result_reference` and `PUBLIC_RUN_EVIDENCE`: verified current-run reference results and checkpoint identities. Neither old nor new checkpoint bodies are included in the public candidate.
* `RELEASE_INPUTS.json`, pinned `code/release_identity.py`, and `VERSION_FILE_SHA256.csv`: package input and file integrity.

Original CT, private physician forms, internal manuscripts, signatures, machine-specific paths and the full local audit archive are excluded. The full local scientific archive is held separately. Portable CSV hashes differ from exact-run CSV hashes; checkpoint hashes and recorded original input fingerprints must never be rewritten to conceal this relocation.

Historical compatibility columns are explained in `frozen_inputs/COLUMN_SCOPE_NOTES.md`. In particular, the preserved pairing-count zeros for patients83/94 describe the earlier phase classification; current eligibility is true and current active-record counts are5/3. These zeros must not be read as current technical failures.

## Environment

Use Windows x64, CPython3.10.11 and the exact `code/environment-python310-executed.lock.txt`. Install the supplied Windows CPython3.10 wheel for pyradiomics-fix before the lock. Third-party licences are retained in code. R with `survival` is required only for complete model-run proportional-hazards diagnostics; pass an actual Rscript executable. Other platforms are not verified.

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install .\code\wheels\pyradiomics_fix-3.0.1.dev4-cp310-cp310-win_amd64.whl
.\.venv\Scripts\python.exe -m pip install -r .\code\environment-python310-executed.lock.txt
.\.venv\Scripts\python.exe -m pip check
```

## Matrix route

Run from any working directory, quoting the actual absolute script path if necessary. Output must be a new empty directory outside this package. Default entry checks/materialises inputs and fits no model.

```powershell
python code/run_revision.py --from-matrix --project-root "new matrix input check" --preflight-only
python code/run_revision.py --from-matrix --project-root "new matrix model run" --train --workers 2 --rscript PATH_TO_RSCRIPT
```

Use `--workers 2` as separate command-line arguments. The explicit training command reruns complete development and verification; resulting checkpoints belong to a new run with its own hashes. No old checkpoint is copied or resumed. COMPLETE requires matching actual matrix/configuration/source identities, all planned tasks/models/folds and paired results, not merely checkpoint file counts.

## Official CT route

Obtain the official WAW-TACE NIfTI CT and original tumour masks from https://zenodo.org/records/12741586. The source data DOI is not a code DOI. Observe `DATA_ATTRIBUTION.md` and source terms. Recursively extracted source folders may be named freely; filenames plus SHA256 in `frozen_inputs/SOURCE_FILE_MAP.csv` identify files. This includes all eight active 83/94 source masks. No network download is performed by this entry.

```powershell
python code/run_revision.py --source-root "official source" --project-root "new raw preflight" --preflight-only
python code/run_revision.py --source-root "official source" --project-root "new raw extraction" --extract-only --workers 2
```

Raw extraction freshly calculates 107 features for all 210 active masks, then checks exact identities, source clinical values, representatives, physical-union volumes and numerical agreement using the frozen tolerance. No case is silently dropped. Feature provenance text correctly distinguishes new extraction from a verified cache; it is the only additional metadata field ignored during raw numerical comparison. Scientific values and identity columns are never ignored. Raw patient row order is explicitly restored from the frozen roster. To train afterwards, use another empty directory and replace `--extract-only` with `--train --rscript PATH_TO_RSCRIPT`.

The source audit folder contains labelled historical recovery/legacy pipeline records for traceability; these are not current model results. Files in `current_result_reference` must all come from the same new scientific run.

## Acceptance, licensing and remaining limits

Run the focused public-input tests documented in `ACCEPTANCE_SUMMARY.json` when available. Input preflight, cached-feature reconstruction, synthetic failures and actual scientific runs are reported separately. A test PASS does not supply medical evidence or establish predictive improvement. Missing DICOM metadata, uncertain medical provenance fields and lack of external validation remain limitations. Secondary metric uncertainty is included only after actual current-run calculation; prior cohort uncertainty tables are not copied.

Code is MIT licensed. Source/derived data and third-party software retain their own licences; MIT does not relicense them. The author has authorised this dedicated GitHub release. Publication, anonymous download verification and submission-system attachment visibility are distinct operations; their success is not inferred from a local package or test report. No journal submission is performed by this code.

## Secondary uncertainty reproduction

The current parameterised module is `secondary_uncertainty/code/run_new_secondary.py`. Follow `secondary_uncertainty/README.md`: freeze the supplement plan from the new matrix and exact runtime, then compute and audit after all new primary checkpoints pass. It reads stored fitted Cox baselines/probabilities and does not refit regression coefficients. This candidate includes no old 118-patient secondary tables; current reference results are attached from actual completed computation.

## Verified current results

The verified current run includes 120 patients, 85 recorded deaths and 210 active ROI records. The primary paired comparison adds training-developed radiomics to clinical variables, maximum3D diameter and mean CT attenuation: delta C-index -0.011579818031430933,95% percentile CI [-0.08051608944641396, 0.03963575576919167]. Actual new validation contains 10 point repetitions and 500 patient-bootstrap repetitions, with five outer folds each. Additional predictive value was not established; a zero-crossing interval is not proof of no effect. Current secondary estimates, common evaluable bootstrap counts and statistical NA reasons are supplied in secondary_uncertainty/reference_results. They do not replace the primary paired comparison.

All 13 secondary endpoint intervals were estimable under the frozen rules; 13 paired intervals include zero. Twelve-month AUC uses 480/500 jointly evaluable replicates (20 NA), and 36-month grouped calibration errors use 489/500 (11 NA). The current per-model endpoint records retain 64 NA rows. Statistics reuse the actual newly fitted baselines and probabilities; no regression coefficient is refitted in this supplementary stage. The archived first primary attempt contains ten audit-gate failures, and a separate secondary CLI startup import failure is preserved. Success of the completed attempts is not an overall zero-exception claim.
