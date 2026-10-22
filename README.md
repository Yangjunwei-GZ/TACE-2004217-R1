# Manuscript 2004217: portable frozen-input reproduction

Version **v1.0.0** for manuscript 2004217. Repository: [https://github.com/Yangjunwei-GZ/TACE-2004217-R1](https://github.com/Yangjunwei-GZ/TACE-2004217-R1). Use the [fixed version release](https://github.com/Yangjunwei-GZ/TACE-2004217-R1/releases/tag/v1.0.0) and its `TACE_2004217_R1_Code_v1.0.0.zip` asset rather than a moving `latest` link. Cite the release URL and version; no code DOI is claimed.

The author explicitly authorised public distribution of this reviewed candidate on 22 October 2026, replacing the earlier local-only decision for this payload. The original candidate was formed and tested before publication. Its exact ZIP SHA256 is `53a580794817a4226741fcb05aed6e0007d34b59e5e9dc2025578cb204719e3d`; the new release archive has its own checksum. `RELEASE_ORIGIN.json` records this distinction. Actual commit, release-asset hash and anonymous-access verification are recorded outside the committed payload to avoid circular identities. See `CITATION.cff` and the release notes.

The preserved scientific run has 118 patients (83 deaths, 35 censored), 202 active ROI records, 107 Original features per ROI, 17 final masks, 21 model specifications and six paired comparisons. Its ten point repetitions and 500 bootstrap tasks each used five outer folds. The primary clinical + diameter + mean-intensity comparison has delta C-index -0.033737024221453304, 95% interval [-0.09207584210690886, 0.05457958159472272]. Added predictive value was not established; a zero-crossing interval does not establish absence of effect. This is prognosis after TACE, not treatment-benefit prediction.

## Contents and scope

This public version includes the reviewed de-identified patient-level derived data, including survival outcomes and feature matrices, frozen lists and mappings, 17 final masks, reference results and checkpoint identities. It excludes original CT, private physician records, internal revision materials and the complete local audit archive.

* `code/data/analysis_matrix.csv`: portable 118-patient frozen matrix. CT and mask paths are relocated; every non-path field, public ID and row order was independently checked against the original matrix. Original matrix SHA256: `75ceb3fdd9f77831c2a2845b20e409054928bfd6885792a4a3111b33a824b183`.
* `frozen_inputs/`: 126-candidate disposition, 230-source-record image/alias roster, 118-patient/representative roster, source filenames and hashes, exact final-mask application map, and predeclared numerical comparison protocol.
* `final_masks/`: all 17 adopted masks, unchanged. The entry applies these files before extraction; it never substitutes the source mask for a packaged final mask.
* `RELEASE_INPUTS.json` and `code/release_identity.py`: reviewed release-input hashes and an independently fixed manifest fingerprint. A missing, modified or mismatched release fails. The legacy audit CSV differs from its original only in eight source-path cells; `frozen_inputs/PATH_TRANSFORM_PROVENANCE.json` records both identities and the unchanged non-path projection.
* `PUBLIC_RUN_EVIDENCE/` and `current_result_reference/`: identities and numerical reference results of the completed scientific run. **This public release does not contain the 510 original compressed checkpoint files.** Those bytes are in the separate complete local `Code_and_Reproducibility.zip`, under `code/reporting_inputs/04_validation/checkpoints`.

This public release contains code and inputs to generate new checkpoints with new input fingerprints. Original checkpoint bindings are never rewritten to claim they used the portable path-transformed matrix. No raw CT, private review form, signature or author-machine source path is distributed here.

## Tested environment and installation

Actual acceptance used Windows x64, CPython **3.10.11**, and `code/environment-python310-executed.lock.txt`. Use that exact interpreter for extraction or training; the entry checks the installed lock. Other operating systems have not been tested. The supplied `pyradiomics_fix-3.0.1.dev4-cp310-cp310-win_amd64.whl` is for CPython 3.10 on Windows x64 only, not Linux, macOS, ARM or other Python ABIs. Its third-party licence is retained. The radiomics import reports `0+unknown`; the distribution version, wheel and lock identify the actual build.

Example PowerShell setup after unpacking, substituting your actual CPython 3.10.11 path:

```powershell
& 'C:\Python31011\python.exe' -m venv '.venv'
& '.\.venv\Scripts\python.exe' -m pip install '.\code\wheels\pyradiomics_fix-3.0.1.dev4-cp310-cp310-win_amd64.whl'
& '.\.venv\Scripts\python.exe' -m pip install -r '.\code\environment-python310-executed.lock.txt'
& '.\.venv\Scripts\python.exe' -m pip check
```

Install R and the `survival` package for `ph_global.R`; pass the Rscript executable with `--rscript` when training. The current host capability check found R 4.5.1 and survival 3.8-3 (locale warnings were recorded). This round did not refit R diagnostics or any scientific models. R is not needed for matrix identity preflight or raw feature-only extraction. Python lock installation/network availability on another computer is not claimed to have been tested.

## A. Frozen matrix route (no CT required)

Use a new empty output directory outside this package. Commands also work from another working directory when the script path is absolute; quote paths containing spaces.

```powershell
python 'code/run_revision.py' --from-matrix --project-root 'new matrix input check'
```

The default verifies and materialises frozen inputs only; it does not fit models. `--preflight-only` performs the earlier identity check only. To explicitly request a new full scientific run:

```powershell
python 'code/run_revision.py' --from-matrix --project-root 'new matrix model run' --train --workers 2 --rscript 'C:/Program Files/R/R-4.5.1/bin/Rscript.exe'
```

Training executes the unchanged training-only feature ranking/filtering, imputation, standardisation, inner tuning, model fitting and full validation protocol. This route does not re-extract image features. The output is a new run, not a byte-identical copy of the archived checkpoint identities. Input mutation, task errors, incomplete models/folds/repeats or paired comparisons cause a nonzero exit and prevent COMPLETE. The original seeds and planned/attempted/successful/failed task records are retained. No unsuccessful resample may be selectively dropped.

## B. Official CT route to frozen analysis inputs

Obtain original CT NIfTI files and source tumour masks from [WAW-TACE, Zenodo record 12741586](https://zenodo.org/records/12741586), subject to its data terms. That data DOI is **not** this code's DOI. Extract the official archives into a source directory. Only the active files enumerated by `frozen_inputs/SOURCE_FILE_MAP.csv` are read; no download is performed by the entry. Their public filenames and SHA256 hashes determine identity. For example:

```text
official source/
  CT/101/101_1_scan.nii.gz
  masks/101/101_1_0_tumor_seg.nrrd
  masks/101/101_1_1_tumor_seg.nrrd
  ... other active filenames in SOURCE_FILE_MAP.csv
```

Subdirectory names may differ: the resolver finds the listed filenames recursively and requires their expected hashes. Identical duplicate cached files are recorded and deterministically resolved; incorrect bytes cannot satisfy a mapping. Clinical/phase inputs needed to reproduce the frozen analysis are packaged. The 126/230 recorded decisions are applied directly, including exclusions and aliases; the entry does not rerun historical automatic eligibility rules. P112's frozen dispositions, union-volume calculation and source clinical lesion count remain unchanged.

```powershell
python 'code/run_revision.py' --source-root 'D:/official source' --project-root 'new raw preflight' --preflight-only
python 'code/run_revision.py' --source-root 'D:/official source' --project-root 'new raw extraction' --extract-only --workers 2
```

The second command uses all 202 active masks with the exact 17 replacements and freshly extracts all 107 features. It recreates the patient matrix using the unchanged physical-volume/union and representative rules, then verifies exact IDs, row order, clinical counts, representative lesions, complete finite features and numerical agreement. It does not train. The tolerance was fixed before examining extraction differences: absolute 1e-8 plus relative 1e-10, with relative-error reporting denominator max(abs(reference), 1e-12). It is not adjusted to make results pass. The complete comparison is `INPUT_REPRODUCTION_COMPARISON.json`.

To request extraction followed by full training, use another new empty directory and replace `--extract-only` with `--train --rscript PATH_TO_RSCRIPT`. Training starts only after extraction and frozen-input comparison pass. There is no bypass for missing or failed ROIs and no silent cohort reduction. The default never starts training automatically.

## Outputs, failures and descriptive stability

Runs contain `INPUT_PREFLIGHT.json`, input bindings, source-file resolution (raw route), `ACTUAL_ENVIRONMENT.json`, extraction/task audits and `ACTUAL_EXECUTION_LOG.json` as applicable. Stage failures exit nonzero and are recorded in `EXECUTION_FAILURE.json` once an output directory exists; earlier invalid-package errors are printed to stderr and also exit nonzero. Validation creates its own structured task/failure/integrity records in `revision_R1_20261014/04_validation`. Completion requires actual complete checkpoint contents and matching matrix/configuration identities, not only file counts.

The public entry generates numerical results and, for a raw run, image-dependent figures. It does not run a historical Word prose generator. Current clean/marked Word documents and exact material snapshots are in the full local revision package. Never use old report prose to overwrite them.

After a new complete model run, selection recurrence can be independently matched to its newly generated checkpoints:

```powershell
python 'code/summarize_selection_stability.py' --validation 'new matrix model run/revision_R1_20261014/04_validation' --output 'new stability summary'
```

This command requires actual checkpoints plus `feature_recurrence_records.csv`, not their hash list. Alternatively use the original checkpoint directory in the full local archive. No private author path is required. Selection means abs(fitted coefficient) > 1e-8; overlapping training fits are dependent, and joint selection is not reproduction of an identical complete model. The post-review-frozen 0.01 selection window is unchanged and is not represented as prospectively prespecified.

## Acceptance and limits

Run `python code/tests/run_release_acceptance.py --report NEW_TEST_REPORT.json` from an unpacked copy. The acceptance report separately identifies synthetic fault injection, real-input preflight, feature re-extraction and scientific training. A software PASS is not new medical evidence or 510 completed scientific runs. See `ACCEPTANCE_SUMMARY.json` for the final actual scope. The input-reproduction acceptance preceded publication; this release operation performs neither extraction nor model training.

The physician completed the later 57 QC assessments and work on the 17 adopted masks, including two contour changes, as clarified by the author on 21 October. Computational resampling/saving, extraction, statistics, audit and code/text assistance retain their computational attribution. The correction does not supply an otherwise unknown later QC date, specific blinding or independence. No new signature, ICC or repeat-reader experiment is claimed.

Code is MIT licensed; source and derived data remain subject to `DATA_ATTRIBUTION.md`. MIT does not relicense the source data or third-party binaries. The public scope is limited to this versioned package; the complete local archive and internal submission materials are not part of the public upload. Journal-portal access to submission files remains unverified. Public access does not resolve missing DICOM metadata, feature-selection instability, absent external validation or secondary metrics without confidence intervals.
