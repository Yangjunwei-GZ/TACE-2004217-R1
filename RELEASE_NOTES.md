TACE manuscript 2004217 — reproducibility release v1.0.0

This fixed version contains analysis code and reviewed de-identified derived inputs for 118 patients, 83 deaths and 202 valid ROI records, including 17 adopted masks. Use `--from-matrix` without CT, or `--source-root` to apply frozen decisions/masks to official CT and reproduce the analysis inputs. Training requires explicit `--train`.

Acceptance before publication: 47 focused software/synthetic tests; actual matrix-route input checks; real extraction of 202 ROIs × 107 features and reconstruction of the 118-patient matrix, with maximum non-path absolute/relative error 0 under the predeclared tolerance. The final candidate was separately unpacked for actual source preflight and numerical comparison, with extraction-kernel/input identity checks. These are distinct from the existing 10 point/500 bootstrap scientific runs. No extraction or model training was rerun for this publication operation.

The primary paired delta C-index is -0.033737024221453304 (95% CI [-0.09207584210690886, 0.05457958159472272]); additional predictive value was not established. The study concerns prognosis after TACE, not comparative treatment benefit. There is no external validation.

Tested platform: Windows x64, CPython 3.10.11, bundled dependency lock and Windows cp310 PyRadiomics wheel; R with survival is needed for PH diagnostics. Other platforms have not been tested. Obtain original CT/masks from WAW-TACE: https://zenodo.org/records/12741586. This source-data DOI is not a code DOI.

This public release includes checkpoint identity/hash records and code/inputs to generate new checkpoints; it does not include the 510 original compressed checkpoint files. Full local audit archives, private physician forms, internal revision materials and original CT are excluded. Code MIT licensing, source-data terms and third-party licences remain separate. Physician medical work and computational/AI assistance retain their documented attribution.

Original candidate ZIP SHA256: 53a580794817a4226741fcb05aed6e0007d34b59e5e9dc2025578cb204719e3d
Formal release ZIP: TACE_2004217_R1_Code_v1.0.0.zip (new identity after documentation changes).
Actual commit, asset checksum and anonymous-download verification are recorded separately; no code DOI is assigned.
