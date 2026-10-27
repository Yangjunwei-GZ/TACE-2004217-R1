# Historical compatibility fields and current decision provenance

The eight original masks from P120/P125 were retained without contour changes following a physician decision communicated by the author.

`AUTHOR_INPUT_ADOPTION.csv` retains its historical filename for compatibility. In its eight rows, `decision_kind=PHYSICIAN_DECISION_RELAYED_BY_AUTHOR`, `decision_maker_role=physician`, and `communication_source_role=author` distinguish the decision-maker from the person who communicated the decision. The physician's actual decision date, signature, decision-specific blinding and independence were not provided. None is inferred from the author's communication or the correction date.

* `decision_date=2026-10-26` is preserved as the existing communication/adoption record date, not the physician's actual decision date. `decision_date_scope` and `author_communication_record_date` make that scope explicit. `physician_decision_date=not_provided` remains separate.
* `attribution_correction_date=2026-10-27` dates this factual attribution correction only. It is not a clinical review date.
* `physician_retention_decision_relayed=True` records the clarified retention decision. The legacy `physician_final_decision_added=False` is retained as a compatibility flag for the absence of an additional direct physician record; it must not be interpreted as "the physician did not decide". Its explicit scope is `legacy_direct_record_flag_not_decision_maker_indicator`. No signature, new independent review or complete medical-record claim is created by this correction.
* P120L03/L04 remain distinct original segmentation records. Retaining both records is not evidence that they are two independent clinical lesions. No merge, alias or clinical-count change was performed.

`code/reporting_inputs/01_data_audit/phase_analysis_eligibility.csv` retains the original schema required by the analysis entry. Its columns have different time scopes:

* `correctly_paired_arterial_masks` is the historical pre-correction pairing/eligibility count. For patients83/94, the preserved value0 reflects the earlier unresolved-phase classification. It is not a current CT/mask geometry failure or a count of current usable masks. In presentation tables this field can be labelled `historical_pre_correction_pair_count`; the source compatibility column is preserved.
* `record_scope` retains inherited historical classification wording, including the21October attribution-correction text. It is not the date or authority of the later retention decision.
* `phase_analysis_eligible=True` and `final_active_unique_arterial_masks=5` and3 remain the current patient83/94 input values. The prior23October author-provided physician phase correction is unchanged. The later use of eight original masks follows the physician decision communicated by the author; its actual physician decision date is not supplied.
* Current explanatory attribution/reason fields describe the clarified physician decision and the author's communication. Historical record/correction dates do not fill unknown clinical review dates, independence or blinding. The previous57QC items and17final masks retain their physician attribution and original evidence scope.

`ATTRIBUTION_CORRECTION_20261027.json` records the affected fields and their previous public values. `PUBLIC_AUXILIARY_METADATA_PROVENANCE.json` preserves the v1.1.0 compatibility-change history, including its previous interpretation and historical hashes; that history is superseded for current attribution and is not rewritten as if the clarification had always been known. `FROZEN_INPUT_ORIGIN.json` keeps the original run status and identity, with a separate current attribution note.

The authoritative scientific identities remain the FINAL patient/image/representative rosters and actual analysis matrix. Clinical tumour count remains4 for each of patients83/94, independent of their5/3 segmentation records. Source/final masks, all computational roster fields and all frozen scientific numeric cells are unchanged. Metadata-file hashes may change; original run/checkpoint hashes are preserved and do not claim a new scientific run.
