# Historical compatibility fields and current eligibility

`code/reporting_inputs/01_data_audit/phase_analysis_eligibility.csv` retains the original schema required by the analysis entry. Its columns have different time scopes:

* `correctly_paired_arterial_masks` is the historical pre-correction pairing/eligibility count. For patients83/94, the preserved value0 reflects the earlier unresolved-phase classification. It is not a current CT/mask geometry failure or a count of current usable masks. In presentation tables this field can be labelled `historical_pre_correction_pair_count`; the source compatibility column is preserved.
* `record_scope` retains inherited historical classification wording, including the21October attribution-correction text. It is not the actual date or authority of the later input-adoption decision.
* `phase_analysis_eligible=True` and `final_active_unique_arterial_masks=5` and3 are the current patient83/94 input values. They reflect the prior23October author-provided physician phase correction and26October author adoption of the eight original masks. This adoption is not new physician contour/clinical-lesion adjudication.
* `current_medical_work_attribution`, `recorded_eligibility_decision_reason` and `correction_date` in the public compatibility copy are transcribed from the current input freeze. Unknown actual physician review dates, independence and blinding are not filled from these record/correction dates.

The authoritative current scientific identities are the FINAL patient/image/representative rosters and the actual analysis matrix. Clinical tumour count remains4 for each of patients83/94, independent of their5/3 segmentation records. All original source files and frozen scientific numeric cells are unchanged. No historical zero is overwritten merely to make a current report look complete.
