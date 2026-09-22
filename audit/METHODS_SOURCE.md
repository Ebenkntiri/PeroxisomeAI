# Repaired workflow schematic provenance

Record: WF-01. Evidence level: L0 (direct software inspection). Status: EXECUTED_CODE_INSPECTION_AND_FIGURE_RENDER; biological evaluation SCRIPTED_NOT_RUN.

Source: `src/repair_model.py`, inspected directly on 22 September 2026. The source hash is recorded in `docs/Workflow_Figure_provenance.json`.

The workflow schematic describes the implemented `validate_sequences`, `checked_splits`, `nested_cv`, `BiProfileEncoder` and `fit_sigmoid` functions. It separates required scientific preflight and subsequent biological analyses from implemented fitting logic. It does not describe the website implementation, claim performance, or replace missing data-dependent result figures.

Rendering: `docs/create_workflow_figure.py`, Python with Pillow 12.3.0. Outputs: editable SVG and 3600 × 2880 pixel PNG tagged 300 dpi. The PNG was visually inspected: labels fit, text is legible, and arrows connect the intended workflow steps. The schematic and caption were embedded into the Markdown review draft. No biological analysis was executed during figure generation.
