# Reproducibility notes

## What is preserved

- The archived MATLAB file is copied without modification.
- The refactored scripts preserve the same governing equations, zone layouts, finite-difference structure, boundary-condition approach, and submitted-script parameter values.
- Report-generated figures are included as reference outputs.

## Refactoring changes

- The two remediation strategies are separated into independent functions.
- Only two time layers are retained during simulation, rather than full space-by-time matrices.
- Selected time series and spatial snapshots are recorded for export.
- Figures and CSV outputs are generated automatically.
- The Method 2 refactored run ends at 700 days because the report presents results only through day 620.
- Source mass is explicitly set to zero when the 0.1% stopping threshold is reached.
- Concentration values below `1e-12 mg/L` are treated as numerical zero for CSV reporting only.

## Verification status

The refactored MATLAB code was prepared from the submitted report and script, but MATLAB was not available in the packaging environment. The scripts should therefore be run locally and checked against the included report figures before the repository is presented as fully reproducible.

## Known source inconsistencies

See `config/parameter_crosscheck.md` for historical report/script differences.

## Authorship

The submitted report states that MATLAB coding was shared among all three group members. Public posting should be approved by all contributors.
