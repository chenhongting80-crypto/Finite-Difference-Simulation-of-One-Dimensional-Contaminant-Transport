# Parameter cross-check

The course report and the submitted MATLAB script contain several differences that should be resolved before the model is used for quantitative interpretation.

| Item | Report | Submitted MATLAB script | Repository treatment |
|---|---:|---:|---|
| Method 1 cell density | `9.2 × 10^9 cells/L` | `2 × 10^9 cells/L` | Refactored script preserves `2 × 10^9` |
| Method 2 hydraulic conductivity | `46.4 m/day` | `17 m/day` | Refactored script preserves `17` |
| Method 1 spatial plot times | 10, 50, 150, 200, 250, 300 days | 10, 50, 200, 600 days | Refactored script exports report plot times |
| Method 2 simulation horizon | Results discussed through 620 days | 5000 days | Refactored script runs 700 days because all reported plots occur by day 620 |
| Spatial step | Not independently resolved | 0.03 m | Refactored scripts use `5/167 m` (approximately 0.02994 m) so every zone endpoint and interface lies on one uniform grid |

The archived script is included unchanged in `archive/Term_Project_Bioremediation_original.m`.

Before publication, the authors should confirm which values generated the submitted figures. Once confirmed, update both the scripts and this document.
