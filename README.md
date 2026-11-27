# RIFT-Map: reproduction material

Material for the peer review of the manuscript

> **RIFT-Map: Metric underwater mapping from monocular images using flexible reference constraints for robotic simulation**
> Jiaqi Wang, Jiaqing Wang, Xinlei Shao, Zheng Zeng, Katsunori Mizuno
> Manuscript under review, Advanced Engineering Informatics

Project page: https://panzerjagerwang.github.io/RIFT-Map/

This repository provides the material for checking the reported numbers: the measurements behind Tables 3 to 6 and a
script that recomputes the tables from them. The manuscript defines the method, the rows of the tables and the
experimental protocol. They are not repeated here, and the files label rows and arms by their order in the manuscript
(`row01` is the first row of Table 4, `tape1` and `arm1` the first tape and the first arm of Table 6).

## How to run

The script needs Python 3 with NumPy and SciPy and runs in a few seconds. It was tested with Python 3.12,
NumPy 1.26.4 and SciPy 1.11.4 and 1.13.1.

```bash
git clone https://github.com/panzerjagerWang/RIFT-Map.git
cd RIFT-Map/reproduce
pip install -r requirements.txt
python reproduce_tables.py
```

It prints every recomputed entry beside the value printed in the manuscript and ends with a summary.
With `--quiet` it prints the summary alone:

```text
table3    4 of 4 entries equal the manuscript
table4    11 of 11 entries equal the manuscript
text_5_3  5 of 5 entries equal the manuscript
table5    7 of 7 entries equal the manuscript
table6    9 of 9 entries equal the manuscript
text_5_5  3 of 3 entries equal the manuscript
ALL VALUES REPRODUCED
```

The exit status is 0 when every value equals the manuscript and 1 otherwise. The bootstrap uses a fixed seed, so the
intervals are reproduced to the printed digit.

## What is recomputed

| In the manuscript | Recomputed quantities | Input in `reproduce/data/` |
|---|---|---|
| Table 3: Registration of the middle-camera frames | registered frames, number of connected models and size of the largest model, from the model sizes of each run | `table3_registration.csv` |
| Table 4: Discrepancy of the 31 withheld intervals | MAE, RMSE, 95% bootstrap interval of the RMSE (20,000 resamples), largest discrepancy, one-sided Wilcoxon signed-rank test and the number of intervals with a smaller absolute discrepancy, from the signed discrepancy of every interval under each of the eleven rows | `table4_withheld_intervals.csv` |
| Table 5: Agreement of the renderings with six withheld photographs | image coverage and structural similarity per frame and surface, their means and the share of tape points inside the tape mask | `table5_rendering_agreement.csv` |
| Table 6: Simulation results over 100 seeds | mean and standard deviation of the scale error, camera error, MAE and RMSE of every tape and arm, from the result of every seed | `table6_simulation.csv` |
| Sections 5.3 and 5.5: Numbers quoted in the text | the ranges and counts that the text states for single rows, intervals and seeds (eight entries) | the files of Tables 4 and 6 |

## Files

```text
reproduce/
  reproduce_tables.py      recomputes the tables and compares them with the manuscript
  expected_tables.json     the values printed in the manuscript
  requirements.txt
  data/                    the four input files listed above
docs/                      the project page
```

## Not part of this release

The reconstruction software, the survey images and the reconstructed scene are not part of this release. Reviewers
who need them for their assessment can request them through the handling editor.

© 2026 the authors. All rights reserved. The material is provided for the peer review of the manuscript.
