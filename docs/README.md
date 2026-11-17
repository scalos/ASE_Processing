# Documentation

| page | read it for |
|---|---|
| [pipeline.md](pipeline.md) | how the pipeline works: diagrams of preprocessing and fitting, frame / regime layout, signal model, optional corrections, quality criteria |
| [settings.md](settings.md) | every config key with its default and meaning, and how to override it |
| [outputs.md](outputs.md) | every file written by preprocessing, a fit run, painting and group statistics, including the ROI table columns |
| [decisions.md](decisions.md) | method decisions (what changed, when, and why) and open questions |
| [quality-mask-plan.md](quality-mask-plan.md) | the original work order (action plan) for the November 2026 quality-mask build |
| [quality-mask-report.md](quality-mask-report.md) | detailed lab log of the November 2026 quality-mask work and review: all numbers, simulations and checks behind decisions.md |
| `../AnLab_Comp/results/report.md` | comparison with the AnLab pipeline on the same subject |

## Quick start

```bash
# 1. preprocess (convert, header, motion correction, registration, segmentations, quality masks)
bash BashFiles/asePreProcess.sh --singleSubjID <id>          # or all subjects in datalist_ase

# 2. fit (defaults from BashFiles/ase_config.json; override fit options with --fit)
bash BashFiles/aseFit.sh --singleSubjID <id>
bash BashFiles/aseFit.sh --singleSubjID <id> --fit dBCorr=true --fit 'quality={"enabled":true}'

# 3. look at a run
bash BashFiles/aseView.sh <subject>/session_ase/fits/<timestamp>

# 4. group statistics per ROI (subjects from cov_table_ase), painted onto a template
bash BashFiles/aseFit.sh --tag calib1                          # tag the runs to use
bash BashFiles/aseGroupStats.sh --fitTag calib1 --select group=ctrl,mdt --covariates age \
    --paint --templateNames <template seg(s)> --LUTs <painting LUT>
```

Each fit run writes `config_resolved.json`, so every result can be traced to its exact settings.

## Tests

```matlab
addpath(genpath('CoreMatFiles'), genpath('Aux_Mat_Functions'));
runtests('Tests')
```
