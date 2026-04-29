# AGENTS

## Scope
This repository is a research codebase for resolution/relevance analysis of CA1 and SUB neural recordings.
Keep edits narrow and local. Prefer changing the existing script or notebook that owns the behavior instead of introducing new abstractions.

## Stack
- Python is the main execution path.
- MATLAB is used in `data_extract/` and `data_model/`.
- Julia mirrors some utilities in `lib/`, but Python is the primary workflow.
- Data is stored as JSON and CSV.
- There is no package manifest and no formal automated test suite.

## Working Directories
Many scripts depend on relative paths and should be run from their own directory, not from the repo root.

- `cd data_gen && python ResRel_gen.py`
- `cd data_bool && python BoolOP_ResRel_gen.py`
- `cd data_jitter && python jitter_gen.py`
- Open notebooks from the folder they live in.

Use cheap syntax checks before long runs:

- `python -m py_compile lib/*.py data_gen/*.py data_bool/*.py data_jitter/*.py`

## Repo Map
- [`lib/func_MSR.py`](lib/func_MSR.py): core MSR and spike-train utilities.
- [`lib/func_Info.py`](lib/func_Info.py): spatial and head-direction information metrics.
- [`data_extract/`](data_extract/): MATLAB extraction scripts that produce source JSON.
- [`data_gen/`](data_gen/): baseline per-rat ResRel generation.
- [`data_bool/`](data_bool/): Boolean-operator analyses across CA1 and SUB units.
- [`data_jitter/`](data_jitter/): jitter perturbation experiments.
- [`src/`](src/): figure-generation and exploratory notebooks/scripts.
- [`figures_plotdata/`](figures_plotdata/): derived CSVs used for plots.
- [`data_model/`](data_model/): MATLAB model outputs and plotting helpers.

## Conventions
- Preserve the current script-first style; this repo is not organized as a Python package.
- Most Python scripts import shared code with `sys.path.append('../lib')`.
- Common parameters are intentionally hard-coded, especially `chunksize=20` and the fixed rat ID list.
- Outputs are usually written to the current working directory with names that match existing JSON files.
- If a data schema changes, update the direct consumers in the same edit slice.

## Known Pitfalls
- [`data_bool/BoolOP_ResRel_gen.py`](data_bool/BoolOP_ResRel_gen.py) uses `../../Codes_5/data_extract/` instead of the in-repo data path. Verify paths before assuming a fresh clone is runnable.
- Long-running scripts use multiprocessing and full-rat loops. Prefer syntax checks or small-scope validation first.
- Notebooks are the closest thing to integration tests; there is no unit-test safety net.

## Good Entry Points
- [`data_gen/ResRel_gen.py`](data_gen/ResRel_gen.py)
- [`lib/func_MSR.py`](lib/func_MSR.py)
- [`data_bool/BoolOP_support_functions.py`](data_bool/BoolOP_support_functions.py)
- [`data_jitter/jitter_gen.py`](data_jitter/jitter_gen.py)
- [`src/Fig_test_info.ipynb`](src/Fig_test_info.ipynb)