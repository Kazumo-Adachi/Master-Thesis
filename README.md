# Master thesis code

This repository contains the numerical simulation codes, input files, and representative output files used in my master's thesis.

## Overview

The codes are used to simulate erosion and deposition processes of debris flows.
The repository is organized by the experimental cases used for model validation and comparison.

## Repository structure

```text
master-thesis-code/
├── README.md
├── Tazawa_experiment/
│   ├── README.md
│   ├── source_code/
│   ├── input_txt/
│   ├── experimental_data_txt/
│   └── output_example_txt/
├── Suzuki_experiment/
│   ├── README.md
│   ├── source_code/
│   ├── input_txt/
│   ├── experimental_data_txt/
│   └── output_example_txt/
└── docs/
```

## Contents

### `Tazawa_experiment/`

This folder contains the source codes, input files, experimental data, and representative output files used for the reproduction analysis of the Tazawa experiment.

### `Suzuki_experiment/`

This folder contains the source codes, input files, experimental data, and representative output files used for the reproduction analysis of the Suzuki experiment.

### `docs/`

This folder contains supplementary documentation, including notes on the numerical model, calculation conditions, and reproduction procedures.

## Unit system

The source codes use the cm-g-sec unit system unless otherwise noted.

## How to run

Compile each C source file using a C compiler.

Example:

```bash
gcc source_code/FIX_Egashira.c -o FIX_Egashira
./FIX_Egashira
```

The exact file name and execution procedure may differ depending on the experimental case and model.

## Notes

This repository corresponds to the version of the code used for the master's thesis.

Large raw data files, videos, and complete output files are not included in this repository.
Only the source codes, input files, and representative output examples necessary to understand and reproduce the main calculations are provided.
