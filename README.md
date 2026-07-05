# Master thesis code

This repository contains the numerical simulation codes, input files, and representative output files used in my master's thesis.

## Overview

The codes are used to simulate erosion and deposition processes of debris flows.
The repository is organized by the experimental cases used for model validation and comparison.

## Repository structure

```text
master-thesis-code/
├── README.md
├── 田澤実験/
│   ├── txtファイル/
│   ├── expファイル/
│   ├── output修正/
│   ├── output江頭/
│   ├── Egashira.C/
│   └── Hotta.C/
├── 鈴木実験/
│   ├── txtファイル/
│   ├── 水路形状/
│   ├── output江頭堆積/
│   ├── output修正堆積/
│   ├── output鈴木堆積/
│   ├── output江頭侵食/
│   ├── output修正侵食/
│   ├── output鈴木侵食/
│   ├── D_Egashira.C/
│   ├── D_Hotta.C/
│   ├── D_Suzuki.C/
│   ├── E_Egashira.C/
│   ├── E_Hotta.C/
│   └── E_Suzuki.C/
└── 經隆実験/
│   ├── txtファイル/
│   ├── 水路形状/
│   ├── outputFIX江頭/
│   ├── outputFIX修正/
│   ├── outputFIX鈴木/
│   ├── FIX_Egashira.C/
│   ├── FIX_Hotta.C/
│   └── FIX_Suzuki.C/
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
