# courbes 

[![DOI](https://zenodo.org/badge/781865659.svg)](https://doi.org/10.5281/zenodo.14003027)



> Automated Statistics Extraction from (Multi-Replica) MD Simulations with Curves+

> [!NOTE]
> Please note that this project is still under development and is not yet ready for production use. We are working hard
> to make it available as soon as possible. If you are interested in contributing to this project, please refer to
> the [Contributing Guidelines](CONTRIBUTING.md).

## Motivation

**The analysis of DNA dynamics through MD simulations** is crucial for understanding its behavior and function. However,
extracting meaningful insights from these simulations can be a complex and time-consuming process. The existing DNA
analyzer, Curves+, provides a powerful tool for characterizing DNA structure, but its integration into a streamlined
workflow may require significant effort. To address this challenge, **we introduce Courbes, a Python package designed to
simplify and accelerate DNA dynamics analysis using Curves+**. Courbes acts as a user-friendly wrapper for Curves+,
automating the extraction of key descriptors from multiple molecular dynamics replicas, generating basic statistical
reports, and facilitating efficient data analysis. This package empowers researchers to gain deeper insights into DNA
dynamics with minimal effort, enabling them to focus on scientific discoveries rather than technical hurdles.

## Installation

courbes will be soon released as a Python package in PyPI, that will allow users to install it
with `pip install courbes`. In the meantime, you can install it by cloning this repository and running the following
commands in your terminal:

```bash
git clone https://github.com/rglez/courbes.git
cd courbes
conda env create -f environment.yml
conda activate courbes
pip install -e .
```

## Usage

The main functionality of courbes is to extract the DNA descriptors reported by curves+ and its statistics from DNA
simulations. An easy-to-prepare input file is required to run the analysis (see the Documentation section). Once the
input file is ready, you can run the analysis with the following command:

```bash
courbes path-to-config-file.cfg
```

The previous command, will generate a report with the value of each descriptor for each frame, as well as the statistics
of the DNA descriptors extracted from the simulations. The report directory will contain five sub-folders:

- **axis:** Base pair-axis parameters
- **intra:** Intra-base pair parameters
- **inter:** Inter-base pair parameters & Helical rise and twist
- **groove:** Groove geometry
- **backbone:** Backbone parameters

## Documentation

The most detailed and updated documentation can be found [in the Wiki](https://github.com/rglez/courbes/wiki).

## Contributing

We welcome contributions to this project. Please read the [Contributing Guidelines](CONTRIBUTING.md) for more
information.

## Project Organization

```
├── AUTHORS.md              <- List of developers and maintainers.
├── CHANGELOG.md            <- Changelog to keep track of new features and fixes.
├── CONTRIBUTING.md         <- Guidelines for contributing to this project.
├── Dockerfile              <- Build a docker container with `docker build .`.
├── LICENSE.txt             <- License as chosen on the command-line.
├── README.md               <- The top-level README for developers.
├── configs                 <- Directory for configurations of model & application.
├── data
│   ├── external            <- Data from third party sources.
│   ├── interim             <- Intermediate data that has been transformed.
│   ├── processed           <- The final, canonical data sets for modeling.
│   └── raw                 <- The original, immutable data dump.
├── docs                    <- Directory for Sphinx documentation in rst or md.
├── environment.yml         <- The conda environment file for reproducibility.
├── models                  <- Trained and serialized models, model predictions,
│                              or model summaries.
├── notebooks               <- Jupyter notebooks. Naming convention is a number (for
│                              ordering), the creator's initials and a description,
│                              e.g. `1.0-fw-initial-data-exploration`.
├── pyproject.toml          <- Build configuration. Don't change! Use `pip install -e .`
│                              to install for development or to build `tox -e build`.
├── references              <- Data dictionaries, manuals, and all other materials.
├── reports                 <- Generated analysis as HTML, PDF, LaTeX, etc.
│   └── figures             <- Generated plots and figures for reports.
├── scripts                 <- Analysis and production scripts which import the
│                              actual PYTHON_PKG, e.g. train_model.
├── setup.cfg               <- Declarative configuration of your project.
├── setup.py                <- [DEPRECATED] Use `python setup.py develop` to install for
│                              development or `python setup.py bdist_wheel` to build.
├── src
│   └── courbes             <- Actual Python package where the main functionality goes.
├── tests                   <- Unit tests which can be run with `pytest`.
├── .coveragerc             <- Configuration for coverage reports of unit tests.
├── .isort.cfg              <- Configuration for git hook that sorts imports.
└── .pre-commit-config.yaml <- Configuration of pre-commit git hooks.
```
