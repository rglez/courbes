# courbes 

[![DOI](https://zenodo.org/badge/781865659.svg)](https://doi.org/10.5281/zenodo.14003027)



> Automated Statistics Extraction from (Multi-Replica) MD Simulations with Curves+

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
