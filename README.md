# deg_tools

[![GitHub stars](https://img.shields.io/github/stars/ebareke/deg_tools.svg)](https://github.com/ebareke/deg_tools/stargazers) [![GitHub license](https://img.shields.io/github/license/ebareke/deg_tools.svg)](https://github.com/ebareke/deg_tools/blob/master/LICENSE)

This repository contains a set of ready-to-run R scripts, that automate differential gene expression and enrichment analyses and generate various QC and summary plots. The scripts are be combined or run separated and are easy to integrate into pipelines.

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Documentation](#documentation)
- [Support](#support)
- [Contributing](#contributing)
- [License](#license)

## Introduction

A detailed introduction to your project. This should give the reader an idea of what your project is all about and why it is useful.

## Features

- Perform Differential Expression analysis using 3 different approaches (DESeq2, edgeR & limma+voom).
- Do Enrichment Analysis (using gprofiler2 and webgestalt package) and Provide summary results and plots.
- Can be easily integrated in dashboard reports or pipelines (runnable as a CLI command).

## Installation

There is no specific installation required to run tools in this rep.
Just download them and install below packages:

``` r
# install the packages (if you haven't already)
pkgs1 <- c("gplots", "gprofiler2", "ggplot2")
install.packages(pkgs1)

install.packages("BiocManager")
pkgs2 <- c("DESeq2", "edgeR", "limma", "voom", "pcaMethods")
BiocManager::install(pkgs2)
```

## Usage

E.g. to run cnt2deg.R script, you will need to open a terminal and navigate to the directory where the script is located. Then, you can run the script using the following command:

``` r
Rscript cnt2deg.R count_matrix_file sample_metadata_file lfc_cutoff adjp_cutoff
```

Replace `count_matrix_file`, `sample_metadata_file`, `lfc_cutoff`, and `adjp_cutoff` with the paths to your count matrix file, sample metadata file, and the desired cutoffs for log fold change and adjusted p-value, respectively.


Alternatively, you can also run the above script from within the R console by using the `source()` function. For example:

``` r
args <- c("count_matrix_file", "sample_metadata_file", "lfc_cutoff", "adjp_cutoff")
source("script.R", local = list(args = args))
```


## Documentation

Coming soon (when project complete).

## Support

If you have any questions or need help with the project, here are a few options:

- Search for solutions on [Stack Overflow](https://stackoverflow.com/) or other online communities.
- If you've found a bug or have a feature request, please [open an issue](https://github.com/ebareke/deg_tools/issues).
- If you want to contribute to the project, please see the [Contributing](#contributing) section below.

We try to respond to all issues and pull requests as quickly as possible. However, please keep in mind that the maintainers of this project have full-time jobs and may not be able to respond immediately. We appreciate your patience and understanding.


## Contributing

We welcome contributions to this project! If you're interested in contributing, here are a few ways you can help:

- **Report bugs**: If you find a bug, please [open an issue](https://github.com/ebareke/deg_tools/issues) and let us know. Be sure to include steps to reproduce the bug, and any error messages or other output that might be relevant.

- **Fix bugs**: If you have expertise in debugging and want to help fix existing bugs, please feel free to [open a pull request](https://github.com/ebareke/deg_tools/pulls) with your fix.

- **Implement new features**: If you have an idea for a new feature, please [open an issue](https://github.com/ebareke/deg_tools/issues) to discuss it. If the feature is accepted, you can then [open a pull request](https://github.com/ebareke/deg_tools/pulls) with your implementation.

- **Improve documentation**: If you see an opportunity to improve the documentation, please feel free to [open a pull request](https://github.com/ebareke/deg_tools/pulls) with your changes.

We appreciate your help! Please be sure to read our [contribution guidelines](https://github.com/ebareke/deg_tools/blob/master/CONTRIBUTING.md) before submitting your pull request.

## License

Copyright (c) 2022 Eric Bareke
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

