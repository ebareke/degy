# degy 

[![GitHub stars](https://img.shields.io/github/stars/ebareke/deg_tools.svg)](https://github.com/ebareke/degy/stargazers) [![GitHub license](https://img.shields.io/github/license/ebareke/degy.svg)](https://github.com/ebareke/degy/blob/master/LICENSE)

This repository contains a set of ready-to-run R scripts, that automates differential gene/transcript expression and related enrichment analyses and generates various (ready-to-publish) summary plots. 

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## Introduction

A companion sets of R scripts for post-processing RNA-seq analyses.

## Features

This set allows to:

- Perform Differential Expression analysis using 3 different approaches (DESeq2, edgeR & limma+voom).
- Do Enrichment Analysis (using gprofiler2 and webgestalt package) and Provide summary results and plots.
- Can be easily integrated in dashboard reports or pipelines (runnable as a CLI command).

## Installation

There is no specific installation required to run tools in this repository.

## Usage

E.g. to run cnt2deg.R script, you will need to open a terminal and navigate to the directory where the script is located. Then, you can run the script using the following command:

``` r
Rscript DEA.r -f root_folder -c condition_1,condition_2 -r reference_condition -fc fold_change_threshold -p adjusted_pval_threshold -m min_avg_coverage
```

### Replace:

- `root_folder` with a folder containing your count matrix and your sample information (`rawCounts.t*` and `sampleInfo.t*` files), 
- `condition_1` with a meaningful description of the condition you consider to be one of your conditions to do differential expression analysis and should be a value in the `Condition` column of the `sampleInfo.t*` file, 
- `condition_2`,with a meaningful description of the condition you consider to be one of your conditions to do differential expression analysis and should be a value in the `Condition` column of the `sampleInfo.t*` file, 
- `reference_condition` with a meaningful description of the condition you consider to be baseline condition to do differential expression analysis and should be a value in the `Condition` column of the `sampleInfo.t*` and be equal either to `condition_1` or  `condition_2`, 
- `fold_change_threshold` with a desired level of fold change cutoff, 
- `adjusted_pval_threshold` with a desired level of adjusted p value cutoff, and 
- `min_avg_coverage` with a desired minimumn coverage one would need to filter low counts.

## Contributing

We welcome contributions to this project! If you're interested in contributing, here are a few ways you can help:

- **Report bugs**: If you find a bug, please [open an issue](https://github.com/ebareke/deg_tools/issues) and let us know. Be sure to include steps to reproduce the bug, and any error messages or other output that might be relevant.

- **Fix bugs**: If you have expertise in debugging and want to help fix existing bugs, please feel free to open a pull request with your fix.

- **Implement new features**: If you have an idea for a new feature, please open an issue to discuss it. If the feature is accepted, you can then open a pull request with your implementation.

- **Improve documentation**: If you see an opportunity to improve the documentation, please feel free to provide suggestions.

We appreciate your help and we will try to respond to all issues and pull requests as quickly as possible. However, please keep in mind that the maintainers of this project have full-time jobs and may not be able to respond immediately. We appreciate your patience and understanding.


## License

MIT License

Copyright (c) 2023 Eric BAREKE (PhD), Emma Carlson, and Prof. Majewski Lab (McGill University)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

