# snpArcher -- standalone-qc branch

Howdy y'all,  the `standalone-qc` branch of Eric's fork of this project just has a few little modifications
that Eric C. Anderson made so that the super-awesome, quick-and-dirty qc module in the
snpArcher workflow can be run on any VCF.gz or BCF file that you happen to have, without
having to run the whole snpArcher pipeline. It is also setup so that it should be easy
to process multiple VCF files (different species, different reference genomes) with only
minimal configuration overhead (preparing a csv file with the paths to the VCF and the `.fai`
files.

This branch has a directory, `standalone-qc-example` that has a vcf.gz
file in it with 16 chinook salmon at 4 pieces of chromosome, and then some scaffolds, which
should be useful for testing and demonstration.

To test this out on your system, the steps are:

1.  Clone this fork of the snpArcher repo:
    ```sh
    # this is using the SSH address
    git clone git@github.com:eriqande/snpArcher.git
    ```
2.  `cd` into the resulting snpArcher repo and switch to the `standalone-qc` branch:
    ```sh
    cd snpArcher
    git checkout standalone-qc
    ```
3.  Dry-run the standalone qc example by running the Snakefile in the `workflow/qc` directory, like this:
    ```sh
    snakemake -np --cores 8 --use-conda -s workflow/modules/qc/Snakefile --configfile standalone-qc-example/config.yaml
    ```
    It is probably worth explaining that line a little.  It is doing a dry-run `-n` and printing
    out the shell commands `-p`, allowing for up to 8 cores of use, and it is using conda for
    managing software.  The Snakefile is specified using the `-s` option to be
    `workflow/modules/qc/Snakefile` and the config file for this run is in
    `standalone-qc-example/config.yaml`.
    

## standalone-qc configuration

Here we have a little explanation of how to configure things. We discuss this in
the context of the examples in the `standalone-qc-example` directory.

### Main config file

The main config file is a simple YAML file.  The example version can be viewed
[here](https://github.com/eriqande/snpArcher/blob/standalone-qc/standalone-qc-example/config.yaml).
You probably should open that in a new tab so you can refer back to it.

That file has the following YAML keys:

- The `vcf_list` is the path to a CSV file that holds information about all the VCF files
and all the reference genomes that you are using.  We will discuss it below.
- The `fai_path` is the path to a directory where any fai files from the reference genomes
  you have used are stored.  They must be named according to the entries in the
  `genome` column of the CSV file, `vcf_list`.
- `scaffolds_to_exclude` is the path to a file that contains scaffold names (one per line)
  that you with to exclude from these qc analyses.


### The `vcf_list` file

This file is a CSV file. The example version can be viewed [here](https://github.com/eriqande/snpArcher/blob/standalone-qc/standalone-qc-example/vcf-list.csv).

This file must be a CSV file and it should contain a header row and then an
additional row for each VCF file that you wish to process.







# snpArcher

<img src="./docs/img/logo.png" alt="snpArcher logo" height="300"/>


snpArcher is a reproducible workflow optimized for nonmodel organisms and comparisons across datasets, built on the [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html#) workflow management system. It provides a streamlined approach to dataset acquisition, variant calling, quality control, and downstream analysis.

### Usage
For usage instructions and complete documentation, please visit our [docs](https://snparcher.readthedocs.io/en/latest/).

### Datasets generated by snpArcher
A number of resequencing datasets have been run with snpArcher generating consistent variant calls, available via [Globus](https://www.globus.org/) in the [Comparative Population Genomics Data collection](https://app.globus.org/file-manager?origin_id=a6580c44-09fd-11ee-be16-195c41bc0be4&origin_path=%2F). Details of data processing are described [in our manuscript](https://www.biorxiv.org/content/10.1101/2023.06.22.546168v1). If you use any of these datasets in your projects, please cite both the [snpArcher paper](https://www.biorxiv.org/content/10.1101/2023.06.22.546168v1) and the original data producers.

### Citing snpArcher
- Cade D Mirchandani, Allison J Shultz, Gregg W C Thomas, Sara J Smith, Mara Baylis, Brian Arnold, Russ Corbett-Detig, Erik Enbody, Timothy B Sackton, A fast, reproducible, high-throughput variant calling workflow for population genomics, Molecular Biology and Evolution, 2023;, msad270, https://doi.org/10.1093/molbev/msad270
- Also, make sure to cite the tools you used within snpArcher.
