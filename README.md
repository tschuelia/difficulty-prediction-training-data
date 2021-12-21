# Difficulty Prediction Training Data Pipeline
This pipeline generates training data for predicting the difficulty of molecular datasets.
The goal is to predict the difficulty of a dataset in order to improve the runtime and ressource usage of phylogenetic inference.

### Requirements
1. Setup the conda environment:
    ```
    conda env create -f environment.yml
    ```
2. Install RAxML-NG.
   Note that you need my version of RAxML-NG in order to generate parsimony trees for protein datasets.
   ```
   git clone --recursive https://github.com/tschuelia/raxml-ng.git
   cd raxml-ng
   mkdir build && cd build
   cmake ..
   make -j
   ```
3. Install Parsimonator:
   ```
   git clone https://github.com/stamatak/Parsimonator-1.0.2.git
   cd Parsimonator-1.0.2
   make -f Makefile.gcc
   ```
4. Install IQ-Tree by following the instructions on [their website](http://www.iqtree.org).

### Running the pipeline
1. Configure the pipeline by changing the `config.yaml` file:
   * Provide the paths to the MSA files in `msa_paths`. Note that all MSAs need to be of the same type.
   * Set `data_type` to the respective type. This can be either `DNA` or `AA` for protein data.
   * You can change the `outdir` variable to tell snakemake where to store the output files. The default is a folder called `data` in the current workdir.
   * In the `software` section provide the paths to executables of RAxML-NG, IQ-Tree and Parsimonator from the above installs.
   * For RAxML-NG and IQ-Tree: specify the model of evolution in the respective `model` variable.   
   For DNA data, we recommend `GTR+G` (RAxML-NG) and `GTR+FO+G4` (IQ-Tree)   
   For protein data `LG+G` (RAxML-NG) and `LG+FO+G4` (IQ-Tree)

2. Run `snakemake -n --quiet` for a dry run of snakemake. Snakemake will print a summary of all tasks to the console.
3. Finally, start the pipeline with `snakemake --cores [num_cores]`.    
If you intend to run snakemake on a slurm cluster, you might want to check out [my instructions](https://github.com/tschuelia/snakemake-on-slurm-clusters) on how to setup snakemake for slurm.


### Output
In the output directory you will find a subdirectory for each MSA you provided as input. Each subdirectory contains a `data.sqlite3` SQLite database.
The following image depicts its database schema: [TODO DB Schema]

The subdirectory `output_files` contains the intermediate files and logs of all snakemake steps. 