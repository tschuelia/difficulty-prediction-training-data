# Difficulty Prediction Training Data Pipeline
This pipeline generates training data for predicting the difficulty of phylogenetic data.
The goal is to predict the difficulty of a dataset in order to improve the runtime and resource usage of phylogenetic inference.
This pipeline is used to generate the training data for [Pythia](https://github.com/tschuelia/PyPythia). For details on Pythia, see the linked repository or the preprint publication below. 

### Requirements
1. Setup the conda environment:
    ```
    conda env create -f environment.yml
    ```
2. Install RAxML-NG.
   ```
   git clone --recursive https://github.com/amkozlov/raxml-ng.git
   cd raxml-ng
   mkdir build && cd build
   cmake ..
   make -j
   ```
3. Install IQ-Tree by following the instructions on [their website](http://www.iqtree.org).
4. Install the required R-Package RPANDAS by running `python install_rpanda.py`

### Running the pipeline
1. Configure the pipeline by changing the `config.yaml` file:
   * Provide the paths to the MSA files in `msa_paths`.
   * You can either provide partition files for each of the MSAs, or let the pipeline automatically set the model. Note that you have to either provide partitions for all MSAs or for none.
   In this case, all DNA MSAs will be analyzed using the `GTR+G` (RAxML-NG) and `GTR+FO+G4` (IQ-Tree) model, respectively for protein data `LG+G` (RAxML-NG) and `LG+FO+G4` (IQ-Tree) will be used. 
   For morphological data the pipeline uses `MULTI{num_states}_GTR` in RAxML-NG and `MK` in IQ-Tree.
   * You can change the `outdir` variable to tell snakemake where to store the output files. The default is a folder called `data` in the current workdir.
   * In the `software` section provide the paths to executables of RAxML-NG, IQ-Tree from the above installs.

2. Run `snakemake -n --quiet` for a dry run of snakemake. Snakemake will print a summary of all tasks it will execute to the console.
3. Finally, start the pipeline with `snakemake --cores [num_cores]`.    
If you intend to run snakemake on a slurm cluster, you might want to check out [my instructions](https://github.com/tschuelia/snakemake-on-slurm-clusters) on how to setup snakemake for slurm.


### Output
In the output directory you will find a subdirectory for each MSA you provided as input. Each subdirectory contains a `data.sqlite3` SQLite database and two parquet files:
- `training_data.parquet` contains MSA attributes and the ground-truth difficulty label
- `raxmlng_tree_data.parquet` contains tree characteristics and the results of the statistical tests for all inferred trees. There is also a column labeled `is_best` that indicates whether this tree is the best found ML tree.

The subdirectory `output_files` contains the intermediate files and logs of all snakemake steps. 

#### Collecting the MSA and tree characteristics of the best tree
Instead of including this step in the snakemake pipeline, I provide an extra script for collecting the MSA attributes and tree characteristics for all results in a given directory.
The reason for this is that snakemake will not execute this final step in case of an error during the pipeline execution of a single MSA and some datasets fail with RAxML-NG or IQ-Tree for no apparent reason...

You can execute this final data collection by calling: `python final_data_collection.py --dir [directory]`. Pass the directory you specified in the `outdir` variable in the `config.yaml` file.


### Training Data
This repository also contains the current set training data as parquet file. To open the file, you need `pyarrow` or `fastparquet` installed in your environment. 
Then you can open the file with pandas using `pd.read_parquet("training_data.parquet")`.

Careful, the training data does not only contain all features we used during our prediction experiments, but also the features we used to quantify the difficulty. 
The file `features.py` contains further explanation of the features and lists of features and respective labels that can be used for experimenting with the training data.


## Preprint Publication
The paper explaining the details of Pythia is available as preprint on BioRxiv:   
Haag, J., Höhler, D., Bettisworth, B., & Stamatakis, A. (2022). **From Easy to Hopeless - Predicting the Difficulty of Phylogenetic Analyses.** BioRxiv. [https://doi.org/10.1101/2022.06.20.496790](https://doi.org/10.1101/2022.06.20.496790)
