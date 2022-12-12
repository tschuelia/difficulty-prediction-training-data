# Difficulty Prediction Training Data Pipeline
This pipeline generates training data for predicting the difficulty of phylogenetic data.
The goal is to predict the difficulty of a dataset in order to improve the runtime and resource usage of phylogenetic inference.
This pipeline is used to generate the training data for [Pythia](https://github.com/tschuelia/PyPythia). For details on Pythia, see the linked repository or the preprint publication below. 

### Installation
1. Clone this repo: `git clone https://github.com/tschuelia/difficulty-prediction-training-data.git`
2. Install RAxML-NG by following the instructions in [the GitHub repo](https://github.com/amkozlov/raxml-ng).
3. Install IQ-Tree by following the instructions on [their website](http://www.iqtree.org).
4. Setup the conda environment:
    ```
    conda env create -f environment.yml
    ```
   Tipp: use mamba for faster dependency solving. You can simply install mamba by running `conda install mamba -c conda-forge` and then replace `conda` in this command with `mamba`.
5. Activate the new conda environment:
   ```
   conda activate difficulty
   ```
6. Install the required R-Package RPANDAS by running `python install_rpanda.py`. This might take a few minutes to finish. Please check the output of this script for any errors. 
Installing R packages from python can be a bit messy. If it complains that some R package, e.g. `randomPackage` is missing, try installing it using mamba and conda-forge and adding the prefix `r-`: `mamba install r-randomPackage -c conda-forge`  If there is any issue feel free to contact me.

#### System Requirements
The pipeline in this state currently only works on x86 unix machines. It does especially not run on OSX-ARM machines. 
This is due to missing cross-compilations of required R packages in conda-forge.


### Running the pipeline
1. Configure the pipeline by changing the `config.yaml` file:
   * Provide the paths to the MSA files in `msa_paths`.
   * You can either provide partition files for each of the MSAs, or let the pipeline automatically set the model. Note that you have to either provide partitions for all MSAs or for none.
   In this case, all DNA MSAs will be analyzed using the `GTR+G` (RAxML-NG) and `GTR+FO+G4` (IQ-Tree) model, respectively for protein data `LG+G` (RAxML-NG) and `LG+FO+G4` (IQ-Tree) will be used. 
   For morphological data the pipeline uses `MULTI{num_states}_GTR` in RAxML-NG and `MK` in IQ-Tree.
   * You can change the `outdir` variable to tell snakemake where to store the output files. The default is a folder called `results` in the current workdir.
   * In the `software` section provide the paths to executables of RAxML-NG, IQ-Tree from the above installs.
   * The pipeline will also infer some tree metrics for all inferred trees. You can decide whether these metrics should include bootstrap support values by setting `bootstrap_based_metrics`. 
   Due to the runtime overhead of bootstrapping, we recommend setting this to True only for small MSAs (e.g. morphological data). Note that it is turned ON per default!
2. Run `snakemake -n --quiet` for a dry run of snakemake. Snakemake will print a summary of all tasks it will execute to the console.
3. Finally, start the pipeline with `snakemake --cores [num_cores]`.    
If you intend to run snakemake on a slurm cluster, you might want to check out [my instructions](https://github.com/tschuelia/snakemake-on-slurm-clusters) on how to set up snakemake for slurm.


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


## Publication
The paper explaining the details of Pythia is published in MBE:    
Haag, J., Höhler, D., Bettisworth, B., & Stamatakis, A. (2022). **From Easy to Hopeless - Predicting the Difficulty of Phylogenetic Analyses.** *Molecular Biology and Evolution*, 39(12). [https://doi.org/10.1093/molbev/msac254](https://doi.org/10.1093/molbev/msac254)
