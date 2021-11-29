Input: List of paths to MSA files

For each MSA:
-[x] run 100 tree inferences with raxml-ng (default settings) = search trees 
-[x] reevaluate each tree with raxml-ng = candidate trees 
-[x] filter tree topologies + iqtree statistical tests = plausible trees
-[ ] store all trees + llhs + runtimes in a database in case we need anything else! 
-[ ] compute label metrics:
  -[ ] average rfdistance per set
  -[ ] number of unique topologies per set
  -[ ] mean llh per set
  -[ ] store all llhs per set
  -[ ] stdev llh per set
  -[ ] fraction of trees that ended up in the plausible tree set
-[ ] compute msa metrics:
  - [x] number of taxa
  - [x] number of sites
  - [x] percentage of invariant sites
  - [x] percentage of gaps
  - [x] number of patterns
  - [x] msa entropy
  - [x] all column entropies 
  - [x] bollback multinomial
  - [x] treelikeness score
  - [x] char frequencies
  - [ ] sum-of-pairs score

-[ ] compute metrics for a single tree inferece:
-[ ] parsimony trees:
-[ ] generate 100 trees (+ store them in a database)
-[ ] compute parsimony scores and store them
-[ ] average rfdistance


Snakemake Workflow per Dataset:
- [x] Infer 100 trees with raxml-ng (separately, seeds 0-99, 50 pars, 50 random)
  - [ ] Extract Information
- [x] Collect all search trees in one file
- [x] Compute rfdistance + number of unique topologies for the 100 search trees
- [x] Evaluate all trees with raxml-ng (seed = 0 for each)
  - [ ] Extract Information
- [x] Collect all eval trees in one file
- [x] Compute rfdistance + number of unique topologies for the 100 eval trees
- [x] Perform significance tests on the 100 eval trees
  - [ ]  Extract Information
- [ ] Compute rfdistance + number of unique topologies for the plausible trees
- [ ] Compute MSA Metrics
- [ ] Infer a single tree and compute the tree metrics