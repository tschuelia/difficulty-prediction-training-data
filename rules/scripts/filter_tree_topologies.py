import pickle

from raxml_parser import get_raxml_num_unique_topos, read_rfdistances

def get_rfdist_clusters(rfdistances, trees):
    # instead of indexing via the eval_tree.id
    # use the newick string
    # this is slower but more robust
    clusters = []
    for (t1, t2), dist in rfdistances.items():
        tree1 = trees[t1].strip()
        tree2 = trees[t2].strip()
        seen_t1 = False
        seen_t2 = False
        for s in clusters:
            if tree1 in s:
                seen_t1 = True
                if dist == 0:
                    s.add(tree2)
                    seen_t2 = True

            if tree2 in s:
                seen_t2 = True
                if dist == 0:
                    s.add(tree1)
                    seen_t1 = True

        if not seen_t1:
            if dist == 0:
                clusters.append({tree1, tree2})
                seen_t1 = True
                seen_t2 = True
            else:
                clusters.append({tree1})
                seen_t1 = True

        if not seen_t2:
            clusters.append({tree2})
            seen_t2 = True

    # remove duplicates
    removed_duplicates = []
    for s in clusters:
        if s not in removed_duplicates:
            removed_duplicates.append(s)

    # check if sets are disjoint
    union = set().union(*removed_duplicates)
    n = sum(len(s) for s in removed_duplicates)
    assert n == len(union)

    return removed_duplicates

eval_trees = [l.strip() for l in open(snakemake.input.all_eval_trees).readlines()]
log_file = snakemake.input.eval_trees_rfdistances_log
rfdists_file = snakemake.input.eval_trees_rfdistances

num_trees = len(eval_trees)

if num_trees > 1:
    num_topos = get_raxml_num_unique_topos(log_file)
    abs_pairwise, _ = read_rfdistances(rfdists_file)

    if num_topos == 1:
        clusters = [set(eval_trees)]

    else:
        clusters = get_rfdist_clusters(abs_pairwise, eval_trees)

else:
    clusters = [set(eval_trees)]
    num_topos = 1

# for each cluster: keep only one tree as representative of the cluster
unique_trees = [next(iter(cluster)) for cluster in clusters]

# sanity checks
assert len(clusters) == num_topos
assert len(unique_trees) == num_topos
assert sum([len(s) for s in clusters]) <= num_trees

open(snakemake.output.filtered_trees, "w").write("\n".join(unique_trees))

with open(snakemake.output.clusters, "wb") as f:
    pickle.dump(clusters, f)