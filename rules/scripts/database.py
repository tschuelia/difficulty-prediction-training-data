import peewee as P
from playhouse.sqlite_ext import JSONField

db = P.SqliteDatabase(None)


class Dataset(P.Model):
    uuid = P.UUIDField()
    verbose_name = P.TextField(null=True)
    data_type = P.TextField(null=True)

    # Label features
    num_searches = P.IntegerField(null=True)

    avg_rfdist_search = P.FloatField(null=True)
    num_topos_search = P.IntegerField(null=True)
    mean_llh_search = P.FloatField(null=True)
    std_llh_search = P.FloatField(null=True)

    avg_rfdist_eval = P.FloatField(null=True)
    num_topos_eval = P.IntegerField(null=True)
    mean_llh_eval = P.FloatField(null=True)
    std_llh_eval = P.FloatField(null=True)

    avg_rfdist_plausible = P.FloatField(null=True)
    num_topos_plausible = P.IntegerField(null=True)
    mean_llh_plausible = P.FloatField(null=True)
    std_llh_plausible = P.FloatField(null=True)
    num_trees_plausible = P.IntegerField(null=True)
    proportion_plausible = P.FloatField(null=True)

    # Single inference features
    num_slow_spr_rounds = P.IntegerField(null=True)
    num_fast_spr_rounds = P.IntegerField(null=True)
    llh_starting_tree = P.FloatField(null=True)
    llh_final_tree = P.FloatField(null=True)
    rfdistance_starting_final = P.FloatField(null=True)
    llh_difference_starting_final = P.FloatField(null=True)
    rate_heterogeneity_final = JSONField(null=True)
    eq_frequencies_final = JSONField(null=True)
    substitution_rates_final = JSONField(null=True)
    average_branch_length_final = P.FloatField(null=True)
    std_branch_length_final = P.FloatField(null=True)
    total_branch_length_final = P.FloatField(null=True)
    minimum_branch_length_final = P.FloatField(null=True)
    maximum_branch_length_final = P.FloatField(null=True)
    newick_starting = P.TextField(null=True)
    newick_final = P.TextField(null=True)

    # MSA Features
    num_taxa = P.IntegerField(null=True)
    num_sites = P.IntegerField(null=True)
    num_patterns = P.IntegerField(null=True)
    proportion_gaps = P.FloatField(null=True)
    proportion_invariant = P.FloatField(null=True)
    entropy = P.FloatField(null=True)
    column_entropies = JSONField(null=True)
    bollback = P.FloatField(null=True)
    treelikeness = P.FloatField(null=True)

    # Parsimony Trees Features
    avg_rfdist_parsimony = P.FloatField(null=True)
    num_topos_parsimony = P.IntegerField(null=True)
    mean_parsimony_score = P.FloatField(null=True)
    std_parsimony_score = P.FloatField(null=True)

    class Meta:
        database = db


class RaxmlNGTree(P.Model):
    uuid = P.UUIDField()
    dataset = P.ForeignKeyField(Dataset)
    dataset_uuid = P.UUIDField()
    starting_type = P.CharField(choices=[("random", "random"), ("parsimony", "parsimony")])
    newick_search = P.TextField(null=True)
    llh_search = P.FloatField(null=True)
    compute_time_search = P.FloatField(null=True)
    newick_eval = P.TextField(null=True)
    llh_eval = P.FloatField(null=True)
    compute_time_eval = P.FloatField(null=True)
    is_best = P.BooleanField(null=True)

    # significance tests
    plausible = P.BooleanField(null=True)

    bpRell = P.FloatField(null=True)
    bpRell_significant = P.BooleanField(null=True)
    pKH = P.FloatField(null=True)
    pKH_significant = P.BooleanField(null=True)
    pSH = P.FloatField(null=True)
    pSH_significant = P.BooleanField(null=True)
    pWKH = P.FloatField(null=True)
    pWKH_significant = P.BooleanField(null=True)
    pWSH = P.FloatField(null=True)
    pWSH_significant = P.BooleanField(null=True)
    cELW = P.FloatField(null=True)
    cELW_significant = P.BooleanField(null=True)
    pAU = P.FloatField(null=True)
    pAU_significant = P.BooleanField(null=True)

    class Meta:
        database = db


class ParsimonyTree(P.Model):
    uuid = P.UUIDField()
    dataset = P.ForeignKeyField(Dataset)
    dataset_uuid = P.UUIDField()
    newick_tree = P.TextField(null=True)
    parsimony_score = P.FloatField(null=True)
    compute_time = P.FloatField(null=True)

    class Meta:
        database = db

