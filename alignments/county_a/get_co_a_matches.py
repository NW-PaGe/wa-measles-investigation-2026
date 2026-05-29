import baltic as bt
from matplotlib import pyplot as plt
import matplotlib.lines as mlines
import pandas as pd

def main():
    tree_data = bt.loadJSON(
        "../../phylogenetic/auspice/measles_genome.json",
        json_translation={
            "name": "name",
            "absoluteTime": "num_date",
            "divergence": "div",
        },
        verbose=False,
        sort=True,
        stats=True,
    )
    tree = tree_data[0]

    metadata = pd.read_csv("../../tree_fig/outbreak_metadata.tsv", sep="\t")

    # the auspice .json does not create keys where there are no values. So, if there is no county
    # information, then there will be key errors when filtering data by county. Create those keys
    # in the dictionary and assign a default NA value.
    for l in tree.Objects:
        if l.is_leaf():
            l.traits.setdefault("node_attrs", {}).setdefault("division", {}).setdefault(
                "value", "NA"
            )

    for l in tree.Objects:
        if l.is_leaf():
            l.traits.setdefault("node_attrs", {}).setdefault("county", {}).setdefault(
                "value", "NA"
            )

    # Add county metadata
    def add_county_metadata():
        for l in tree.Objects:
            if l.is_leaf() and l.name in metadata["accession"].values:
                row = metadata[metadata["accession"] == l.name]
                if not row.empty and "county" in metadata.columns:
                    l.traits["node_attrs"]["county"]["value"] = row["county"].values[0]

    add_county_metadata()

    ###  The tree file is good now. Make a sub-tree for just the  ###
    ###  Washington sequences analyseed from 2026 (for paper)     ###

    # Find the most recent common ancestor to the Washington sequences from 2026
    def filter_outbreak_leaves(tree):
        outbreak_leaves = []
        for l in tree.Objects:
            if l.is_leaf():
                if l.traits["node_attrs"]["county"]["value"] in [
                    "County C",
                    "County D",
                    "County A",
                ]:
                    outbreak_leaves.append(l)
        return outbreak_leaves

    outbreak_leaves = filter_outbreak_leaves(tree=tree)

    # Find the MRCA for outbreak sequences
    mrca = tree.commonAncestor(outbreak_leaves)

    print(f"the mrca node for the outbreak sequences is {mrca} \n")

    subtree = tree.subtree(k=mrca)
    subtree.sortBranches()

    def get_divergence_stats(tree):
        print(f"checking for missing divergence values in {tree} \n")
        divergences = []
        for thing in tree.Objects:
            divergences.append(thing.divergence)
        result = dict(
            {
                "minimum": min(divergences),
                "maximum": max(divergences),
                "mean": sum(divergences) / len(divergences),
            }
        )
        return result
    divergence_stats = get_divergence_stats(tree=subtree)

    print(f"divergence stats for Wa outbreak tree: {divergence_stats} \n")

    def divergence_recalc(tree: bt.tree):
        root_div = tree.root.divergence
        for ob in tree.Objects:
            if ob == tree.root:
                ob.divergence=0
            else:
                ob.divergence = ob.divergence - root_div
        try:
            tree.root.parent.divergence = tree.root.parent.divergence - root_div
        except NameError:
            print("This tree does not have a stem, everything else adjusted.")
        print(f"tree root divergence is {tree.root.divergence}")
        return tree
    
    print(f"Updating the divergence values with respect to new root \n")
    subtree = divergence_recalc(subtree)
    recalc_div_stats = get_divergence_stats(subtree)
    print(f"The new divergence stats are {recalc_div_stats}")

    county_a_matches = [l for l in subtree.Objects if l.is_leaf() and l.divergence==0]
    n_matches=len(county_a_matches)

    accessions = [l.name for l in county_a_matches]

    with open("county_a_matches.txt", "w") as f:
        for name in accessions:
            f.write(f"{name}\n")

    print(f"There were {n_matches} matches, including the Wa sequence")

if __name__ == "__main__":
    main()