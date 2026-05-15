import baltic as bt
from matplotlib import pyplot as plt
import matplotlib.lines as mlines
import pandas as pd

def main():
    tree_data = bt.loadJSON(
        "../phylogenetic/auspice/measles_genome.json",
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

    metadata = pd.read_csv("./outbreak_metadata.tsv", sep="\t")

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

    ###  Outbreak sub-tree is extracted. Now     ###
    ###  Collapse the large polytomies from the  ###
    ###  South Carolina sequences into single    ###
    ###  leaves for plotting as a summary        ###

    def merge_polytomy_layers(tree, attribute="divergence", verbose=False):
        """
        Remove intermediate nodes that share the same attribute value as their parent.
        This grafts children upward, combining polytomy layers under the same value.
        """
        if verbose:
            print(f"Merging polytomy layers based on {attribute}...")
        nodes_to_remove = []
        
        for node in tree.getInternal():
            if node.parent is not None:
                node_val = node.divergence
                parent_val = node.parent.divergence
                
                if node_val == parent_val and node_val is not None:
                    # This node has the same value as its parent - mark for removal
                    nodes_to_remove.append(node)
        
        # Remove marked nodes by grafting their children to their parent
        for node in nodes_to_remove:
            if verbose:
                print(f"  Removing {node.name} (keeping its {len(node.children)} children)")
            # Reassign node's children to its parent
            for child in list(node.children):
                child.parent = node.parent
                node.parent.children.append(child)
            # Remove node from its parent's children
            node.parent.children.remove(node)
            # Remove from tree objects
            tree.Objects.remove(node)
        
        if verbose:
            print(f"  Removed {len(nodes_to_remove)} intermediate nodes")
        return tree

    def rake_leaves(input_tree: bt.tree, 
                    threshold=2, 
                    verbose=False):
        input_tree.fixHangingNodes()
        if verbose is True:
            print(f"input tree {input_tree} stats are: \n")
            input_tree.treeStats()
            print(f"subsetting tree one node at a time")
        nodes = [k for k in input_tree.getInternal() if k.is_node()]
        if verbose is True:
            print(f"All nodes in tree: {nodes} \n")
        for node in nodes:
            if verbose is True:
                print(f"node: {node.name}, height: {node.height}")
            leaves_crit = [
                k
                for k in node.children
                if k.is_leaf()
                and k.traits["node_attrs"]["division"]["value"] == "South Carolina"
                and k.divergence == node.divergence
            ]
            # get the divergence values for the leaves meeting the filter criteria:
            leaves_div = []
            for l in leaves_crit:
                leaves_div.append(l.divergence)  # get all divergence values for leaves
            # get the counts for each unique value
            counts = {}  
            for l in leaves_crit:
                if l.divergence in counts:
                    counts[l.divergence] += 1
                else:
                    counts[l.divergence] = 1
            # make a dict to log what was removed
            removed = {}
            # loop through leaves and remove if it's divergence count > 1
            for l in leaves_crit:
                if counts[l.divergence] <= threshold-1:
                    pass
                else:
                    # start here - if not yet in log, keep and assign count from counts dict
                    # else remove
                    if l.divergence not in removed.keys():
                        removed.update({l.divergence: 0})
                        if verbose is True:
                            print(f"leaf attributes are {l.traits}")
                        l.traits = {}
                        if verbose is True:
                            print(f"confirm the leaf traits were cleared: {l.traits}")
                        l.traits.setdefault("node_attrs", {}).setdefault(
                            "division", {}
                        )["value"] = "South Carolina"
                        l.traits.update({"count": counts[l.divergence], "combo": True})
                        if verbose is True:
                            print(f"The final leaf traits are: {l.traits}")
                    else:
                        removed[l.divergence] += 1
                        input_tree.Objects.remove(l)
                        node.children.remove(l)
        # assert here that the number removed + 1 equals the number in counts table for the div value
        input_tree.sortBranches()
        for l in input_tree.Objects:
            if l.is_leaf():
                l.traits.setdefault("count", None)
                l.traits.setdefault("combo", False)
        # print(f"the new tree's stats are:")
        input_tree.treeStats()
        # print(County C_subtree.Objects)
        return input_tree

    # Apply polytomy layer merging before rake_leaves
    subtree = merge_polytomy_layers(subtree, attribute="divergence", verbose=False)
    subtree_summarized = rake_leaves(subtree, verbose=False, threshold=3)

    def color_node(node):
        """
        Color nodes based on state and/or Washington County
        """
        try:
            division = node.traits["node_attrs"]["division"]["value"]
            county = node.traits["node_attrs"]["county"]["value"]
            if division == "Washington":
                if county == "County C":
                    return "#fb8072"
                elif county == "County D":
                    return "#8dd3c7"
                elif county == "County A":
                    return "#ffffb3"
                else:
                    return "white"
            elif division == "South Carolina":
                return "#bebada"
        except (KeyError, TypeError):
            pass
        return "grey"

    fig, ax = plt.subplots(figsize=(7.5, 10.5), facecolor="w", layout="constrained")
    x_attr = lambda k: k.divergence
    
    subtree_summarized.plotTree(ax, x_attr=x_attr, colour="gray")
    #plotting just the combo leaves
    subtree_summarized.plotPoints(
         ax,
         x_attr=x_attr,
         size=20,
         colour="#bebada",
         marker='o',
         zorder=77,
         target=lambda k: k.is_leaf() and k.traits["combo"] is True,
         outline_colour='black'
    )
    #plot just the non-combination leaves
    subtree_summarized.plotPoints(
        ax,
        x_attr=x_attr,
        size=15,
        colour=color_node, # callable defined above
        marker='o',
        zorder=78,
        target=lambda k: k.is_leaf() and k.traits["combo"] is False,
    )
    kwargs={'va': 'center', 'ha': 'left', 'size': 10}
    subtree_summarized.addText(
        ax,
        x_attr=x_attr,
        offset=[.01,0],
        color="black",
        target= lambda k: k.is_leaf() and k.traits["combo"] is True,
        text= lambda k: '-'+str(k.traits['count']),
        zorder=81,
        **kwargs
    )
    ax.xaxis.tick_bottom()
    ax.yaxis.tick_left()

    [ax.spines[loc].set_visible(False) for loc in ["top", "right", "left"]]

    ## These next mlines artists don't actually plot;
    ## it's just to trick the legend into working, since
    ## coloring is multilevel/complex based on a callable
    county_c = mlines.Line2D(
        [],
        [],
        color="#fb8072",
        mec="black",
        marker="o",
        markersize=10,
        linestyle="none",
        label="County C",
    )
    county_a = mlines.Line2D(
        [],
        [],
        color="#ffffb3",
        mec="black",
        marker="o",
        markersize=10,
        linestyle="none",
        label="County A",
    )
    county_d = mlines.Line2D(
        [],
        [],
        color="#8dd3c7",
        mec="black",
        marker="o",
        markersize=10,
        linestyle="none",
        label="County D",
    )
    sc = mlines.Line2D(
        [],
        [],
        color="#bebada",
        mec="black",
        marker="o",
        markersize=10,
        linestyle="none",
        label="South Carolina",
    )

    # Pass the list of custom handles and labels to the legend function
    ax.legend(handles=[county_a, county_c, county_d, sc], frameon=False, loc="center right")
    ax.tick_params(axis="y", size=0)
    ax.set_yticklabels([])
    #ax.set_ylim(-5, subtree.ySpan + 5)
    #plt.tight_layout()
    plt.xlabel("Mutations Compared to Root")
    plt.savefig("full.svg")
    # plt.show()

    ##########################################################
    ###   Full tree is looking good. Start to prune for    ###
    ###   the final figure                                 ###
    ##########################################################

    # remove SC clades (nodes) that don't have Wa leaves or child nodes
    def prune_branches(tree: bt.tree, strip = False):
        def collect_keepers(tree: bt.tree):
            keep=[]
            for k in tree.getInternal():
                leaves = [l for l in k.children if l.is_leaf()]
                leaf_states = [l.traits['node_attrs']['division']['value'] for l in leaves]
                child_nodes=[n for n in k.children if n.is_node()]
                if 'Washington' in leaf_states or len(child_nodes)>0:
                    keep.extend(leaves)
                else:
                    continue
            return keep
        keep = collect_keepers(tree)
        subtree_reduced = tree.reduceTree(keep=keep)
        subtree_reduced.sortBranches()
        if strip is True:
            while len(subtree_reduced.getExternal()) > len(collect_keepers(subtree_reduced)):
                keep = collect_keepers(subtree_reduced)
                subtree_reduced = tree.reduceTree(keep=keep)
                subtree_reduced.sortBranches(sortByHeight=True)
        return subtree_reduced
    
    subtree_reduced = prune_branches(subtree_summarized, strip = True)
    # modify the stock sorting function to sort based on divergence rather than length, since this is a divergence tree
    mod=-1
    sort_function=lambda k: (k.is_node(),-len(k.leaves)*mod,k.divergence*mod) if k.is_node() else (k.is_node(),k.divergence*mod)
    # use the new sort function to order the branches and make them pretty
    subtree_reduced.sortBranches(sort_function=sort_function)
    print(f"The reduced tree stats are {subtree_reduced.treeStats()}")

    fig, ax = plt.subplots(figsize=(7.5, 10.5), facecolor="w", layout="constrained")
    x_attr = lambda k: k.divergence
    subtree_reduced.plotTree(ax, x_attr=x_attr, colour="slategrey")
    #plotting just the combo leaves
    subtree_reduced.plotPoints(
         ax,
         x_attr=x_attr,
         size=55,
         colour="#bebada",
         marker='o',
         zorder=77,
         target=lambda k: k.is_leaf() and k.traits["combo"] is True,
         outline_colour='black'
    )
    #plot just the non-combination leaves
    subtree_reduced.plotPoints(
        ax,
        x_attr=x_attr,
        size=55,
        colour=color_node, # callable defined above
        marker='o',
        zorder=78,
        target=lambda k: k.is_leaf() and k.traits["combo"] is False,
    )
    kwargs={'va': 'center', 'ha': 'left', 'size': 12}
    subtree_reduced.addText(
        ax,
        x_attr=x_attr,
        offset=[.05,0],
        color="black",
        target= lambda k: k.is_leaf() and k.traits["combo"] is True,
        text= lambda k: '-'+str(k.traits['count']),
        zorder=81,
        **kwargs
    )
    ax.xaxis.tick_bottom()
    ax.yaxis.tick_left()

    [ax.spines[loc].set_visible(False) for loc in ["top", "right", "left"]]

    ax.legend(handles=[county_a, county_c, county_d, sc], frameon=False, loc="center right")
    ax.tick_params(axis="y", size=0)
    ax.set_yticklabels([])
    ax.set_xticks([0, 1, 2, 3, 4])
    #ax.set_ylim(-5, subtree.ySpan + 5)
    #plt.tight_layout()
    plt.xlabel("Mutations Compared to Root of Clade")
    plt.savefig("reduced.svg")
    plt.savefig("reduced.png")
    #plt.show()

    ##########################################################
    ###    For nodes with Washington children, remove      ###
    ###    SC leaves that are longer than the max Wa       ###
    ###    leaf height                                     ###
    ##########################################################

    def nodes_to_keep(tree: bt.tree):
        keep=[]
        for k in tree.getInternal():
            leaves = [l for l in k.children if l.is_leaf()]
            leaf_states = [l.traits['node_attrs']['division']['value'] for l in leaves]
            #child_nodes=[n for n in k.children if n.is_node()]
            if 'Washington' in leaf_states:
                wa_leaf_divs=[]
                for l in leaves:
                    if l.traits['node_attrs']['division']['value'] == 'Washington':
                        keep.append(l)
                        wa_leaf_divs.append(l.divergence)
                for l in leaves:
                    if l.traits['node_attrs']['division']['value'] != 'Washington' and l.divergence <= max(wa_leaf_divs):
                        keep.append(l)
        return keep
    
    keep = nodes_to_keep(subtree_reduced)
    print(f"keeping the following leaves: {keep}")
            
    subtree_reduced_again = subtree_reduced.reduceTree(keep=keep)
    subtree_reduced_again.sortBranches()
    print(f"The second reduced tree stats are {subtree_reduced_again.treeStats()}")

    fig, ax = plt.subplots(figsize=(4, 6.5), facecolor="w", layout="constrained")
    x_attr = lambda k: k.divergence
    subtree_reduced_again.plotTree(ax, x_attr=x_attr, colour="slategrey")
    #plotting just the combo leaves
    subtree_reduced_again.plotPoints(
         ax,
         x_attr=x_attr,
         size=20,
         colour='lightgrey',
         marker='o',
         zorder=77,
         target=lambda k: k.is_leaf() and k.traits["combo"] is True,
         outline_colour='black'
    )
    #plot just the non-combination leaves
    subtree_reduced_again.plotPoints(
        ax,
        x_attr=x_attr,
        size=17,
        colour=color_node, # callable defined above
        marker='o',
        zorder=78,
        target=lambda k: k.is_leaf() and k.traits["combo"] is False,
    )
    kwargs={'va': 'center', 'ha': 'left', 'size': 10}
    subtree_reduced_again.addText(
        ax,
        x_attr=x_attr,
        offset=[.01,0],
        color="black",
        target= lambda k: k.is_leaf() and k.traits["combo"] is True,
        text= lambda k: '-'+str(k.traits['count']),
        zorder=81,
        **kwargs
    )
    ax.xaxis.tick_bottom()
    ax.yaxis.tick_left()

    [ax.spines[loc].set_visible(False) for loc in ["top", "right", "left"]]

    ax.legend(handles=[county_a, county_c, county_d, sc], frameon=False, loc="upper right")
    ax.tick_params(axis="y", size=0)
    ax.set_yticklabels([])
    #ax.set_ylim(-5, subtree.ySpan + 5)
    #plt.tight_layout()
    plt.xlabel("Mutations Compared to Root")
    plt.savefig("reduced_again.svg")
    #plt.show()

    ##########################################################
    ###      Extract the County C tree                     ###
    ########################################################## 

    # Find the most recent common ancestor to the County C seqs
    def filter_outbreak_leaves(tree):
        county_c_leaves = []
        for l in tree.Objects:
            if l.is_leaf():
                if l.traits["node_attrs"]["county"]["value"] in ["County C"]:
                    county_c_leaves.append(l)
        return county_c_leaves

    county_c_leaves = filter_outbreak_leaves(tree=tree)
    # Find the MRCA for County C sequences
    mrca_county_c = tree.commonAncestor(county_c_leaves)
    print(f"the mrca node for the County C sequences is {mrca_county_c}")
    county_c_subtree = tree.subtree(k=mrca_county_c)
    county_c_subtree.fixHangingNodes()
    county_c_subtree.treeStats()

    # plot County C tree
    fig, ax = plt.subplots(figsize=(4, 7), facecolor="w")
    county_c_subtree.plotTree(ax, x_attr=x_attr, colour="black", label="MeV origin")
    county_c_subtree.plotPoints(
        ax, x_attr=x_attr, size=20, colour=color_node, zorder=100
    )

    [ax.spines[loc].set_visible(False) for loc in ["top", "right", "left"]]
    ax.tick_params(axis="y", size=0)
    ax.tick_params(axis="x", labelrotation=90)
    ax.set_yticklabels([])
    ax.xaxis.tick_bottom()
    ax.yaxis.tick_left()
    plt.savefig("county_c.svg")
    # plt.show()

    # Find the County D Co. leaves
    def filter_outbreak_leaves(tree):
        county_d_leaves = []
        for l in tree.Objects:
            if l.is_leaf():
                if l.traits["node_attrs"]["county"]["value"] in ["County D"]:
                    county_d_leaves.append(l)
        return county_d_leaves
    county_d_leaves = filter_outbreak_leaves(tree=tree)

    # Find the MRCA for County D sequences
    mrca_county_d = tree.commonAncestor(county_d_leaves)
    print(f"the mrca node for the County D sequences is {mrca_county_d} \n")
    county_d_subtree = tree.subtree(k=mrca_county_d)
    county_d_subtree.fixHangingNodes()
    print(f"The County D County subtree stats are: \n")
    county_d_subtree.treeStats()

    # plot County D tree
    fig, ax = plt.subplots(figsize=(4, 6.5), facecolor="w")
    county_d_subtree.plotTree(ax, x_attr=x_attr, colour="black")
    county_d_subtree.plotPoints(
        ax, x_attr=x_attr, size=20, colour=color_node, zorder=100
    )

    [ax.spines[loc].set_visible(False) for loc in ["top", "right", "left"]]
    ax.tick_params(axis="y", size=0)
    ax.tick_params(axis="x", labelrotation=90)
    ax.set_yticklabels([])
    ax.xaxis.tick_bottom()
    ax.yaxis.tick_left()
    plt.savefig("county_d.svg")
    #plt.show()


if __name__ == "__main__":
    main()
