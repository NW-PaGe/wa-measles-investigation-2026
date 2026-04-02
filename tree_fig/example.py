import matplotlib as mpl
from matplotlib import pyplot as plt

import requests
from io import StringIO as sio

import baltic as bt

address='https://raw.githubusercontent.com/evogytis/fluB/master/data/mcc%20trees/InfB_NPt_ALLs1.mcc.tre' ## address of example tree
fetch_tree = requests.get(address) ## fetch tree
treeFile=sio(fetch_tree.text) ## stream from repo copy

ll=bt.loadNexus(treeFile,tip_regex='_([0-9\-]+)$') ## treeFile here can alternatively be a path to a local file
ll.treeStats() ## report stats about tree
print(ll.Objects[4].traits)
fig,ax = plt.subplots(figsize=(10,15),facecolor='w')

traitName='PB1'
tipSize=20

cumulative_y=0
keep_tips=[k for k in ll.Objects if k.is_leaf() and k.traits[traitName]=='V'] ## list of leaf objects that will remain in the tree - here only those whose PB1 is Victoria lineage
reduced_tree=ll.reduceTree(keep_tips) ## retrieve a reduced tree (multitype by default)

x_attr=lambda k: k.absoluteTime
c_func=lambda k: 'indianred' if k.traits[traitName]=='V' else 'steelblue' ## Victoria lineage red, Yamagata blue

reduced_tree.plotTree(ax,x_attr=x_attr,colour=c_func)
reduced_tree.plotPoints(ax,x_attr=x_attr,size=tipSize,colour=c_func,zorder=100)

ax.xaxis.tick_bottom()
ax.yaxis.tick_left()

[ax.spines[loc].set_visible(False) for loc in ['top','right','left']]

ax.tick_params(axis='y',size=0)
ax.set_yticklabels([])
ax.set_ylim(-5,reduced_tree.ySpan+5)

plt.show()