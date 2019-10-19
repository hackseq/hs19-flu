from ete3 import Tree
from ete3 import TreeStyle

t = Tree('../clade_assignments/trees/flutree2018_5.nwk')
ts = TreeStyle()
ts.show_branch_length = True
ts.show_branch_support = True
t.show(tree_style=ts)