from ete3 import Tree, NodeStyle, TreeStyle
import random
import csv
results = list(csv.reader(open("../visualization/result1.csv")))


t = Tree("../clade_assignments/trees/flutree2018_5.nwk")

# Basic tree style
ts = TreeStyle()
ts.show_leaf_name = True
ts.branch_vertical_margin = 99
ts.scale = 9999
# Creates an independent node style for each node, which is
# initialized with a red foreground color.

green = "#00ff00"
red = "#ff0000"
grey = "#9c9c9c"
i = 0
for n in t.traverse():
    if n.is_leaf():
        tipPredicted = False;
        for result in results:
            if result[1] == n.name:
                tipPredicted = True;
                if float(result[2]) > 0.8:
                    color = green
                else:
                    color = red
                nstyle = NodeStyle()
                nstyle["fgcolor"] = color
                nstyle["size"] = 500
                n.set_style(nstyle)
                i += 1;
                break
        if not tipPredicted:
            n.delete()
    else:
        for result in results:
            if result[1] == n.name:
                tipPredicted = True;
                if float(result[2]) > 0.8:
                    color = green
                else:
                    color = red
                nstyle = NodeStyle()
                nstyle["fgcolor"] = color
                nstyle["size"] = 99
                n.set_style(nstyle)
                i += 1;
                break



t.show(tree_style=ts)
