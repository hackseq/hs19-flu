from ete3 import Tree, NodeStyle, TreeStyle

t = Tree("../clade_assignments/trees/flutree2018_5.nwk")

# Basic tree style
ts = TreeStyle()
ts.show_leaf_name = True
ts.scale =  120
# Creates an independent node style for each node, which is
# initialized with a red foreground color.
r = 0;
g = 0;
b = 0;
for n in t.traverse():
    if(n.is_leaf()):
        print(n.name)
        print(n.dist)
        nstyle = NodeStyle()
        print("#" + str(hex(r))[2:]+str(hex(g))[2:]+str(hex(b))[2:])
        nstyle["fgcolor"] ="#" + str(hex(r))[2:]+str(hex(g))[2:]+str(hex(b))[2:]
        nstyle["size"] = 15
        n.set_style(nstyle)
        if(r < 255):
            r +=1
        elif(g < 255):
            g += 1
        elif(b < 255):
            b +=1
        elif(b >= 255):
            r = 0
            g = 0
            b = 0

# Let's now modify the aspect of the root node
t.img_style["size"] = 30
t.img_style["fgcolor"] = "blue"

t.show(tree_style=ts)
