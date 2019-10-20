from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace, CircleFace
import csv
from colour import Color


def visualize_tree(csv_path, newick, select_clade = "default"):
    """
    Shows the tree in an ETE window
    :param select_clade: Which WHO clade you are visualising. Defaults to the first clade seen in the csv
    :param csv_path: Information about nodes/tips that need to be displayed.
    :param newick: Phylo tree in newick format.
    :param threshold: Threshold of when to show as green, otherwise red
    :return: null
    """
    results = list(csv.reader(open(csv_path)))

    if select_clade == "default":
        select_clade = results[1][1]

    blue = Color("blue")
    colors = list(blue.range_to(Color("red"), 1000))

    # Load tree with ETE
    t = Tree(newick)

    # Stylize the entire tree
    ts = TreeStyle()
    ts.layout_fn = layout

    # We will add node names manually
    ts.show_leaf_name = False
    # Show branch data
    ts.show_branch_length = True
    ts.show_branch_support = True
    ts.branch_vertical_margin = 99
    ts.scale = 999

    # Define colours
    green = "#00ff00"
    red = "#ff0000"
    grey = "#9c9c9c"

    # i is just a variable used for debugging
    i = 0

    # Traverse through the entire tree, processing all tips
    for n in t.traverse():
        if n.is_leaf():
            isDeleteNode = True
            for result in results:
                if result[2] == n.name:
                    if result[1] != select_clade:
                        break
                    n.add_feature("color", colors[round(float(result[4]) * 999)].get_hex())
                    isDeleteNode = False
                    i += 1
                    break
            if isDeleteNode:
                n.delete()

    t.show(tree_style=ts)

def layout(node):
    if node.is_leaf():
        # Add node name to laef nodes
        N = AttrFace("name", fsize=14, fgcolor="black")
        faces.add_face_to_node(N, node, 0)
        # Creates a sphere face whose size is proportional to node's
        # feature "weight"
        C = CircleFace(radius=250, color=node.color)
        # Let's make the sphere transparent
        C.opacity = 0.3
        # And place as a float face over the tree
        faces.add_face_to_node(C, node, 0, position="float")


visualize_tree("../visualization/result2.csv", "../clade_assignments/trees/flutree2018_5.nwk", select_clade="3c2")
