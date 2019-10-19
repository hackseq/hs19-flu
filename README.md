# hs19-flu

## Summary and outline
Project summary: We aim to predict which variants of influenza virus will succeed and circulate in the next 1-2 years. 
We have previously established a "proof of principle" that we can do this task fairly well using phylogenetic trees reconstructed from influenza HA sequences.

### Previous approach
Briefly, in that work, we selected subtrees from large-scale phylogenies. These were of varying size (10-35). Each subtree contains a set of closely-related flu sequences. We looked to see whether each subtree grew into the future, removing subtrees that arose so recently that we don't know yet whether they grew into their future. For example, if a subtree only started in April 2018 and our original dataset ended in May 2018 then we don't know whether the subtree was successful. Success was categorical: 1 if the size of the subtree increased by a factor of at least 1.1 (for example) in a fixed time frame, and 0 otherwise. 

Next, we computed a suite of features of each subtree, seen in its own right as a phylogenetic tree. Features included tree imbalance (measures of asymmetry in the tree), features derived from network science (centrality, spectral features) and features derived from the trees' branch lengths. We also computed a few simple feature of the HA sequences corresponding to the tips. 

We then trained machine learning methods to distinguish successful from unsuccessful trees, holding back some data to test the predictions. We found a reasonably good accuracy, did various comparisons to previous work, and performed a range of testing (for exmaple, training on one strain of influenza and testing on others). 

Here, we are expanding that work in 3-4 themes. 

## Theme 1 - design vaccines with our methods, and compare them to WHO-selected seasonal vaccines

Align HA sequences from human influenza A H3N2, influenza A H1N1 and influenza B to reference sequences provided by nextflu.org and the 'augur' pipeline. 

Next, use those alignments together with phylogenetic trees to determine which WHO-defined clade each sequence belongs to. 

This is a key next step for the flu prediction project, because it will allow us to compare our predictions to the current and recent seasonal influenza vaccines and determine whether our designed vaccines would have been a good match to circulating strains in the next season.

## Theme 2 - improve on the central model in the preprint 

There are a number of relatively straightforward improvements to make on the code we already used in the preprint. These include using regression instead of classification, modifying the threshold used for defining "success" and expanding previous predictions in to 2019. 

## Theme 3 - explore alternative learning methods 

Maryam suggested exploring whether we could use graph neural networks instead of creating data frames with features that summarise the structure of the subtrees. We think that solving this problem is not feasible during this weekend, but perhaps setting up the problem is. We have team members reading about the required inputs for graph neural networks and thinking through how we would set up our problem in that kind of structure. 

There are two ways to employ graph neural networks to solve this problem.First approach performs graph level inference on the phylogenetic subtrees defined previously. The model would learn high level features from the nodes and structure of the graph and use them to predict the success of the subtree.
The second approach learns higher level representations for each node in the phylogenetic tree based on its features and features of its neighbours. The model then uses these representations to predict the success of the node.

## Theme 4 - use simpler models and larger clades of trees 

Here, instead of dividing the phylogeny into many small subtrees of size approximately 10-35 we will focus on fewer, larger clades. Now, we won't have sufficient data to train machine learning models, but we can use the number of lineages in the trees per unit time and other simple summaries of tree growth to see if we can build a predictive model that signifies whether a clade is going to grow into the future. 


