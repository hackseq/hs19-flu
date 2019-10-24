![Logo](logo.png)

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

Current status: Theme 1 has successfully completed the WHO clade assignments. Yay! They are moving on to linking those clade assignments to the predictions from Theme 2 so that we can select strains for putative vaccines and then see how our choice of vaccine would perform in our data. 

## Theme 2 - improve on the central model in the preprint 

There are a number of relatively straightforward improvements to make on the code we already used in the preprint. These include using regression instead of classification, modifying the threshold used for defining "success" and expanding previous predictions in to 2019. 

Current status: Theme 2 folks have tested regression, checked the quality of predictions on 2019 data, explored lasso regression (for feature selection) and are linking with Theme 1 re WHO clades and vaccine strain selection. 

Our preduction for the success of the subtrees in 2019 was currect in 70% of the cases when we consider the majority vote between our general model and 10 models from 10-folds cross validation. In 82% of the cases, at least one of our model could curectly predict the success of the subtrees. We also found this paper "Antigenicsites of H1N1 influenza virus hemagglutinin revealed by natural isolates and inhibition assays" which identified the epitope sites of the H1N1 sequences. 

## Theme 3 - explore alternative learning methods 

Maryam suggested exploring whether we could use graph neural networks instead of creating data frames with features that summarise the structure of the subtrees. We think that solving this problem is not feasible during this weekend, but perhaps setting up the problem is. We have team members reading about the required inputs for graph neural networks and thinking through how we would set up our problem in that kind of structure. 

There are two ways to employ graph neural networks to solve this problem. The first approach performs graph-level inference on the phylogenetic subtrees defined previously. The model would learn high level features from the nodes and structure of the graphs and use them to predict the success of the subtree.  The second approach learns higher level representations for each node in the phylogenetic tree (the one, big, global tree) based on the node's features and features of its neighbours. The model then uses these representations to predict the success of the node.

Setting up the first of these structures should be reasonably straightforward: obtain subtrees from Theme 2, represent them as graphs; use the success outcomes from Theme 2 and train GNNs to predict these outcomes. It would be pretty straightforward, but a bit more time-consuming, to develop a suite of features to place on the nodes of the subtrees. These should be derived from the NA and HA sequences for the tips in the subtree -- for example, the epitope features would be a great starting point. 

## Theme 4 - use simpler models and larger clades of trees 

Here, instead of dividing the phylogeny into many small subtrees of size approximately 10-35 we will focus on fewer, larger clades. Now, we won't have sufficient data to train machine learning models, but we can use the number of lineages in the trees per unit time and other simple summaries of tree growth to see if we can build a predictive model that signifies whether a clade is going to grow into the future. 

Current status: We have functions that extract larger clades from the global phylogeny and profile how these have evolved over time. These profiles look at the diversification of the tree in time and in genetic distance space. We have also implemented a definition of "success" in this broader context. The next step is to combine visualisations and get intuition for whether there are signatures of which clades will be highly successful. 


