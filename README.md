# ARNG
Automated Reaction Network Generation (ARNG) methodology is designed to automatize application of chemical reactions to the molecules. Current sample of code automatizes reconstruction of the reaction network of ethyl linoleate polymerizing at the room temperature. The reaction network is further transformed into a kinetic model - system of ODEs describing how the concentration of molecular species changes in time. 

INPUT: reactant molecules (defined as molecular graphs) that are needed to start the reaction process, library of reaction rules for a given reaction mechanism, initial concentrations of the starting reactant molecules, time 

The computational methodology is split into 4 main parts, which have to be consecutively executed:
- main_ARNG.m - takes the input, initialized all main variables that are going to be updated used through all 4 parts of the code. This part of the code consecutively applies the transformations to the input molecular graphs (the reactants which start the reaction process). The methodology creates more intermediate species which react further till no unreacted species are left. The output of this part is a list of all possible reactions that can happen between allpossible molecular structures present in the chemical system. 

- main_Reaction_Network.m - starts with the list of all the reactions and transforms this list into a reaction network. Reaction network is a bipartite graph with two types of nodes, reactions and molecules. Reaction network captures all the interactions between molecules in a systematic way. Reaction network is an output of this section.

- main_Rate_Parameters.m - estimates rate parameters of all the reactions in the reaction network using Arhenius equation. Pre-exponential factor and the intrinsic energy are taken from the literature, while the activation energy is estimated using Bell–Evans–Polanyi approximation.

- main_Kinetic_Model.m - this part of the code takes a reaction network as an input and automatically transforms it into a kinetic model. Kinetic model is a system of differential equations describing the temporal changes in concentration of all the generated molecular species. The ODEs are solved numerically and the result can be visualized in various way depending on the needs of the user. 

Some terminology:

- Molecules are represented as molecular graphs and defined by their adjacency matrix, which captures the connectivity between atoms, and list of labels, which defines a chemical element of the atom. 

- Patterns are subgraphs of molecular graphs, which correspond to the reactive sites of the molecules. 

- Reaction rules act on patterns of reactant molecules and transform them into patterns of product molecules.  

- By repeatedly applying reaction rules to the molecular graphs, we are able to systematically uncover the reaction network of the chemical system.

- Reaction network is a bipartite graph that describes interactions between molecules. 

