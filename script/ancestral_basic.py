import pandas as pd
from Bio import Phylo
from collections import defaultdict

def load_matrix(path):
    return pd.read_csv(path, index_col=0)

def load_tree(path):
    return next(Phylo.parse(path, "newick"))


### Load trees 
matrix_path = "/Users/hugo/Desktop/ancestral/gene_copy_number_matrix.csv"        
tree_path = "/Users/hugo/Desktop/ancestral/combined_sequences_aln.tree"  

# Load files
matrix = load_matrix(matrix_path)
tree = load_tree(tree_path)

print("\nMatrix with orthologues genes:")
print(matrix)

print("\nTree summary:")
Phylo.draw_ascii(tree)

print(f"\nTotal terminal nodes: {len(tree.get_terminals())}")
print(f"\nTotal internal nodes: {len(tree.get_nonterminals())}")
print(f"\nTree root: {tree.root}")
print(f"\nRoot has {len(tree.root.clades)} children")


#  Parsimony of Fitch 

def run_fitch_algorithm(matrix, tree):
    fitch = defaultdict(dict)
    states = defaultdict(dict)

    def postorder(node):
        if node.is_terminal():
            for gene in matrix.columns:
                val = int(matrix.loc[node.name, gene])
                fitch[node][gene] = {val}
                states[node][gene] = val
        else:
            for child in node.clades:
                postorder(child)
            for gene in matrix.columns:
                child_sets = [fitch[child][gene] for child in node.clades]
                intersect = set.intersection(*child_sets)
                fitch[node][gene] = intersect if intersect else set.union(*child_sets)

    def preorder(node, parent_state=None):
        for gene in matrix.columns:
            options = fitch[node][gene]
            if parent_state and parent_state[gene] in options:
                states[node][gene] = parent_state[gene]
            else:
                states[node][gene] = min(options)
        for child in node.clades:
            preorder(child, states[node])

    postorder(tree.root)
    preorder(tree.root)
    return states

# Run the parsimony 
state_matrix = run_fitch_algorithm(matrix, tree)
ancestral = pd.Series(state_matrix[tree.root], name="Ancestral (root)")

print("\n Ancestral genome, node root is:")
print(ancestral)