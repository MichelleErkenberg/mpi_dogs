from ete3 import Tree
from Bio import SeqIO
import sys
from collections import defaultdict

def get_immediate_relatives(node):
    relatives = set()
    
    # Parent
    parent = node.up
    if parent:
        relatives.add(parent.name)
        
        # Siblings and their children (nephews/nieces)
        for sibling in parent.children:
            if sibling != node:
                relatives.add(sibling.name)
                relatives.update(child.name for child in sibling.children)
        
        # Aunts/Uncles and their children (cousins)
        grandparent = parent.up
        if grandparent:
            for aunt_uncle in grandparent.children:
                if aunt_uncle != parent:
                    relatives.add(aunt_uncle.name)
                    relatives.update(cousin.name for cousin in aunt_uncle.children)
    
    # Children
    relatives.update(child.name for child in node.children)
    
    return relatives

def extract_and_organize_sequences(nwk_file, alignment_file, output_file):
    # Load the phylogenetic tree
    tree = Tree(nwk_file)

    # Dictionary to store 's_all_' sequences and their related sequences
    sequence_groups = defaultdict(set)

    # Find all 's_all_' nodes and their immediate relatives
    for node in tree.traverse():
        if node.name and node.name.startswith('s_all_'):
            relatives = get_immediate_relatives(node)
            sequence_groups[node.name] = relatives.union({node.name})

    # Load all sequences from the alignment file
    all_sequences = {record.id: record for record in SeqIO.parse(alignment_file, "fasta")}

    # Write organized sequences to output file
    with open(output_file, 'w') as out_handle:
        for s_all_seq in sorted(sequence_groups.keys()):
            if s_all_seq in all_sequences:
                SeqIO.write(all_sequences[s_all_seq], out_handle, "fasta")
            
            related_seqs = sorted(sequence_groups[s_all_seq] - {s_all_seq})
            for seq_id in related_seqs:
                if seq_id in all_sequences:
                    SeqIO.write(all_sequences[seq_id], out_handle, "fasta")

    print(f"Organized sequences have been saved to {output_file}.")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <nwk_file> <alignment_file> <output_file>")
        sys.exit(1)
    
    nwk_file = sys.argv[1]
    alignment_file = sys.argv[2]
    output_file = sys.argv[3]
    
    extract_and_organize_sequences(nwk_file, alignment_file, output_file)
