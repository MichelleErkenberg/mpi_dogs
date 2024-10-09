from ete3 import Tree
from Bio import SeqIO
import sys
from collections import defaultdict

def get_immediate_relatives(node, exclude_names):
    """
    Find immediate relatives of a given node in the phylogenetic tree.
    This includes parent, siblings, children, aunts/uncles, and cousins.
    Excludes names present in exclude_names set.
    """
    relatives = set()
    
    # Check for parent
    parent = node.up
    if parent and parent.name not in exclude_names:
        relatives.add(parent.name)
        
        # Add siblings and their children (nieces/nephews)
        for sibling in parent.children:
            if sibling != node and sibling.name not in exclude_names:
                relatives.add(sibling.name)
                relatives.update(child.name for child in sibling.children if child.name not in exclude_names)
        
        # Check for grandparent
        grandparent = parent.up
        if grandparent:
            # Add aunts/uncles and their children (cousins)
            for aunt_uncle in grandparent.children:
                if aunt_uncle != parent and aunt_uncle.name not in exclude_names:
                    relatives.add(aunt_uncle.name)
                    relatives.update(cousin.name for cousin in aunt_uncle.children if cousin.name not in exclude_names)
    
    # Add children
    relatives.update(child.name for child in node.children if child.name not in exclude_names)
    
    return relatives

def extract_and_organize_sequences(nwk_file, alignment_file, names_file, output_file):
    """
    Main function to extract and organize sequences based on phylogenetic relationships.
    """
    # Load the phylogenetic tree
    tree = Tree(nwk_file)

    # Dictionary to store sequences and their related sequences
    sequence_groups = defaultdict(set)

    # Read names from the names file
    with open(names_file, 'r') as f:
        names = [line.strip() for line in f]
    
    # Create a set of names to exclude
    exclude_names = set(names)

    # Find immediate relatives for each name
    for name in names:
        node = tree.search_nodes(name=name)
        if node:
            node = node[0]  # Take the first match if there are multiple
            relatives = get_immediate_relatives(node, exclude_names)
            sequence_groups[name] = relatives.union({name})

    # Load all sequences from the alignment file
    all_sequences = {record.id: record for record in SeqIO.parse(alignment_file, "fasta")}

    # Write organized sequences to output file
    with open(output_file, 'w') as out_handle:
        for main_seq in sorted(sequence_groups.keys()):
            if main_seq in all_sequences:
                SeqIO.write(all_sequences[main_seq], out_handle, "fasta")
            
            related_seqs = sorted(sequence_groups[main_seq] - {main_seq})
            for seq_id in related_seqs:
                if seq_id in all_sequences:
                    SeqIO.write(all_sequences[seq_id], out_handle, "fasta")

    print(f"Organized sequences have been saved to {output_file}.")

if __name__ == "__main__":
    # Check if the correct number of command-line arguments is provided
    if len(sys.argv) != 5:
        print("Usage: python script.py <nwk_file> <alignment_file> <names_file> <output_file>")
        sys.exit(1)
    
    # Assign command-line arguments to variables
    nwk_file = sys.argv[1]
    alignment_file = sys.argv[2]
    names_file = sys.argv[3]
    output_file = sys.argv[4]
    
    # Run the main function
    extract_and_organize_sequences(nwk_file, alignment_file, names_file, output_file)
