from ete3 import Tree
from Bio import SeqIO
import sys

def extract_related_sequences(nwk_file, fasta_file, output_file):
    # Laden des phylogenetischen Baums
    tree = Tree(nwk_file)

    # Finden aller 's_all_' Knoten und ihrer Verwandten
    related_sequences = set()
    for node in tree.traverse():
        if node.name and node.name.startswith('s_all_'):
            # Füge den 's_all_' Knoten hinzu
            related_sequences.add(node.name)
            # Füge alle Schwesterknoten hinzu
            for sister in node.get_sisters():
                related_sequences.update(sister.get_leaf_names())

    # Laden und Filtern der FASTA-Sequenzen
    with open(output_file, 'w') as out_handle:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in related_sequences:
                SeqIO.write(record, out_handle, "fasta")

    print(f"Verwandte Sequenzen wurden in {output_file} gespeichert.")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Verwendung: python script.py <nwk_file> <fasta_file> <output_file>")
        sys.exit(1)
    
    nwk_file = sys.argv[1]
    fasta_file = sys.argv[2]
    output_file = sys.argv[3]
    
    extract_related_sequences(nwk_file, fasta_file, output_file)
