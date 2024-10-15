import sys
import csv
from Bio import AlignIO
from Bio import SeqIO

def process_alignment(ref_file, aln_file, output_file):
    # Lese die Referenzsequenz
    ref_record = next(SeqIO.parse(ref_file, "fasta"))
    ref_seq = str(ref_record.seq)

    # Lese das Alignment
    alignment = AlignIO.read(aln_file, "fasta")

    # Erstelle ein Dictionary für die Sequenzen
    sequences = {record.id: str(record.seq) for record in alignment}

    # Öffne die Ausgabedatei
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        # Schreibe die Kopfzeile
        header = ['Position', 'Reference'] + list(sequences.keys())
        writer.writerow(header)

        # Schreibe die Daten
        for i, ref_base in enumerate(ref_seq, start=1):
            row = [i, ref_base]
            for seq_id in sequences:
                if i <= len(sequences[seq_id]):
                    row.append(sequences[seq_id][i-1])
                else:
                    row.append('')
            writer.writerow(row)

    print(f"Output written to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <reference_file> <alignment_file> <output_file>")
        sys.exit(1)

    ref_file = sys.argv[1]
    aln_file = sys.argv[2]
    output_file = sys.argv[3]

    process_alignment(ref_file, aln_file, output_file)
