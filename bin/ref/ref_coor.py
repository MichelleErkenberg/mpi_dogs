import sys
import logging
from Bio import SeqIO
from Bio.Align import PairwiseAligner

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def aln_to_ref(aligned_ref):
    aln_to_ref = {}
    index = 0
    for i, nuc in enumerate(aligned_ref):
        if nuc != "-":
            index += 1
        aln_to_ref[i] = index
    return aln_to_ref

def compare_sequences(ref_seq, query_seq):
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    alignment = aligner.align(ref_seq, query_seq)[0]
    aligned_ref, aligned_query = alignment.target, alignment.query

    coord_map = aln_to_ref(aligned_ref)

    differences = []
    for i, (ref_base, query_base) in enumerate(zip(aligned_ref, aligned_query)):
        if ref_base != query_base and ref_base != '-' and query_base != '-':
            ref_pos = coord_map[i] - 1
            differences.append({
                'ref_pos': ref_pos,
                'ref_base': ref_base,
                'query_base': query_base
            })

    return differences

def process_sequences(ref_file, data_file, output_file):
    try:
        logging.info(f"Reading reference file: {ref_file}")
        ref_seq = next(SeqIO.parse(ref_file, "fasta")).seq

        logging.info(f"Processing data file: {data_file}")
        with open(output_file, 'w') as out_f:
            for record in SeqIO.parse(data_file, "fasta"):
                logging.info(f"Comparing sequence: {record.id}")
                differences = compare_sequences(str(ref_seq), str(record.seq))
                
                out_f.write(f"Differences found for {record.id}:\n")
                for diff in differences:
                    out_f.write(f"  Reference position: {diff['ref_pos']}\n")
                    out_f.write(f"    Reference base: {diff['ref_base']}, Query base: {diff['query_base']}\n")
                out_f.write("\n")
        logging.info(f"Results written to {output_file}")
    except Exception as e:
        logging.error(f"An error occurred: {str(e)}")
        raise

def main(ref_file, data_file, output_file):
    process_sequences(ref_file, data_file, output_file)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        logging.error("Usage: python3 ref_coor.py <reference_file> <data_file> <output_file>")
        sys.exit(1)
    
    ref_file = sys.argv[1]
    data_file = sys.argv[2]
    output_file = sys.argv[3]
    main(ref_file, data_file, output_file)
