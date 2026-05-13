import os
import argparse
from Bio import SeqIO
from Bio.PDB import PDBParser, MMCIFParser, PPBuilder

def extract_seqs_from_structure_files(input_folder, output_file):
    pdb_parser = PDBParser(QUIET=True)
    cif_parser = MMCIFParser(QUIET=True)
    pp_builder = PPBuilder()  # Extract sequences

    with open(output_file, 'w') as fasta_out:
        for structure_file in os.listdir(input_folder):
            if structure_file.endswith('.pdb') or structure_file.endswith('.cif'):
                file_name, file_extension = os.path.splitext(structure_file)
                file_path = os.path.join(input_folder, structure_file)
                try:
                    # Choose the correct parser based on file extension
                    if file_extension == '.pdb':
                        structure = pdb_parser.get_structure(file_name, file_path)
                    elif file_extension == '.cif':
                        structure = cif_parser.get_structure(file_name, file_path)
                    else:
                        continue

                    # Extract sequences
                    for model in structure:
                        for chain in model:
                            peptides = pp_builder.build_peptides(chain)
                            if peptides:
                                sequence = ''.join([str(pp.get_sequence()) for pp in peptides])
                                chain_id = chain.id
                                fasta_header = f">{file_name}_{chain_id}"
                                fasta_out.write(f"{fasta_header}\n{sequence}\n")
                            else:
                                print(f"Warning: No sequence found for {file_name}_{chain_id}")
                except Exception as e:
                    print(f"Error processing {structure_file}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract amino acid sequences from PDB or CIF files and save in FASTA format.")
    parser.add_argument("input_folder", type=str, help="Path to the folder containing PDB or CIF files.")
    parser.add_argument("output_file", type=str, help="Path to the output FASTA file.")

    args = parser.parse_args()

    extract_seqs_from_structure_files(args.input_folder, args.output_file)
