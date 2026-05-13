import os
import sys
import shutil
from Bio import PDB

def average_plddt(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    plddt_scores = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    plddt_scores.append(residue['CA'].get_bfactor())
    if plddt_scores:
        return sum(plddt_scores) / len(plddt_scores)
    else:
        raise ValueError(f"No pLDDT scores found in {pdb_file}")


def main():
    threshold = 70
    if len(sys.argv) < 2:
        print("Usage: python 1_clean_plddt70.py <input_directory>")
        sys.exit(1)
    raw_input_dir =  os.path.abspath(os.path.normpath(sys.argv[1]))
    base_name = os.path.basename(os.path.abspath(raw_input_dir))
    parent_dir = os.path.dirname(raw_input_dir)
    outdir = os.path.join(parent_dir, f'{base_name}_t{threshold}')
    outfile = os.path.join(parent_dir, f'0avg_plddt_{base_name}.txt')
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    with open(outfile, "w") as f:
        f.write("Filename\tAverage_plddt\n")

        for filename in os.listdir(raw_input_dir):
            if filename.endswith(".pdb"):
                pdb_path = os.path.join(raw_input_dir, filename)
                try:
                    avg_plddt = average_plddt(pdb_path)
                    f.write(f"{filename}\t{avg_plddt:.2f}\n")
                    if avg_plddt >= threshold:
                        shutil.copy(pdb_path, outdir)
                    else:
                        print(f"{filename} {avg_plddt}")
                except Exception as e:
                    print(f"Error processing {filename}: {e}")
if __name__ == '__main__':
    main()
