import os
import argparse
import gemmi

def extract_plddt_from_structure_first_chain(file_path):
    """只从第一条链中提取 pLDDT 值并计算平均值。"""
    try:
        structure = gemmi.read_structure(file_path)
        plddt_values = []
        # 只处理 model[0] 的第一条链
        model = structure[0]
        if len(model) == 0:
            print(f"警告：文件 {file_path} 中没有任何链")
            return None
        first_chain = model[0]
        for residue in first_chain:
            atom_ca = None
            for atom in residue:
                if atom.name.upper() == 'CA':
                    atom_ca = atom
                    break
                # gemmi 里 pLDDT 通常存储在 b_iso 字段
            if atom_ca and atom_ca.b_iso > 0:
                plddt_values.append(atom_ca.b_iso)
            else :
                pass
        if plddt_values:
            return sum(plddt_values) / len(plddt_values)
        else:
            print(f"警告：在文件 {file_path} 的第一条链中未找到有效的 CA pLDDT 值。")
            return None
    except Exception as e:
        print(f"处理文件 {file_path} 时出错：{e}")
        return None

def process_structure_files(input_folder, output_file):
    """处理文件夹中的结构文件，计算每个文件第一条链的平均 pLDDT，并写入输出文件。"""
    with open(output_file, 'w') as out_f:
        for filename in os.listdir(input_folder):
            if filename.endswith(('.cif', '.pdb')):
                file_path = os.path.join(input_folder, filename)
                avg_plddt = extract_plddt_from_structure_first_chain(file_path)
                if avg_plddt is not None:
                    out_f.write(f"{filename}\t{avg_plddt:.2f}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="计算文件夹中每个 .cif 或 .pdb 文件第一条链的平均 pLDDT 值。"
    )
    parser.add_argument("input_folder",  help="包含 .cif 或 .pdb 文件的文件夹路径")
    parser.add_argument("output_file",   help="写入结果的输出文件路径")
    args = parser.parse_args()
    process_structure_files(args.input_folder, args.output_file)
