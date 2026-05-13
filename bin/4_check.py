from pathlib import Path


def check_files():
    # 1. 定义路径 (Path.cwd() 获取当前路径)
    base_path = Path.cwd() / "4Multimer_prediction"
    afm_dir = base_path / "nanopore_AFM"
    af3_dir = base_path / "nanopore_AF3"
    id_file = Path("3_nanopore_allminingF_id.txt")

    # 2. 检查目录是否存在，不存在则创建
    for d in [afm_dir, af3_dir]:
        if not d.exists():
            d.mkdir(parents=True)
            print(f"Created: {d}\nPlease put files in it.")

    if not id_file.exists():
        print(f"Error: {id_file} not found.")
        return

    # 3. 读取 ID 并检查
    # read_text().splitlines() 自动去除换行符并过滤空行
    ids = [line.strip() for line in id_file.read_text().splitlines() if line.strip()]

    error_log = []

    for _id in ids:
        # 定义每个 ID 对应的 4 个预期文件
        expected_files = [
            afm_dir / f"{_id}.pdb",
            afm_dir / f"{_id}.json",
            af3_dir / f"{_id}.cif",
            af3_dir / f"{_id}.json"
        ]

        # 找出缺失的文件名
        missing = [f.name for f in expected_files if not f.exists()]

        if missing:
            error_log.append(f"ID {_id}: Missing -> {', '.join(missing)}")

    # 4. 汇总输出
    if not error_log:
        print(f"✅ Pass! All {len(ids)} IDs are consistent.")
    else:
        print(f"❌ Failed! Found {len(error_log)} inconsistent IDs:")
        print("\n".join(error_log))


if __name__ == "__main__":
    check_files()