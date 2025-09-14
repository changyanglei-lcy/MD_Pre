#! / usr / bin / python3
import os

def replace_resname_in_mol2_inplace(file_path):
    base_name = os.path.splitext(os.path.basename(file_path))[0]

    lines_out = []
    in_atom_section = False
    with open(file_path, 'r') as f_in:
        for line in f_in:
            stripped = line.strip()
            if stripped == "@<TRIPOS>ATOM":
                in_atom_section = True
                lines_out.append(line)
                continue
            if stripped.startswith("@<TRIPOS>") and stripped != "@<TRIPOS>ATOM":
                in_atom_section = False
                lines_out.append(line)
                continue

            if in_atom_section:
                parts = line.rstrip('\n').split()
                if len(parts) >= 9:
                    parts[-2] = base_name  # 替换残基名为文件名
                    line = "  ".join(parts) + "\n"
                lines_out.append(line)
            else:
                lines_out.append(line)

    with open(file_path, 'w') as f_out:
        f_out.writelines(lines_out)

    print(f"[成功] 修改完成：{file_path} （残基名改为 {base_name}）")

def recursive_modify_all_mol2_files(root_dir="."):
    mol2_count = 0
    for dirpath, _, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename.lower().endswith(".mol2"):
                file_path = os.path.join(dirpath, filename)
                try:
                    replace_resname_in_mol2_inplace(file_path)
                    mol2_count += 1
                except Exception as e:
                    print(f"[!] 处理失败：{file_path}，原因：{e}")

    if mol2_count == 0:
        print("警告：没有找到任何 .mol2 文件。")
    else:
        print(f"\n成功：总共处理了 {mol2_count} 个 mol2 文件。")

if __name__ == "__main__":
    recursive_modify_all_mol2_files(".")