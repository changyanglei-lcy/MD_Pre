#!/usr/bin/env python3
"""
GROMACS模拟文件自动生成脚本
批量处理分子对并生成GROMACS分子动力学模拟所需的完整文件集

作者：Claude Code
版本：1.0
"""

import os
import sys
import csv
import requests
import subprocess
import shutil
import logging
import time
from datetime import datetime
from pathlib import Path
from tqdm import tqdm
import pandas as pd
import numpy as np


class GromacsProcessor:
    """GROMACS模拟文件自动处理器"""

    def __init__(self, mol_csv_path="Mol.csv", file_dir="File", log_level=logging.INFO):
        """初始化处理器"""
        self.mol_csv_path = mol_csv_path
        self.file_dir = file_dir
        self.base_dir = Path(".")

        # 样本数据
        self.samples = []
        self.completed_samples = []
        self.failed_samples = []
        self.stats = {
            'total': 0,
            'completed': 0,
            'failed': 0,
            'skipped': 0
        }

        # 创建日志目录
        self.log_dir = Path("logs")
        self.log_dir.mkdir(exist_ok=True)

        # 设置日志
        self.setup_logging(log_level)

        # 检查依赖
        self.check_dependencies()

    def setup_logging(self, level):
        """设置日志系统"""
        log_file = self.log_dir / f"gromacs_processor_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

        logging.basicConfig(
            level=level,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file, encoding='utf-8'),
                logging.StreamHandler(sys.stdout)
            ]
        )

        self.logger = logging.getLogger(__name__)
        self.logger.info(f"日志文件：{log_file}")

    def check_dependencies(self):
        """检查依赖工具"""
        self.logger.info("检查依赖工具...")

        required_tools = [
            ('obabel', 'Open Babel'),
            ('acpype', 'ACPype'),
            ('conda', 'Conda')
        ]

        missing_tools = []
        for tool, name in required_tools:
            if not shutil.which(tool):
                missing_tools.append(f"{name} ({tool})")

        if missing_tools:
            self.logger.warning(f"缺少以下工具（脚本仍可继续运行但部分功能会受限）：{', '.join(missing_tools)}")
            self.logger.info("请确保已安装：Open Babel、ACPype、Conda")
        else:
            self.logger.info("所有依赖工具检查通过")

    def read_mol_csv(self):
        """读取Mol.csv文件并验证数据"""
        self.logger.info("读取Mol.csv文件...")

        try:
            # 尝试不同的编码方式读取CSV
            try:
                df = pd.read_csv(self.mol_csv_path, encoding='utf-8-sig')
            except UnicodeDecodeError:
                df = pd.read_csv(self.mol_csv_path, encoding='utf-8')

            self.logger.info(f"读取到 {len(df)} 个样本数据")

            # 验证数据列（支持小写和大写）
            required_columns = ['sample', 'CID_A', 'CID_B']
            missing_columns = [col for col in required_columns if col not in df.columns]

            if missing_columns:
                # 尝试大小写不敏感匹配
                df_columns_lower = [col.lower() for col in df.columns]
                required_columns_lower = [col.lower() for col in required_columns]

                matched_columns = []
                for req_col in required_columns_lower:
                    for df_col in df.columns:
                        if df_col.lower() == req_col:
                            matched_columns.append(df_col)
                            break

                if len(matched_columns) == len(required_columns):
                    # 重命名列以匹配标准格式
                    column_mapping = {}
                    for req_col, match_col in zip(required_columns, matched_columns):
                        column_mapping[match_col] = req_col

                    df = df.rename(columns=column_mapping)
                    self.logger.info(f"列名已映射：{column_mapping}")
                else:
                    error_msg = f"CSV文件缺少必需列：{missing_columns}"
                    self.logger.error(error_msg)
                    raise ValueError(error_msg)

            # 检查数据完整性
            valid_samples = []
            for idx, row in df.iterrows():
                sample_id = row['sample']
                cid_a = row['CID_A']
                cid_b = row['CID_B']

                if pd.isna(sample_id) or pd.isna(cid_a) or pd.isna(cid_b):
                    self.logger.warning(f"第 {idx+2} 行数据不完整，跳过")
                    continue

                if not isinstance(cid_a, (int, float, np.int64)) or not isinstance(cid_b, (int, float, np.int64)):
                    self.logger.warning(f"第 {idx+2} 行CID格式错误：CID_A={cid_a} (type: {type(cid_a)}), CID_B={cid_b} (type: {type(cid_b)})")
                    continue

                valid_samples.append({
                    'Sample': int(sample_id),
                    'CID_A': int(cid_a),
                    'CID_B': int(cid_b)
                })

            self.logger.info(f"有效样本数据：{valid_samples}")

            self.samples = valid_samples
            self.stats['total'] = len(self.samples)
            self.logger.info(f"有效样本数量：{self.stats['total']}")

        except Exception as e:
            self.logger.error(f"读取Mol.csv失败：{e}")
            raise

    def check_completed_samples(self):
        """检查已完成的样本"""
        self.logger.info("检查已完成的样本...")

        completed_count = 0
        for sample in self.samples:
            sample_dir = self.base_dir / str(sample['Sample'])
            if sample_dir.exists():
                # 检查关键文件是否存在
                key_files = [
                    sample_dir / 'MOA_GMX.gro',
                    sample_dir / 'MOA_GMX.itp',
                    sample_dir / 'MOB_GMX.gro',
                    sample_dir / 'MOB_GMX.itp'
                ]

                if all(f.exists() for f in key_files):
                    self.completed_samples.append(sample['Sample'])
                    completed_count += 1
                    continue

        self.stats['skipped'] = completed_count
        self.logger.info(f"已跳过 {completed_count} 个已完成样本")

    def download_pubchem_structure(self, cid, output_path, sample_dir):
        """从PubChem下载分子结构"""
        self.logger.info(f"下载 CID {cid} 的结构...")

        # 尝试下载3D结构
        url_3d = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"

        try:
            response = requests.get(url_3d, timeout=30)
            if response.status_code == 200:
                with open(output_path, 'wb') as f:
                    f.write(response.content)
                self.logger.info(f"成功下载3D结构：{output_path}")
                return True
        except Exception as e:
            self.logger.warning(f"3D结构下载失败：{e}")

        # 如果3D失败，尝试下载2D结构
        url_2d = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF"

        try:
            response = requests.get(url_2d, timeout=30)
            if response.status_code == 200:
                with open(output_path, 'wb') as f:
                    f.write(response.content)
                self.logger.info(f"成功下载2D结构：{output_path}")
                return True
        except Exception as e:
            self.logger.error(f"2D结构下载失败：{e}")

        return False

    def convert_sdf_to_mol2(self, sdf_path, mol2_path):
        """转换SDF格式到MOL2格式"""
        self.logger.info(f"转换 {sdf_path} 到 MOL2 格式...")

        try:
            cmd = f"obabel {sdf_path} -O {mol2_path}"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=60)

            if result.returncode != 0:
                self.logger.error(f"SDF转MOL2失败：{result.stderr}")
                return False

            # 删除SDF文件
            os.remove(sdf_path)
            self.logger.info(f"成功转换并删除SDF文件：{mol2_path}")
            return True

        except Exception as e:
            self.logger.error(f"转换过程中出现错误：{e}")
            return False

    def optimize_molecule_geometry(self, mol2_path, force_field="MMFF94"):
        """优化分子几何结构"""
        self.logger.info(f"优化分子几何结构：{mol2_path}")

        try:
            # 使用Open Babel进行几何优化
            cmd = f"obabel {mol2_path} -O {mol2_path} --minimize --ff {force_field} --steps 1000 --dielectric 78.0"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=120)

            if result.returncode != 0:
                self.logger.error(f"几何优化失败：{result.stderr}")
                return False

            self.logger.info(f"成功优化分子几何结构：{mol2_path}")
            return True

        except Exception as e:
            self.logger.error(f"几何优化过程中出现错误：{e}")
            return False

    def run_replace_resname(self, sample_dir):
        """运行残基重命名脚本"""
        self.logger.info(f"运行残基重命名脚本：{sample_dir}")

        try:
            # 切换到工作目录
            original_dir = os.getcwd()
            os.chdir(self.base_dir)

            cmd = f"python replace_resname.py {sample_dir}"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=60)

            if result.returncode != 0:
                self.logger.error(f"残基重命名失败：{result.stderr}")
                return False

            os.chdir(original_dir)
            self.logger.info("成功完成残基重命名")
            return True

        except Exception as e:
            self.logger.error(f"残基重命名过程中出现错误：{e}")
            return False

    def generate_amber_topology(self, mol2_path, basename, charge=0):
        """生成AMBER拓扑"""
        self.logger.info(f"生成AMBER拓扑：{mol2_path}")

        try:
            # 切换到样本目录
            original_dir = os.getcwd()
            mol_dir = mol2_path.parent
            os.chdir(mol_dir)

            # 激活conda环境（使用conda run命令）
            cmd = f"conda run -n gcc acpype -i {mol2_path.name} -c bcc -a gaff2 -n {charge} -b {basename}"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=3600)

            os.chdir(original_dir)

            if result.returncode != 0:
                self.logger.error(f"AMBER拓扑生成失败：{result.stderr}")
                return False

            self.logger.info(f"成功生成AMBER拓扑：{basename}")
            return True

        except Exception as e:
            self.logger.error(f"AMBER拓扑生成过程中出现错误：{e}")
            return False

    def extract_atomtypes(self, itp_file, output_file):
        """提取原子类型参数（剪切功能）"""
        self.logger.info(f"提取原子类型参数：{itp_file}")

        try:
            with open(itp_file, 'r') as f:
                lines = f.readlines()

            # 读取原文件内容
            original_content = lines.copy()

            # 提取 atomtypes 段
            atomtypes_content = []
            in_atomtypes = False
            atomtypes_start = -1
            atomtypes_end = -1

            for i, line in enumerate(lines):
                if line.strip() == '[ atomtypes ]':
                    in_atomtypes = True
                    atomtypes_start = i
                    atomtypes_content.append(line)
                elif in_atomtypes:
                    if line.strip() == '':  # 空行表示段落结束
                        atomtypes_end = i + 1
                        break
                    atomtypes_content.append(line)

            # 如果找到 atomtypes 段，写入输出文件并从原文件删除
            if atomtypes_start >= 0 and atomtypes_end >= 0:
                # 写入输出文件
                with open(output_file, 'w') as f:
                    f.writelines(atomtypes_content)

                # 从原文件中删除 atomtypes 段
                new_content = lines[:atomtypes_start] + lines[atomtypes_end:]

                # 重写原文件（删除 atomtypes 段）
                with open(itp_file, 'w') as f:
                    f.writelines(new_content)

                self.logger.info(f"成功提取并剪切原子类型参数：{output_file}")
                self.logger.info(f"从 {itp_file} 中删除了 {atomtypes_end - atomtypes_start} 行")
                return True
            else:
                self.logger.warning(f"在 {itp_file} 中未找到 [ atomtypes ] 段")
                return False

        except Exception as e:
            self.logger.error(f"提取原子类型参数失败：{e}")
            return False

    def copy_template_files(self, sample_dir):
        """复制模板文件"""
        self.logger.info(f"复制模板文件到：{sample_dir}")

        try:
            template_dir = Path(self.file_dir)
            target_dir = Path(sample_dir)

            # 复制所有文件
            for file_path in template_dir.iterdir():
                if file_path.is_file():
                    shutil.copy2(file_path, target_dir / file_path.name)

            self.logger.info("成功复制模板文件")
            return True

        except Exception as e:
            self.logger.error(f"复制模板文件失败：{e}")
            return False

    def process_single_sample(self, sample_data, progress_bar):
        """处理单个样本"""
        sample_id = sample_data['Sample']
        cid_a = sample_data['CID_A']
        cid_b = sample_data['CID_B']

        self.logger.info(f"开始处理样本：{sample_id}")
        start_time = time.time()

        try:
            # 创建样本目录
            sample_dir = Path(str(sample_id))
            sample_dir.mkdir(exist_ok=True)

            # 步骤1：下载分子A结构
            moa_sdf = sample_dir / "MOA.sdf"
            if not self.download_pubchem_structure(cid_a, moa_sdf, sample_dir):
                raise RuntimeError(f"无法下载分子A结构 (CID: {cid_a})")

            # 步骤2：下载分子B结构
            mob_sdf = sample_dir / "MOB.sdf"
            if not self.download_pubchem_structure(cid_b, mob_sdf, sample_dir):
                raise RuntimeError(f"无法下载分子B结构 (CID: {cid_b})")

            # 步骤3：格式转换
            moa_mol2 = sample_dir / "MOA.mol2"
            mob_mol2 = sample_dir / "MOB.mol2"

            if not self.convert_sdf_to_mol2(moa_sdf, moa_mol2):
                raise RuntimeError("分子A格式转换失败")

            if not self.convert_sdf_to_mol2(mob_sdf, mob_mol2):
                raise RuntimeError("分子B格式转换失败")

            # 步骤4：分子几何优化
            if not self.optimize_molecule_geometry(moa_mol2):
                raise RuntimeError("分子A几何优化失败")

            if not self.optimize_molecule_geometry(mob_mol2):
                raise RuntimeError("分子B几何优化失败")

            # 步骤5：残基名称标准化
            if not self.run_replace_resname(str(sample_dir)):
                raise RuntimeError("残基名称标准化失败")

            # 步骤6：AMBER拓扑生成
            if not self.generate_amber_topology(moa_mol2, "MOA", 0):
                raise RuntimeError("分子A拓扑生成失败")

            if not self.generate_amber_topology(mob_mol2, "MOB", 0):
                raise RuntimeError("分子B拓扑生成失败")

            # 步骤7：文件整理
            moa_acpype_dir = sample_dir / "MOA.acpype"
            mob_acpype_dir = sample_dir / "MOB.acpype"

            if moa_acpype_dir.exists():
                shutil.copy2(moa_acpype_dir / "MOA_GMX.gro", sample_dir)
                shutil.copy2(moa_acpype_dir / "MOA_GMX.itp", sample_dir)
                shutil.rmtree(moa_acpype_dir)

            if mob_acpype_dir.exists():
                shutil.copy2(mob_acpype_dir / "MOB_GMX.gro", sample_dir)
                shutil.copy2(mob_acpype_dir / "MOB_GMX.itp", sample_dir)
                shutil.rmtree(mob_acpype_dir)

            # 步骤8：原子类型参数提取
            moa_itp = sample_dir / "MOA_GMX.itp"
            mob_itp = sample_dir / "MOB_GMX.itp"
            moa_prm = sample_dir / "MOA_GMX_prm.itp"
            mob_prm = sample_dir / "MOB_GMX_prm.itp"

            if moa_itp.exists():
                self.extract_atomtypes(moa_itp, moa_prm)

            if mob_itp.exists():
                self.extract_atomtypes(mob_itp, mob_prm)

            # 步骤9：模拟文件部署
            if not self.copy_template_files(sample_dir):
                raise RuntimeError("模板文件部署失败")

            # 记录完成时间
            end_time = time.time()
            processing_time = end_time - start_time

            self.logger.info(f"样本 {sample_id} 处理完成，耗时：{processing_time:.2f}秒")
            self.completed_samples.append(sample_id)
            self.stats['completed'] += 1

            return True

        except Exception as e:
            end_time = time.time()
            processing_time = end_time - start_time

            self.logger.error(f"样本 {sample_id} 处理失败：{e}")
            self.logger.error(f"处理耗时：{processing_time:.2f}秒")
            self.failed_samples.append({
                'Sample': sample_id,
                'Error': str(e),
                'Time': processing_time
            })
            self.stats['failed'] += 1

            return False

    def generate_report(self):
        """生成处理报告"""
        self.logger.info("生成处理报告...")

        report_file = self.log_dir / f"processing_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"

        try:
            with open(report_file, 'w', encoding='utf-8') as f:
                f.write("=== GROMACS模拟文件处理报告 ===\n")
                f.write(f"生成时间：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"\n")

                f.write("=== 统计信息 ===\n")
                f.write(f"总样本数：{self.stats['total']}\n")
                f.write(f"成功处理：{self.stats['completed']}\n")
                f.write(f"处理失败：{self.stats['failed']}\n")
                f.write(f"跳过样本：{self.stats['skipped']}\n")
                f.write(f"成功率：{(self.stats['completed'] / self.stats['total'] * 100):.1f}%\n")
                f.write(f"\n")

                if self.completed_samples:
                    f.write("=== 成功处理的样本 ===\n")
                    for sample in self.completed_samples:
                        f.write(f"[成功] 样本 {sample}\n")
                    f.write(f"\n")

                if self.failed_samples:
                    f.write("=== 处理失败的样本 ===\n")
                    for failed in self.failed_samples:
                        f.write(f"[失败] 样本 {failed['Sample']}: {failed['Error']}\n")
                    f.write(f"\n")

                f.write("=== 日志文件 ===\n")
                log_files = list(self.log_dir.glob("gromacs_processor_*.log"))
                for log_file in sorted(log_files):
                    f.write(f"{log_file}\n")

            self.logger.info(f"处理报告已生成：{report_file}")

        except Exception as e:
            self.logger.error(f"生成处理报告失败：{e}")

    def run(self):
        """运行处理流程"""
        try:
            # 读取和验证数据
            self.read_mol_csv()

            # 检查已完成的样本
            self.check_completed_samples()

            if self.stats['total'] == 0:
                self.logger.error("没有找到有效的样本数据")
                return False

            # 过滤待处理的样本
            samples_to_process = [
                sample for sample in self.samples
                if sample['Sample'] not in self.completed_samples
            ]

            if not samples_to_process:
                self.logger.info("所有样本已完成，无需处理")
                return True

            self.logger.info(f"开始处理 {len(samples_to_process)} 个样本")

            # 创建进度条
            with tqdm(
                total=len(samples_to_process),
                desc="处理进度",
                unit="sample"
            ) as progress_bar:

                # 处理每个样本
                for sample_data in samples_to_process:
                    success = self.process_single_sample(sample_data, progress_bar)
                    progress_bar.update(1)

                    if success:
                        progress_bar.set_postfix({"状态": "成功"})
                    else:
                        progress_bar.set_postfix({"状态": "失败"})

            # 生成报告
            self.generate_report()

            # 输出最终统计
            self.logger.info("=== 处理完成 ===")
            self.logger.info(f"总样本数：{self.stats['total']}")
            self.logger.info(f"成功处理：{self.stats['completed']}")
            self.logger.info(f"处理失败：{self.stats['failed']}")
            self.logger.info(f"跳过样本：{self.stats['skipped']}")
            self.logger.info(f"成功率：{(self.stats['completed'] / self.stats['total'] * 100):.1f}%")

            return self.stats['failed'] == 0

        except Exception as e:
            self.logger.error(f"处理流程出现错误：{e}")
            return False


def main():
    """主函数"""
    import argparse

    parser = argparse.ArgumentParser(description='GROMACS模拟文件自动生成脚本')
    parser.add_argument('--csv', default='Mol.csv', help='Mol.csv文件路径')
    parser.add_argument('--file-dir', default='File', help='模板文件目录')
    parser.add_argument('--log-level', default='INFO',
                       choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                       help='日志级别')

    args = parser.parse_args()

    try:
        processor = GromacsProcessor(
            mol_csv_path=args.csv,
            file_dir=args.file_dir,
            log_level=getattr(logging, args.log_level)
        )

        success = processor.run()

        if success:
            print("\n[成功] 所有样本处理成功完成！")
            sys.exit(0)
        else:
            print("\n[失败] 部分样本处理失败，请查看日志文件")
            sys.exit(1)

    except Exception as e:
        print(f"\n[错误] 运行失败：{e}")
        sys.exit(1)


if __name__ == "__main__":
    main()