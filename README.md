# GROMACS分子动力学模拟前处理自动化系统

## 项目概述

本项目是一个在WSL Ubuntu环境下运行的自动化Python脚本，用于批量处理分子对并生成GROMACS分子动力学模拟所需的完整文件集。该系统实现了从分子结构下载到拓扑生成的全自动化流程。

## 🛠️ 系统环境

**操作系统**：WSL2 Ubuntu
**必需软件**：
- Python 3.x (pandas, requests, subprocess, tqdm, logging)
- Open Babel
- acpype
- conda (已配置gcc环境)
- GROMACS (用于生产模拟)

## 📁 项目结构

```
MD_Pre/
├── README.md                     # 项目说明文档                   # Claude Code项目配置
├── Mol.csv                      # 样本清单文件
├── replace_resname.py           # 残基重命名脚本
├── main.py                      # 主处理脚本
├── File/                        # GROMACS模拟模板文件夹
│   ├── topol.top               # 系统拓扑文件
│   ├── run_all.sh              # 自动运行脚本
│   ├── em.mdp                  # 能量最小化参数
│   ├── nvt.mdp                 # NVT平衡参数
│   ├── npt.mdp                 # NPT平衡参数
│   ├── md.mdp                  # 生产模拟参数
│   ├── ions.mdp                # 离子化参数
│   ├── tip3p.gro               # TIP3P水模型
│   └── [其他配置文件]
└── [生成的样本目录]
    ├── 1/
    │   ├── MOA.mol2            # 优化后的分子A结构
    │   ├── MOB.mol2            # 优化后的分子B结构
    │   ├── MOA_GMX.gro         # GROMACS坐标文件
    │   ├── MOA_GMX.itp         # 完整拓扑文件
    │   ├── MOA_GMX_prm.itp     # 原子类型参数
    │   ├── MOB_GMX.gro
    │   ├── MOB_GMX.itp
    │   ├── MOB_GMX_prm.itp
    │   └── [File/中的所有文件] # 模拟配置文件
    └── ...
```

## 📋 输入文件格式

### Mol.csv文件

```csv
sample,CID_A,CID_B
1,3715,2179
2,3715,11556711
3,3715,445643
```

**列说明**：
- **sample**：样本编号（将作为输出目录名）
- **CID_A**：分子A的PubChem化合物ID (CID)
- **CID_B**：分子B的PubChem化合物ID (CID)

## 🚀 核心处理流程

### 1. 分子结构下载
从PubChem下载分子结构文件：
- 优先获取3D SDF格式：`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{CID}/SDF?record_type=3d`
- 若3D不可用，获取2D格式：`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{CID}/SDF`

文件保存：
- 分子A → `{Sample}/MOA.sdf`
- 分子B → `{Sample}/MOB.sdf`

### 2. 格式转换 (SDF → MOL2)
```bash
obabel MOA.sdf -O MOA.mol2
obabel MOB.sdf -O MOB.mol2
rm MOA.sdf MOB.sdf
```

### 3. 分子几何优化
使用MMFF94力场优化两个mol2文件：
- 力场：MMFF94
- 介电常数：78 (模拟水环境)
- 其他参数：默认值

### 4. 残基名称标准化
```bash
cd MD_Pre/
python replace_resname.py
```
- 自动将样本文件夹中MOA.mol2的残基名改为"MOA"
- 自动将样本文件夹中MOB.mol2的残基名改为"MOB"

### 5. 设置默认电荷
所有分子统一使用默认电荷值：0

### 6. AMBER拓扑生成
```bash
cd {Sample}/
conda activate gcc
acpype -i MOA.mol2 -c bcc -a gaff2 -n 0 -b MOA
acpype -i MOB.mol2 -c bcc -a gaff2 -n 0 -b MOB
```

**参数说明**：
- `-c bcc`：AM1-BCC电荷计算方法
- `-a gaff2`：GAFF2力场
- `-n 0`：分子净电荷默认为0
- `-b {basename}`：输出文件基础名

### 7. 文件整理
复制关键文件到样本目录并清理临时文件夹：
```bash
cp MOA.acpype/MOA_GMX.gro ./
cp MOA.acpype/MOA_GMX.itp ./
rm -rf MOA.acpype/
```

### 8. 原子类型参数提取
从`.itp`文件中提取`[ atomtypes ]`部分到独立文件：
- 输出：`MOA_GMX_prm.itp`和`MOB_GMX_prm.itp`

### 9. 模拟文件部署
```bash
cp -r ../File/* ./
```
将预配置的GROMACS模拟文件复制到样本目录。

## 📊 产出结果

每个样本目录将包含完整的GROMACS模拟文件：

### 分子文件
- `MOA.mol2` - 优化后的分子A结构
- `MOB.mol2` - 优化后的分子B结构
- `MOA_GMX.gro` - GROMACS坐标文件
- `MOA_GMX.itp` - 完整拓扑文件
- `MOA_GMX_prm.itp` - 原子类型参数
- `MOB_GMX.gro` - GROMACS坐标文件
- `MOB_GMX.itp` - 完整拓扑文件
- `MOB_GMX_prm.itp` - 原子类型参数

### 模拟配置文件
- `topol.top` - 系统拓扑文件
- `run_all.sh` - 自动运行脚本
- `em.mdp` - 能量最小化参数
- `nvt.mdp` - NVT平衡参数
- `npt.mdp` - NPT平衡参数
- `md.mdp` - 生产模拟参数
- `ions.mdp` - 离子化参数
- `tip3p.gro` - TIP3P水模型

## 🔧 安装与配置

### 1. 环境依赖安装

```bash
# 更新系统
sudo apt update && sudo apt upgrade

# 安装Python依赖
pip install pandas requests tqdm

# 安装Open Babel
sudo apt install openbabel

# 安装conda（如未安装）
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# 配置conda环境
conda create -n gcc gcc gfortran openbabel acpype -c conda-forge -c bioconda
conda activate gcc
```

### 2. GROMACS安装

```bash
# 安装GROMACS
sudo apt install gromacs

# 验证安装
gmx --version
```

## 🚀 使用方法

### 基本运行

```bash
# 进入项目目录
cd MD_Pre

# 运行主脚本
python main.py
```
### 运行GROMACS模拟

```bash
# 进入样本目录
cd 1/

# 运行模拟脚本
bash run_all.sh
```

该脚本将自动执行以下步骤：
1. 创建模拟盒子并插入分子
2. 溶剂化
3. 添加离子
4. 能量最小化
5. NVT平衡
6. NPT平衡
7. 生产模拟

## 📝 日志系统

系统自动生成详细的日志文件：

- **主日志**：`MD_Pre_processing_YYYYMMDD_HHMMSS.log`
- **样本日志**：每个样本目录中的`processing_log.txt`
- **错误日志**：`error_samples.txt`（记录失败样本）

日志内容包括：
- 处理进度和状态
- 错误详情和调试信息
- 处理时间统计
- 成功率统计

## 🛡️ 健壮性特性

### 错误处理
- **下载失败**：记录错误，跳过该样本，继续处理下一个
- **格式转换失败**：记录详细错误信息，提供调试建议
- **拓扑生成失败**：记录acpype命令输出，分析失败原因
- **文件操作失败**：提供详细的错误信息和恢复建议

### 进度监控
- 实时显示处理进度条
- 显示当前处理的样本信息
- 提供剩余时间预估

## 📈 输出统计

处理完成后，系统会生成统计报告，包括：
- 总样本数
- 成功样本数
- 失败样本数
- 跳过样本数
- 平均处理时间
- 成功率

## 🔍 故障排除

### 常见问题

1. **PubChem下载失败**
   - 检查网络连接
   - 验证CID的有效性
   - 查看日志中的错误详情

2. **acpype失败**
   - 确保conda环境已激活
   - 检查GAFF2力场是否正确安装
   - 验证mol2文件格式

3. **GROMACS模拟失败**
   - 检查GROMACS版本兼容性
   - 验证力场参数
   - 查看详细的错误日志


## 📊 性能优化

### 批处理建议
- 建议100个样本为一批进行处理
- 定期备份处理进度
- 监控系统资源使用情况

### 存储管理
- 删除临时文件以节省空间
- 定期归档已完成的项目
- 使用压缩存储大型输出文件

---

## 🚀 快速开始

1. **准备环境**
   ```bash
   # 确保所有依赖已安装
   conda activate gcc
   gmx --version
   ```

2. **配置输入文件**
   ```bash
   # 编辑Mol.csv文件，添加你的样本数据
   nano Mol.csv
   ```

3. **运行脚本**
   ```bash
   # 启动处理流程
   python main.py
   ```

4. **查看结果**
   ```bash
   # 检查生成的样本目录
   ls -la
   ```

5. **运行模拟**
   ```bash
   # 进入样本目录并运行模拟
   cd 1/
   bash run_all.sh
   ```

**祝您使用愉快！** 🎉