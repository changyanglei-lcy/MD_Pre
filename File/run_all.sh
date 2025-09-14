#!/bin/bash
set -e  # 出错即退出
set -o pipefail  # 防止管道中的错误被忽略

# Step 1: 创建模拟盒子并插入分子
echo ">>> Step 1: 创建模拟盒子并插入分子"
gmx editconf -f 	MOA_GMX.gro -o box.gro -box 6 6 6 -bt cubic -d 1.2 -c
gmx insert-molecules -ci MOA_GMX.gro -nmol 9 -f box.gro -o system1.gro -try 200 -rot xyz -seed 123
gmx insert-molecules -ci MOB_GMX.gro -nmol 10 -f system1.gro -o system.gro -try 500 -rot xyz -seed 456

# Step 2: 溶剂化
echo ">>> Step 2: 溶剂化"
gmx solvate -cp system.gro -cs tip3p.gro -o solvated.gro -p topol.top

# Step 3: 添加离子
echo ">>> Step 3: 添加离子"
gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 100
echo "SOL" | gmx genion -s ions.tpr -o ions.gro -p topol.top -pname NA -nname CL -neutral

# Step 4: 能量最小化（steep）
echo ">>> Step 4: 能量最小化（steep）"
gmx grompp -f em.mdp -c ions.gro -p topol.top -o em.tpr -maxwarn 100
gmx mdrun -v -deffnm em

# Step 6: 创建索引文件（非交互式）
echo ">>> Step 6: 创建索引文件"
echo -e "2 | 3\nname 7 complex\nq" | gmx make_ndx -f em.gro -o index.ndx

# Step 7: NVT 平衡
echo ">>> Step 7: NVT 平衡"
gmx grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr -n index.ndx -maxwarn 100
gmx mdrun -v -deffnm nvt

# Step 8: NPT 平衡
echo ">>> Step 8: NPT 平衡"
gmx grompp -f npt.mdp -c nvt.gro -p topol.top -o npt.tpr -n index.ndx -maxwarn 100
gmx mdrun -v -deffnm npt

# Step 9: 生产模拟
echo ">>> Step 9: 生产模拟"
gmx grompp -f md.mdp -c npt.gro -p topol.top -o md.tpr -n index.ndx -maxwarn 100
gmx mdrun -v -deffnm md

echo ">>> 所有步骤完成 ✅"