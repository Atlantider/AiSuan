#!/bin/bash
# AiSuan环境变量设置脚本

# Gaussian 16设置
export g16root=/public/software/g16
export GAUSS_EXEDIR=$g16root
export PATH=$g16root:$PATH

# Multiwfn设置
ulimit -s unlimited
export OMP_STACKSIZE=200M
export Multiwfnpath=/public/home/xiaoji/software/Multiwfn_3.7_bin_Linux
# 直接指定完整路径
export PATH=$PATH:$Multiwfnpath
# 建立符号链接到用户bin目录
mkdir -p $HOME/bin
ln -sf $Multiwfnpath/Multiwfn $HOME/bin/Multiwfn
export PATH=$HOME/bin:$PATH

# 验证设置
echo "环境变量已设置:"
echo "===== Gaussian 16 ====="
echo "g16root = $g16root"
echo "GAUSS_EXEDIR = $GAUSS_EXEDIR"
echo "===== Multiwfn ====="
echo "Multiwfnpath = $Multiwfnpath"
echo "OMP_STACKSIZE = $OMP_STACKSIZE"
echo "符号链接: $HOME/bin/Multiwfn -> $Multiwfnpath/Multiwfn"
echo "===== PATH ====="
echo "PATH = $PATH"

# 检查命令是否可用
echo "===== 命令验证 ====="
echo -n "g16: "
which g16 > /dev/null 2>&1 && echo "可用" || echo "不可用"
echo -n "Multiwfn: "
which Multiwfn > /dev/null 2>&1 && echo "可用" || echo "不可用" 