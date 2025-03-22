#!/bin/bash

# 加载环境变量设置脚本
echo "加载环境变量..."
source /public/home/xiaoji/AiSuan/setup_env.sh

# 添加额外的环境变量
export PATH=/public/home/xiaoji/software/redis/src:/public/home/xiaoji/software/bin:$PATH
export LD_LIBRARY_PATH=~/lib:$LD_LIBRARY_PATH

# 输出环境信息
echo "PATH: $PATH"
echo "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"

# 启动Redis服务
echo "启动Redis服务..."
redis-server --daemonize yes

# 检查Redis是否启动成功
echo "检查Redis状态..."
redis-cli ping

# 启动Django后端服务
echo "启动Django后端服务..."
cd /public/home/xiaoji/AiSuan/backend
source /public/software/anaconda3/bin/activate molyte
python manage.py runserver 0.0.0.0:8000 --noreload --nothreading > /tmp/django.log 2>&1 &
cd /public/home/xiaoji/AiSuan

# 启动Celery服务（包含Gaussian和Multiwfn环境变量）
echo "启动Celery服务..."
cd /public/home/xiaoji/AiSuan/backend
source /public/software/anaconda3/bin/activate molyte

# 为Celery环境设置Gaussian和Multiwfn
export g16root=/public/software/g16
export GAUSS_EXEDIR=$g16root
export PATH=$g16root:$PATH
export Multiwfnpath=/public/home/xiaoji/software/Multiwfn_3.7_bin_Linux
export PATH=$Multiwfnpath:$PATH
ulimit -s unlimited
export OMP_STACKSIZE=200M

# 验证计算环境
echo "Celery计算环境:"
echo -n "g16: "
which g16 > /dev/null 2>&1 && echo "可用" || echo "不可用"
echo -n "Multiwfn: "
which Multiwfn > /dev/null 2>&1 && echo "可用" || echo "不可用"

# 启动Celery worker
PYTHONPATH=/public/home/xiaoji/AiSuan celery -A electrolyte_backend worker --loglevel=INFO --concurrency=1 > /tmp/celery.log 2>&1 &
cd /public/home/xiaoji/AiSuan

# 启动前端服务
echo "启动前端服务..."
cd /public/home/xiaoji/AiSuan/materialsci-platform
npm run dev > /tmp/frontend.log 2>&1 &
cd /public/home/xiaoji/AiSuan

echo "所有服务已启动"
echo "Django日志: /tmp/django.log"
echo "Celery日志: /tmp/celery.log"
echo "前端日志: /tmp/frontend.log" 