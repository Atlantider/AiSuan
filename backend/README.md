# 电解质计算后端

这是一个基于Django和Celery的电解质计算后端系统，用于模拟和计算电解质性质。

## 功能特点

- 溶剂和盐的管理
- 电解质配方创建和管理
- 基于Celery的异步计算任务
- RESTful API接口
- JWT认证
- Swagger API文档

## 技术栈

- Django 5.1
- Django REST Framework
- Celery
- Redis
- PostgreSQL (可选)

## 安装和设置

### 前提条件

- Python 3.10+
- Redis服务器

### 安装步骤

1. 克隆仓库

```bash
git clone https://github.com/yourusername/electrolyte-backend.git
cd electrolyte-backend
```

2. 创建虚拟环境

```bash
python -m venv molyte_env
source molyte_env/bin/activate  # Linux/Mac
# 或
molyte_env\Scripts\activate  # Windows
```

3. 安装依赖

```bash
pip install -r requirements.txt
```

4. 配置环境变量（可选）

创建`.env`文件并设置以下变量：

```
DEBUG=True
SECRET_KEY=your-secret-key
DATABASE_URL=sqlite:///db.sqlite3
REDIS_URL=redis://localhost:6379/0
```

5. 运行数据库迁移

```bash
python manage.py migrate
```

6. 加载初始数据

```bash
python manage.py load_initial_data
```

7. 创建超级用户

```bash
python manage.py createsuperuser
```

## 运行服务

1. 启动Django开发服务器

```bash
python manage.py runserver
```

2. 启动Celery工作进程

```bash
celery -A electrolyte_backend worker -l info
```

## API文档

API文档可通过以下URL访问：

- Swagger UI: http://localhost:8000/swagger/
- ReDoc: http://localhost:8000/redoc/

## API端点

### 认证

- `POST /api/token/` - 获取JWT令牌
- `POST /api/token/refresh/` - 刷新JWT令牌
- `POST /api/token/verify/` - 验证JWT令牌

### 溶剂

- `GET /api/v1/solvents/` - 获取所有溶剂
- `POST /api/v1/solvents/` - 创建新溶剂
- `GET /api/v1/solvents/{id}/` - 获取特定溶剂
- `PUT /api/v1/solvents/{id}/` - 更新溶剂
- `DELETE /api/v1/solvents/{id}/` - 删除溶剂

### 盐

- `GET /api/v1/salts/` - 获取所有盐
- `POST /api/v1/salts/` - 创建新盐
- `GET /api/v1/salts/{id}/` - 获取特定盐
- `PUT /api/v1/salts/{id}/` - 更新盐
- `DELETE /api/v1/salts/{id}/` - 删除盐

### 电解质配方

- `GET /api/v1/formulations/` - 获取用户的所有配方
- `POST /api/v1/formulations/` - 创建新配方
- `GET /api/v1/formulations/{id}/` - 获取特定配方
- `PUT /api/v1/formulations/{id}/` - 更新配方
- `DELETE /api/v1/formulations/{id}/` - 删除配方
- `GET /api/v1/formulations/{id}/components/` - 获取配方组件
- `GET /api/v1/formulations/{id}/parameters/` - 获取配方参数

### 计算

- `GET /api/v1/calculations/` - 获取用户的所有计算
- `POST /api/v1/calculations/` - 创建新计算
- `GET /api/v1/calculations/{id}/` - 获取特定计算
- `DELETE /api/v1/calculations/{id}/` - 删除计算
- `GET /api/v1/calculations/{id}/result/` - 获取计算结果
- `POST /api/v1/calculations/{id}/restart/` - 重启计算

## 示例请求

### 创建电解质配方

```json
POST /api/v1/formulations/
{
  "name": "EC-DMC LiPF6",
  "description": "标准锂离子电池电解液",
  "salts": [
    {
      "id": 1,
      "concentration": 1.0
    }
  ],
  "solvents": [
    {
      "id": 1,
      "concentration": 0.5
    },
    {
      "id": 2,
      "concentration": 0.5
    }
  ],
  "temperature": 298.15,
  "pressure": 1.0,
  "time_step": 1.0,
  "equilibration_steps": 10000,
  "production_steps": 50000,
  "cutoff": 12.0
}
```

### 创建计算任务

```json
POST /api/v1/calculations/
{
  "name": "EC-DMC LiPF6计算",
  "description": "标准电解液离子电导率计算",
  "formulation_id": 1
}
```

## 许可证

MIT 