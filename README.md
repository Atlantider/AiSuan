# 计算材料科学平台 (AiSuan)

## 项目简介
计算材料科学平台是一个面向材料科学研究人员的在线计算服务平台，提供电池材料、催化材料等多种计算服务。该平台旨在简化复杂的材料计算流程，使研究人员能够轻松地进行各种材料性能预测和设计。

## 主要功能
- 多种材料计算模块（电池材料、催化材料等）
- 交互式计算工作台
- 计算模板自动生成
- 计算结果可视化
- 材料数据库集成

### 电池材料计算
- 电极电位计算
- 离子迁移能垒计算
- 形成能计算
- 电子结构计算
- 相稳定性分析
- 体积变化计算

### 催化材料计算
- 表面催化计算
- 电催化计算
- 光催化计算
- 催化剂筛选与设计

## 技术栈
### 前端
- React.js (v18+)
- Ant Design (v5+)
- Redux Toolkit
- React Router (v6+)
- ECharts (v5+)
- Three.js 与 React Three Fiber
- 3Dmol.js
- Axios
- Vite

### 后端
- Django 4.2+ 与 Django REST framework
- Celery
- Redis
- Swagger/OpenAPI
- JWT
- Django Admin

### 计算集成
- pymatgen
- ASE
- GPAW
- atomate2
- MDAnalysis
- NumPy, SciPy, pandas
- scikit-learn, PyTorch

### 部署
- Docker, Docker Compose
- Kubernetes
- GitHub Actions
- Nginx + Gunicorn
- Prometheus + Grafana

## 快速开始

### 开发环境设置
1. 克隆仓库
```bash
git clone https://github.com/Atlantider/AiSuan.git
cd AiSuan
```

2. 安装依赖
```bash
cd materialsci-platform
npm install
```

3. 启动开发服务器
```bash
npm run dev
```

### 生产环境部署
1. 构建前端
```bash
npm run build
```

2. 使用 Docker Compose 部署
```bash
docker-compose up -d
```

## 贡献指南
1. Fork 本仓库
2. 创建您的特性分支 (git checkout -b feature/AmazingFeature)
3. 提交您的更改 (git commit -m 'Add some AmazingFeature')
4. 推送到分支 (git push origin feature/AmazingFeature)
5. 创建一个 Pull Request

## 许可证
本项目采用 MIT 许可证 - 查看 [LICENSE](LICENSE) 文件了解详细信息。

## 联系我们
- 项目维护者：[维护者姓名]
- 邮箱：[联系邮箱]
- 项目链接：[https://github.com/Atlantider/AiSuan](https://github.com/Atlantider/AiSuan) 