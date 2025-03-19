# 电解液计算模块开发计划

## 标准前后端开发流程

一个规范的前后端分离项目开发流程通常包括以下阶段：

### 1. 需求分析与规划阶段 ✅
- 确定功能需求和业务逻辑 ✅
- 制定技术方案和架构设计 ✅
- 确定API接口规范 ✅
- 创建项目时间线和里程碑 ✅

### 2. 设计阶段 ✅
- 数据库模型设计 ✅
- API接口设计 ✅
- 前端界面和交互设计 ✅
- 技术选型确认 ✅

### 3. 开发阶段 (进行中)
- **后端开发** (计划中):
  - 搭建项目框架
  - 实现数据模型
  - 开发API接口
  - 集成第三方服务
  
- **前端开发** (部分完成):
  - 搭建UI框架 ✅
  - 实现页面组件 ✅ (电解液计算页面已完成)
  - 集成API调用 (待完成)
  - 状态管理实现 (待完成)

### 4. 测试阶段
- 单元测试
- 集成测试
- API测试
- 用户界面测试
- 压力测试（如需要）

### 5. 部署阶段
- 准备生产环境
- 数据库迁移
- 配置服务器
- 实施CI/CD流程

### 6. 维护与迭代阶段
- 监控系统运行
- 处理用户反馈
- 进行功能迭代
- 性能优化

## 电解液计算模块具体开发计划

基于backend.md文档和标准开发流程，制定以下开发计划：

## 已完成的任务

1. **需求分析与规划** ✅
   - 分析了电解液计算模块的功能需求 ✅
   - 确定了Molyte作为后端计算引擎 ✅
   - 编写了详细的backend.md技术方案文档 ✅

2. **前端开发** (部分完成)
   - 实现了电解液计算页面`ElectrolyteCalculation.tsx` ✅
     - 阴阳离子选择和比例设置 ✅
     - 溶剂选择（包括自定义溶剂）和比例设置 ✅
     - 温度和摩尔浓度设置 ✅
     - 计算类型选择 ✅
   - 使用React + TypeScript + Ant Design实现前端界面 ✅
   - 优化了页面布局和交互体验 ✅

## 待完成的任务

### 第1周：项目初始化与基础架构
| 天数 | 任务 | 说明 | 状态 |
|-----|------|------|------|
| 1-2天 | 创建Django项目结构 | 创建项目、应用并配置基本设置 | 待完成 |
| 3天 | 设计并实现数据模型 | 实现Job、ElectrolyteFormula、SimulationResult等模型 | 待完成 |
| 1-2天 | 配置Celery和Redis | 设置任务队列和消息代理 | 待完成 |
| 1-2天 | 设置认证系统 | 实现JWT认证和权限控制 | 待完成 |

### 第2-3周：核心功能开发
| 天数 | 任务 | 说明 | 状态 |
|-----|------|------|------|
| 2-3天 | 集成Molyte库 | 将Molyte作为Python包集成到项目中 | 待完成 |
| 3-4天 | 开发LAMMPS文件生成器 | 实现从电解液配方生成LAMMPS输入文件的功能 | 待完成 |
| 2-3天 | 实现命令执行器 | 开发调用LAMMPS执行计算的工具类 | 待完成 |
| 2-3天 | 实现Excel处理功能 | 开发上传和下载Excel模板的功能 | 待完成 |
| 3-4天 | 实现Celery任务 | 开发异步计算任务和状态管理 | 待完成 |

### 第4-5周：API和分析功能开发
| 天数 | 任务 | 说明 | 状态 |
|-----|------|------|------|
| 3-4天 | 实现JobViewSet | 开发作业管理相关API | 待完成 |
| 3-4天 | 实现ElectrolyteViewSet | 开发电解液计算相关API | 待完成 |
| 3-4天 | 集成分析模块 | 实现结果分析功能和可视化生成 | 待完成 |
| 2-3天 | 开发文件下载API | 实现计算结果和图表下载功能 | 待完成 |
| 1-2天 | 编写API文档 | 使用Swagger/OpenAPI生成API文档 | 待完成 |

### 第6周：测试与前后端集成
| 天数 | 任务 | 说明 | 状态 |
|-----|------|------|------|
| 2-3天 | 编写单元测试 | 为模型、视图和任务编写测试 | 待完成 |
| 2天 | 进行API测试 | 测试所有API端点的功能 | 待完成 |
| 2-3天 | 与前端集成 | 协助前端集成API调用和数据处理 | 待完成 |
| 1-2天 | 修复问题和调整 | 解决测试中发现的问题 | 待完成 |

### 第7-8周：优化与部署
| 天数 | 任务 | 说明 | 状态 |
|-----|------|------|------|
| 2-3天 | 性能优化 | 优化计算任务和数据处理 | 待完成 |
| 2-3天 | 安全加固 | 加强输入验证和权限控制 | 待完成 |
| 2-3天 | 准备部署配置 | 编写Docker配置和部署脚本 | 待完成 |
| 2-3天 | 部署到测试环境 | 在测试服务器上部署和验证 | 待完成 |
| 1-2天 | 编写部署文档 | 准备部署和维护文档 | 待完成 |

## 具体实施步骤

### 后端开发起步

1. **创建Django项目和应用**:
   ```bash
   # 创建项目
   django-admin startproject aisuan_backend
   cd aisuan_backend
   
   # 创建应用
   python manage.py startapp jobs
   python manage.py startapp electrolytes
   ```

2. **配置settings.py**:
   ```python
   # 添加应用
   INSTALLED_APPS = [
       # Django内置应用
       'django.contrib.admin',
       'django.contrib.auth',
       'django.contrib.contenttypes',
       'django.contrib.sessions',
       'django.contrib.messages',
       'django.contrib.staticfiles',
       
       # 第三方应用
       'rest_framework',
       'rest_framework_simplejwt',
       'drf_yasg',
       'corsheaders',
       
       # 自定义应用
       'jobs',
       'electrolytes',
   ]
   
   # 配置数据库
   DATABASES = {
       'default': {
           'ENGINE': 'django.db.backends.postgresql',
           'NAME': 'aisuan_db',
           'USER': 'aisuan_user',
           'PASSWORD': 'your_password',
           'HOST': 'localhost',
           'PORT': '5432',
       }
   }
   
   # REST框架设置
   REST_FRAMEWORK = {
       'DEFAULT_AUTHENTICATION_CLASSES': (
           'rest_framework_simplejwt.authentication.JWTAuthentication',
       ),
       'DEFAULT_PERMISSION_CLASSES': [
           'rest_framework.permissions.IsAuthenticated',
       ],
   }
   
   # JWT设置
   SIMPLE_JWT = {
       'ACCESS_TOKEN_LIFETIME': timedelta(days=1),
       'REFRESH_TOKEN_LIFETIME': timedelta(days=7),
   }
   
   # CORS设置
   MIDDLEWARE = [
       'corsheaders.middleware.CorsMiddleware',
       # 其他中间件...
   ]
   
   CORS_ALLOWED_ORIGINS = [
       "http://localhost:3000",
       "http://localhost:5173",
   ]
   ```

3. **实现数据模型**:
   ```python
   # jobs/models.py
   from django.db import models
   from django.contrib.auth.models import User
   import uuid
   
   class Job(models.Model):
       """计算作业模型"""
       STATUS_CHOICES = (
           ('queued', '排队中'),
           ('running', '运行中'),
           ('completed', '已完成'),
           ('failed', '失败'),
       )
       
       id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
       user = models.ForeignKey(User, on_delete=models.CASCADE, related_name='jobs')
       name = models.CharField(max_length=255)
       status = models.CharField(max_length=20, choices=STATUS_CHOICES, default='queued')
       created_at = models.DateTimeField(auto_now_add=True)
       updated_at = models.DateTimeField(auto_now=True)
       
       def __str__(self):
           return f"{self.name} ({self.status})"
   
   # 运行迁移
   # python manage.py makemigrations
   # python manage.py migrate
   ```

4. **配置Celery**:
   ```python
   # aisuan_backend/celery.py
   from __future__ import absolute_import, unicode_literals
   import os
   from celery import Celery
   
   os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'aisuan_backend.settings')
   
   app = Celery('aisuan_backend')
   app.config_from_object('django.conf:settings', namespace='CELERY')
   app.autodiscover_tasks()
   
   # settings.py添加Celery配置
   CELERY_BROKER_URL = 'redis://localhost:6379/0'
   CELERY_RESULT_BACKEND = 'redis://localhost:6379/0'
   CELERY_ACCEPT_CONTENT = ['json']
   CELERY_TASK_SERIALIZER = 'json'
   CELERY_RESULT_SERIALIZER = 'json'
   CELERY_TIMEZONE = 'UTC'
   ```

5. **创建Serializers**:
   ```python
   # jobs/serializers.py
   from rest_framework import serializers
   from .models import Job, ElectrolyteFormula, SimulationResult
   
   class JobSerializer(serializers.ModelSerializer):
       class Meta:
           model = Job
           fields = ['id', 'name', 'status', 'created_at', 'updated_at']
           read_only_fields = ['id', 'status', 'created_at', 'updated_at']
   ```

6. **实现API视图**:
   ```python
   # jobs/views.py
   from rest_framework import viewsets, status
   from rest_framework.decorators import action
   from rest_framework.response import Response
   from .models import Job
   from .serializers import JobSerializer
   
   class JobViewSet(viewsets.ModelViewSet):
       """作业管理视图集"""
       queryset = Job.objects.all()
       serializer_class = JobSerializer
       
       def perform_create(self, serializer):
           serializer.save(user=self.request.user)
   ```

7. **配置URL路由**:
   ```python
   # aisuan_backend/urls.py
   from django.contrib import admin
   from django.urls import path, include
   from rest_framework.routers import DefaultRouter
   from jobs.views import JobViewSet
   
   router = DefaultRouter()
   router.register(r'jobs', JobViewSet)
   
   urlpatterns = [
       path('admin/', admin.site.urls),
       path('api/', include(router.urls)),
   ]
   ```

### 前端开发起步

1. **创建React前端项目**: ✅
   ```bash
   # 使用Vite创建
   npm create vite@latest materialsci-platform -- --template react-ts
   cd materialsci-platform
   npm install
   ```

2. **安装必要依赖**: ✅
   ```bash
   npm install antd @ant-design/icons axios react-router-dom redux @reduxjs/toolkit
   ```

3. **实现电解液计算页面组件**: ✅
   - `ElectrolyteCalculation.tsx`页面已实现，包含表单、选择器和步骤导航

4. **配置API服务**:
   ```typescript
   // src/services/api.ts
   import axios from 'axios';
   
   const API_URL = 'http://localhost:8000/api';
   
   const api = axios.create({
     baseURL: API_URL,
     headers: {
       'Content-Type': 'application/json',
     },
   });
   
   // 请求拦截器添加认证token
   api.interceptors.request.use(
     (config) => {
       const token = localStorage.getItem('token');
       if (token) {
         config.headers.Authorization = `Bearer ${token}`;
       }
       return config;
     },
     (error) => {
       return Promise.reject(error);
     }
   );
   
   export default api;
   ```

5. **实现Redux状态管理**:
   ```typescript
   // src/store/slices/authSlice.ts
   import { createSlice, PayloadAction } from '@reduxjs/toolkit';
   
   interface AuthState {
     token: string | null;
     isAuthenticated: boolean;
     user: any | null;
   }
   
   const initialState: AuthState = {
     token: localStorage.getItem('token'),
     isAuthenticated: !!localStorage.getItem('token'),
     user: null,
   };
   
   const authSlice = createSlice({
     name: 'auth',
     initialState,
     reducers: {
       setCredentials: (state, action: PayloadAction<{ token: string; user: any }>) => {
         const { token, user } = action.payload;
         state.token = token;
         state.isAuthenticated = true;
         state.user = user;
         localStorage.setItem('token', token);
       },
       logout: (state) => {
         state.token = null;
         state.isAuthenticated = false;
         state.user = null;
         localStorage.removeItem('token');
       },
     },
   });
   
   export const { setCredentials, logout } = authSlice.actions;
   export default authSlice.reducer;
   ```

6. **配置路由**:
   ```typescript
   // src/App.tsx
   import { BrowserRouter as Router, Routes, Route } from 'react-router-dom';
   import Login from './pages/Login';
   import Dashboard from './pages/Dashboard';
   import ElectrolyteCalculation from './pages/battery/ElectrolyteCalculation';
   
   function App() {
     return (
       <Router>
         <Routes>
           <Route path="/login" element={<Login />} />
           <Route path="/" element={<Dashboard />} />
           <Route path="/battery/electrolyte" element={<ElectrolyteCalculation />} />
         </Routes>
       </Router>
     );
   }
   
   export default App;
   ```

## 后续开发步骤

随着项目进展，将按照上述计划继续开发以下功能：

1. **完善后端模型和API**:
   - 实现ElectrolyteFormula模型和API
   - 实现SimulationResult模型和API
   - 开发文件上传和下载API

2. **实现Celery任务**:
   - 开发电解液模拟计算任务
   - 实现任务状态跟踪和通知

3. **集成Molyte库**:
   - 将Molyte作为Python包集成
   - 实现文件生成和结果分析功能

4. **开发前端功能**:
   - 完善电解液计算页面
   - 实现结果展示和可视化
   - 开发作业管理和监控界面

5. **测试与部署**:
   - 编写单元测试和集成测试
   - 准备Docker配置和部署脚本
   - 部署到测试和生产环境 