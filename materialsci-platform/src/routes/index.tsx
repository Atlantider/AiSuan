import React, { lazy, Suspense } from 'react';
import { Routes, Route, Navigate } from 'react-router-dom';
import { Spin } from 'antd';
import AuthGuard from '../components/common/AuthGuard';
import DashboardLayout from '../layouts/DashboardLayout';

// 懒加载页面组件
const Home = lazy(() => import('../pages/Home'));
const Login = lazy(() => import('../pages/Login'));
const Register = lazy(() => import('../pages/Register'));
const NotFound = lazy(() => import('../pages/NotFound'));
const Documentation = lazy(() => import('../pages/Documentation'));
const About = lazy(() => import('../pages/About'));

// 计算模块选择页面和计算详情页面
const CalculationsIndex = lazy(() => import('../pages/calculations/Index'));
const CalculationDetail = lazy(() => import('../pages/calculations/Detail'));

// 电池计算模块
const BatteryIndex = lazy(() => import('../pages/battery/Index'));
const ElectrodeMaterial = lazy(() => import('../pages/battery/ElectrodeMaterial'));
const ElectrolyteCalculation = lazy(() => import('../pages/battery/ElectrolyteCalculation'));
const FullBatterySystem = lazy(() => import('../pages/battery/FullBatterySystem'));
const BatteryWorkbench = lazy(() => import('../pages/battery/Workbench'));

// 催化计算模块
const CatalysisIndex = lazy(() => import('../pages/catalysis/Index'));

// 半导体材料计算模块
const SemiconductorIndex = lazy(() => import('../pages/semiconductor/Index'));

// 能源材料计算模块
const EnergyIndex = lazy(() => import('../pages/energy/Index'));

// 结构材料计算模块
const StructuralIndex = lazy(() => import('../pages/structural/Index'));

// 用户中心模块
const Dashboard = lazy(() => import('../pages/user/Dashboard'));
const TaskManagement = lazy(() => import('../pages/user/TaskManagement'));

// 加载组件
const LoadingComponent = () => (
  <div style={{ display: 'flex', justifyContent: 'center', alignItems: 'center', height: '100vh' }}>
    <Spin size="large" tip="页面加载中..." />
  </div>
);

// 页面不存在组件
const PageNotFound = () => (
  <div style={{ 
    display: 'flex', 
    flexDirection: 'column',
    justifyContent: 'center', 
    alignItems: 'center', 
    height: '70vh',
    textAlign: 'center'
  }}>
    <h1 style={{ fontSize: '72px', margin: '0', color: '#1C64F2' }}>404</h1>
    <h2>页面不存在</h2>
    <p>您访问的页面不存在或已被移除。</p>
    <button 
      onClick={() => window.location.href = '/'}
      style={{
        background: '#1C64F2',
        color: 'white',
        border: 'none',
        padding: '10px 20px',
        borderRadius: '4px',
        cursor: 'pointer',
        fontSize: '16px',
        marginTop: '20px'
      }}
    >
      返回首页
    </button>
  </div>
);

const AppRoutes = () => {
  return (
    <Suspense fallback={<LoadingComponent />}>
      <Routes>
        {/* 公开路由 */}
        <Route path="/" element={<Home />} />
        <Route path="/login" element={<Login />} />
        <Route path="/register" element={<Register />} />
        <Route path="/documentation" element={<Documentation />} />
        <Route path="/about" element={<About />} />
        
        {/* 计算模块选择页面 */}
        <Route path="/calculations" element={<CalculationsIndex />} />
        <Route path="/calculations/:id" element={<CalculationDetail />} />
        
        {/* 电池计算模块 */}
        <Route path="/battery" element={<BatteryIndex />} />
        <Route path="/battery/electrode-material" element={<ElectrodeMaterial />} />
        <Route path="/battery/electrode-material/workbench" element={<BatteryWorkbench calculationType="electrode" />} />
        <Route path="/battery/electrolyte" element={<ElectrolyteCalculation />} />
        <Route path="/battery/electrolyte/workbench" element={<BatteryWorkbench calculationType="electrolyte" />} />
        <Route path="/battery/full-battery" element={<FullBatterySystem />} />
        <Route path="/battery/full-battery/workbench" element={<BatteryWorkbench calculationType="fullBattery" />} />
        
        {/* 催化计算模块 */}
        <Route path="/catalysis" element={<CatalysisIndex />} />
        
        {/* 半导体材料计算模块 */}
        <Route path="/semiconductor" element={<SemiconductorIndex />} />
        
        {/* 能源材料计算模块 */}
        <Route path="/energy" element={<EnergyIndex />} />
        
        {/* 结构材料计算模块 */}
        <Route path="/structural" element={<StructuralIndex />} />
        
        {/* 用户中心路由 */}
        <Route element={<AuthGuard><DashboardLayout /></AuthGuard>}>
          <Route path="/user/dashboard" element={<Dashboard />} />
          <Route path="/user/tasks" element={<TaskManagement />} />
          <Route path="/user/data" element={<div>数据管理页面</div>} />
          <Route path="/user/workflows" element={<div>工作流管理页面</div>} />
          <Route path="/user/settings" element={<div>账户设置页面</div>} />
        </Route>
        
        {/* 404页面 */}
        <Route path="*" element={<PageNotFound />} />
      </Routes>
    </Suspense>
  );
};

export default AppRoutes; 