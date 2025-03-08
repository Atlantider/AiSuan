import { lazy, Suspense } from 'react';
import { Routes, Route, Navigate } from 'react-router-dom';
import { Spin } from 'antd';
import MainLayout from '@/layouts/MainLayout';
import DashboardLayout from '@/layouts/DashboardLayout';
import AuthGuard from '@/components/common/AuthGuard';

// 懒加载页面组件
const Home = lazy(() => import('@/pages/Home'));
const Login = lazy(() => import('@/pages/Login'));
const Register = lazy(() => import('@/pages/Register'));
const NotFound = lazy(() => import('@/pages/NotFound'));

// 电池计算模块
const BatteryIndex = lazy(() => import('@/pages/battery/Index'));
const ElectrodeMaterial = lazy(() => import('@/pages/battery/ElectrodeMaterial'));
const ElectrolyteCalculation = lazy(() => import('@/pages/battery/ElectrolyteCalculation'));
const FullBatterySystem = lazy(() => import('@/pages/battery/FullBatterySystem'));
const BatteryWorkbench = lazy(() => import('@/pages/battery/Workbench'));

// 催化计算模块
const CatalysisIndex = lazy(() => import('@/pages/catalysis/Index'));
const SurfaceCatalysis = lazy(() => import('@/pages/catalysis/SurfaceCatalysis'));
const ElectroCatalysis = lazy(() => import('@/pages/catalysis/ElectroCatalysis'));
const PhotoCatalysis = lazy(() => import('@/pages/catalysis/PhotoCatalysis'));
const CatalysisWorkbench = lazy(() => import('@/pages/catalysis/Workbench'));

// 用户中心
const UserDashboard = lazy(() => import('@/pages/user/Dashboard'));
const TaskManagement = lazy(() => import('@/pages/user/TaskManagement'));
const DataManagement = lazy(() => import('@/pages/user/DataManagement'));
const WorkflowManagement = lazy(() => import('@/pages/user/WorkflowManagement'));
const AccountSettings = lazy(() => import('@/pages/user/AccountSettings'));

// 文档页面
const Documentation = lazy(() => import('@/pages/Documentation'));
const About = lazy(() => import('@/pages/About'));

// 加载组件
const LoadingComponent = () => (
  <div style={{ display: 'flex', justifyContent: 'center', alignItems: 'center', height: '100vh' }}>
    <Spin size="large" tip="页面加载中..." />
  </div>
);

const AppRoutes = () => {
  return (
    <Suspense fallback={<LoadingComponent />}>
      <Routes>
        {/* 公开路由 */}
        <Route path="/login" element={<Login />} />
        <Route path="/register" element={<Register />} />
        
        {/* 主布局路由 */}
        <Route element={<MainLayout />}>
          <Route path="/" element={<Home />} />
          <Route path="/documentation" element={<Documentation />} />
          <Route path="/about" element={<About />} />
          
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
          <Route path="/catalysis/surface" element={<SurfaceCatalysis />} />
          <Route path="/catalysis/surface/workbench" element={<CatalysisWorkbench calculationType="surface" />} />
          <Route path="/catalysis/electro" element={<ElectroCatalysis />} />
          <Route path="/catalysis/electro/workbench" element={<CatalysisWorkbench calculationType="electro" />} />
          <Route path="/catalysis/photo" element={<PhotoCatalysis />} />
          <Route path="/catalysis/photo/workbench" element={<CatalysisWorkbench calculationType="photo" />} />
        </Route>
        
        {/* 用户中心（需要认证） */}
        <Route 
          element={
            <AuthGuard>
              <DashboardLayout />
            </AuthGuard>
          }
        >
          <Route path="/user/dashboard" element={<UserDashboard />} />
          <Route path="/user/tasks" element={<TaskManagement />} />
          <Route path="/user/data" element={<DataManagement />} />
          <Route path="/user/workflows" element={<WorkflowManagement />} />
          <Route path="/user/settings" element={<AccountSettings />} />
        </Route>
        
        {/* 默认路由和错误页面 */}
        <Route path="/404" element={<NotFound />} />
        <Route path="*" element={<Navigate to="/404" replace />} />
      </Routes>
    </Suspense>
  );
};

export default AppRoutes; 