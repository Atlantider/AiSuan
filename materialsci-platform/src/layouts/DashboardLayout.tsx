import React, { useState, useEffect } from 'react';
import { Layout, Menu, Button, Avatar, Badge, Dropdown, Typography } from 'antd';
import { useNavigate, useLocation, Outlet } from 'react-router-dom';
import {
  MenuFoldOutlined,
  MenuUnfoldOutlined,
  DashboardOutlined,
  ExperimentOutlined,
  DatabaseOutlined,
  ShareAltOutlined,
  SettingOutlined,
  UserOutlined,
  LogoutOutlined,
  BellOutlined,
  CalculatorOutlined
} from '@ant-design/icons';
import { useSelector, useDispatch } from 'react-redux';
import { RootState } from '../store';
import { setSiderCollapsed } from '../store/slices/uiSlice';
import { logout } from '../store/slices/userSlice';
import { globalStyles } from '../theme';
import type { MenuProps } from 'antd';

const { Header, Sider, Content } = Layout;
const { Title } = Typography;

type MenuItem = {
  key: string;
  icon: React.ReactNode;
  label: string;
  path: string;
};

const DashboardLayout: React.FC = () => {
  const dispatch = useDispatch();
  const navigate = useNavigate();
  const location = useLocation();
  
  const { siderCollapsed } = useSelector((state: RootState) => state.ui);
  const { userInfo } = useSelector((state: RootState) => state.user);
  
  const [selectedKeys, setSelectedKeys] = useState<string[]>(['dashboard']);
  const [pageTitle, setPageTitle] = useState('个人控制台');
  
  // 模拟通知数据
  const notifications = [
    {
      id: '1',
      type: 'success',
      message: '计算任务完成',
      description: 'Battery-001 计算已完成',
      read: false,
      createdAt: '2023-03-10T10:30:00Z'
    },
    {
      id: '2',
      type: 'info',
      message: '系统维护通知',
      description: '系统将于明天凌晨2点进行维护',
      read: true,
      createdAt: '2023-03-09T16:45:00Z'
    }
  ];
  
  // 菜单项配置
  const menuItems: MenuItem[] = [
    {
      key: 'tasks',
      icon: <ExperimentOutlined />,
      label: '计算任务管理',
      path: '/user/tasks',
    },
    {
      key: 'calculations',
      icon: <CalculatorOutlined />,
      label: '新建计算任务',
      path: '/calculations',
    },
    {
      key: 'dashboard',
      icon: <DashboardOutlined />,
      label: '个人控制台',
      path: '/user/dashboard',
    },
    {
      key: 'data',
      icon: <DatabaseOutlined />,
      label: '数据管理',
      path: '/user/data',
    },
    {
      key: 'workflows',
      icon: <ShareAltOutlined />,
      label: '工作流管理',
      path: '/user/workflows',
    },
    {
      key: 'settings',
      icon: <SettingOutlined />,
      label: '账户设置',
      path: '/user/settings',
    },
  ];
  
  // 根据当前路径更新选中菜单和页面标题
  useEffect(() => {
    const path = location.pathname;
    let currentKey = 'dashboard';
    
    for (const item of menuItems) {
      if (path === item.path) {
        currentKey = item.key;
        setPageTitle(item.label);
        break;
      }
    }
    
    setSelectedKeys([currentKey]);
  }, [location.pathname]);
  
  // 处理菜单点击
  const handleMenuClick = (key: string) => {
    const item = menuItems.find(item => item.key === key);
    if (item) {
      navigate(item.path);
    }
  };
  
  // 处理用户菜单点击
  const handleUserMenuClick = ({ key }: { key: string }) => {
    if (key === 'logout') {
      // 清除本地存储
      localStorage.removeItem('token');
      localStorage.removeItem('userInfo');
      // 更新状态
      dispatch(logout());
      // 导航到首页
      navigate('/');
    }
  };
  
  // 用户菜单项
  const userMenuItems: MenuProps['items'] = [
    {
      key: 'profile',
      label: '个人信息',
      icon: <UserOutlined />
    },
    {
      key: 'settings',
      label: '账户设置',
      icon: <SettingOutlined />
    },
    {
      type: 'divider'
    },
    {
      key: 'logout',
      label: '退出登录',
      icon: <LogoutOutlined />
    }
  ];
  
  return (
    <Layout style={{ minHeight: '100vh' }}>
      {/* 侧边栏 */}
      <Sider trigger={null} collapsible collapsed={siderCollapsed} width={250}>
        <div style={{ display: 'flex', alignItems: 'center', padding: '16px 24px', borderBottom: '1px solid rgba(255, 255, 255, 0.1)' }}>
          <img src="/logo.png" alt="Logo" style={{ width: 32, height: 32, marginRight: siderCollapsed ? 0 : 12 }} />
          {!siderCollapsed && <h1 style={{ ...globalStyles.logoText, margin: 0 }}>材料计算平台</h1>}
        </div>
        <Menu
          theme="dark"
          mode="inline"
          selectedKeys={selectedKeys}
          onClick={({ key }) => handleMenuClick(key)}
          items={menuItems.map(item => ({
            key: item.key,
            icon: item.icon,
            label: item.label
          }))}
        />
      </Sider>
      
      <Layout>
        {/* 头部 */}
        <Header style={{ 
          padding: '0 24px', 
          background: '#fff', 
          display: 'flex', 
          alignItems: 'center', 
          justifyContent: 'space-between',
          boxShadow: '0 1px 4px rgba(0, 21, 41, 0.08)'
        }}>
          <div style={{ display: 'flex', alignItems: 'center' }}>
            {/* 折叠按钮 */}
            <Button
              type="text"
              icon={siderCollapsed ? <MenuUnfoldOutlined /> : <MenuFoldOutlined />}
              onClick={() => dispatch(setSiderCollapsed(!siderCollapsed))}
              style={{ fontSize: '16px', width: 64, height: 64 }}
            />
            
            {/* 页面标题 */}
            <Title level={4} style={{ margin: 0, marginLeft: 12 }}>{pageTitle}</Title>
          </div>
          
          <div style={{ display: 'flex', alignItems: 'center' }}>
            {/* 通知图标 */}
            <Badge count={notifications.filter(n => !n.read).length}>
              <Button type="text" icon={<BellOutlined />} style={{ marginRight: 16 }} />
            </Badge>
            
            {/* 用户信息 */}
            <Dropdown menu={{ items: userMenuItems, onClick: handleUserMenuClick }}>
              <div style={{ display: 'flex', alignItems: 'center', cursor: 'pointer' }}>
                <Avatar src={userInfo?.avatar} icon={<UserOutlined />} />
                <span style={{ marginLeft: 8 }}>{userInfo?.username || '用户'}</span>
              </div>
            </Dropdown>
          </div>
        </Header>
        
        {/* 内容区域 */}
        <Content style={{ margin: '24px', background: '#fff', padding: 24, minHeight: 280 }}>
          <Outlet />
        </Content>
      </Layout>
    </Layout>
  );
};

export default DashboardLayout; 