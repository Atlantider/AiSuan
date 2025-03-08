import React, { useState, useEffect } from 'react';
import { Outlet, useLocation, useNavigate } from 'react-router-dom';
import { Layout, Menu, Typography, Avatar, Badge, Dropdown, Button } from 'antd';
import {
  DashboardOutlined,
  ExperimentOutlined,
  DatabaseOutlined,
  ShareAltOutlined,
  SettingOutlined,
  UserOutlined,
  BellOutlined,
  MenuFoldOutlined,
  MenuUnfoldOutlined,
  LogoutOutlined
} from '@ant-design/icons';
import { useSelector, useDispatch } from 'react-redux';
import { RootState } from '@/store';
import { setSiderCollapsed } from '@/store/slices/uiSlice';
import { logout } from '@/store/slices/userSlice';
import { globalStyles } from '@/theme';

const { Header, Sider, Content } = Layout;
const { Title } = Typography;

interface MenuItem {
  key: string;
  icon: React.ReactNode;
  label: string;
  path: string;
  children?: MenuItem[];
}

const DashboardLayout: React.FC = () => {
  const dispatch = useDispatch();
  const navigate = useNavigate();
  const location = useLocation();
  
  const { siderCollapsed, notifications } = useSelector((state: RootState) => state.ui);
  const { userInfo } = useSelector((state: RootState) => state.user);
  
  const [selectedKeys, setSelectedKeys] = useState<string[]>(['dashboard']);
  const [pageTitle, setPageTitle] = useState('个人控制台');
  
  // 菜单项配置
  const menuItems: MenuItem[] = [
    {
      key: 'dashboard',
      icon: <DashboardOutlined />,
      label: '个人控制台',
      path: '/user/dashboard',
    },
    {
      key: 'tasks',
      icon: <ExperimentOutlined />,
      label: '计算任务管理',
      path: '/user/tasks',
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
  
  // 用户菜单
  const userMenuItems = [
    {
      key: 'profile',
      label: '个人信息',
      icon: <UserOutlined />,
    },
    {
      key: 'settings',
      label: '账户设置',
      icon: <SettingOutlined />,
    },
    {
      type: 'divider',
    },
    {
      key: 'logout',
      label: '退出登录',
      icon: <LogoutOutlined />,
    },
  ];
  
  return (
    <Layout style={{ minHeight: '100vh' }}>
      {/* 侧边栏 */}
      <Sider
        width={256}
        collapsible
        collapsed={siderCollapsed}
        onCollapse={(collapsed) => dispatch(setSiderCollapsed(collapsed))}
        style={{ boxShadow: '2px 0 8px 0 rgba(29,35,41,0.05)', zIndex: 999 }}
      >
        {/* 侧边栏顶部Logo */}
        <div style={{ height: 64, padding: '16px', textAlign: 'center', display: 'flex', alignItems: 'center', justifyContent: siderCollapsed ? 'center' : 'flex-start' }}>
          <div style={{ display: 'flex', alignItems: 'center', cursor: 'pointer' }} onClick={() => navigate('/')}>
            <div style={{ width: 32, height: 32, background: '#1C64F2', display: 'flex', alignItems: 'center', justifyContent: 'center', borderRadius: 4 }}>
              <span style={{ color: 'white', fontWeight: 'bold' }}>MS</span>
            </div>
            {!siderCollapsed && <span style={{ marginLeft: 12, fontSize: 16, fontWeight: 600, color: '#1C64F2' }}>用户中心</span>}
          </div>
        </div>
        
        {/* 侧边栏菜单 */}
        <Menu
          mode="inline"
          selectedKeys={selectedKeys}
          style={{ borderRight: 0 }}
          onClick={(e) => handleMenuClick(e.key)}
        >
          {menuItems.map(item => (
            <Menu.Item key={item.key} icon={item.icon}>
              {item.label}
            </Menu.Item>
          ))}
        </Menu>
      </Sider>
      
      <Layout>
        {/* 头部 */}
        <Header style={{ background: '#fff', padding: '0 24px', boxShadow: globalStyles.headerShadow, display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
          <div style={{ display: 'flex', alignItems: 'center' }}>
            <Button 
              type="text" 
              icon={siderCollapsed ? <MenuUnfoldOutlined /> : <MenuFoldOutlined />} 
              onClick={() => dispatch(setSiderCollapsed(!siderCollapsed))}
              style={{ fontSize: 16, marginRight: 16 }}
            />
            <Title level={4} style={{ margin: 0 }}>{pageTitle}</Title>
          </div>
          
          <div style={{ display: 'flex', alignItems: 'center' }}>
            {/* 通知图标 */}
            <Badge count={notifications.filter(n => !n.read).length}>
              <Button type="text" icon={<BellOutlined />} style={{ marginRight: 16 }} />
            </Badge>
            
            {/* 用户信息 */}
            <Dropdown menu={{ items: userMenuItems, onClick: handleUserMenuClick }}>
              <div style={{ display: 'flex', alignItems: 'center', cursor: 'pointer' }}>
                <Avatar src={userInfo.avatar} icon={<UserOutlined />} />
                <span style={{ marginLeft: 8 }}>{userInfo.username || '用户'}</span>
              </div>
            </Dropdown>
          </div>
        </Header>
        
        {/* 内容区 */}
        <Content style={{ margin: '24px 24px 0', overflow: 'initial', background: '#f0f2f5', padding: 24, minHeight: 280 }}>
          <Outlet />
        </Content>
      </Layout>
    </Layout>
  );
};

export default DashboardLayout; 