import React, { useState, useEffect } from 'react';
import { Outlet, useLocation, useNavigate } from 'react-router-dom';
import { Layout, Menu, Button, Drawer, Row, Col, Dropdown, Badge, Avatar, Breadcrumb } from 'antd';
import { MenuFoldOutlined, MenuUnfoldOutlined, UserOutlined, BellOutlined, SettingOutlined, LogoutOutlined } from '@ant-design/icons';
import { useSelector, useDispatch } from 'react-redux';
import { RootState } from '@/store';
import { toggleSider, setBreadcrumbs } from '@/store/slices/uiSlice';
import { logout } from '@/store/slices/userSlice';
import logo from '@/assets/logo.svg';
import { globalStyles } from '@/theme';

const { Header, Sider, Content, Footer } = Layout;

const MainLayout: React.FC = () => {
  const dispatch = useDispatch();
  const navigate = useNavigate();
  const location = useLocation();
  
  const { siderCollapsed, breadcrumbs, notifications } = useSelector((state: RootState) => state.ui);
  const { isLoggedIn, userInfo } = useSelector((state: RootState) => state.user);
  
  const [mobileMenuVisible, setMobileMenuVisible] = useState(false);
  const [activeKey, setActiveKey] = useState('home');
  
  // 根据当前路由更新菜单激活状态
  useEffect(() => {
    const path = location.pathname;
    
    if (path === '/') {
      setActiveKey('home');
    } else if (path.startsWith('/battery')) {
      setActiveKey('battery');
    } else if (path.startsWith('/catalysis')) {
      setActiveKey('catalysis');
    } else if (path.startsWith('/documentation')) {
      setActiveKey('documentation');
    } else if (path.startsWith('/about')) {
      setActiveKey('about');
    }
    
    // 更新面包屑
    updateBreadcrumbs(path);
  }, [location]);
  
  // 更新面包屑导航
  const updateBreadcrumbs = (path: string) => {
    const breadcrumbItems = [{ title: '首页', path: '/' }];
    
    if (path.startsWith('/battery')) {
      breadcrumbItems.push({ title: '电池计算', path: '/battery' });
      
      if (path.includes('/electrode-material')) {
        breadcrumbItems.push({ title: '电极材料计算', path: '/battery/electrode-material' });
        
        if (path.includes('/workbench')) {
          breadcrumbItems.push({ title: '工作台', path: '/battery/electrode-material/workbench' });
        }
      } else if (path.includes('/electrolyte')) {
        breadcrumbItems.push({ title: '电解液计算', path: '/battery/electrolyte' });
        
        if (path.includes('/workbench')) {
          breadcrumbItems.push({ title: '工作台', path: '/battery/electrolyte/workbench' });
        }
      } else if (path.includes('/full-battery')) {
        breadcrumbItems.push({ title: '全电池系统计算', path: '/battery/full-battery' });
        
        if (path.includes('/workbench')) {
          breadcrumbItems.push({ title: '工作台', path: '/battery/full-battery/workbench' });
        }
      }
    } else if (path.startsWith('/catalysis')) {
      breadcrumbItems.push({ title: '催化计算', path: '/catalysis' });
      
      if (path.includes('/surface')) {
        breadcrumbItems.push({ title: '表面催化计算', path: '/catalysis/surface' });
        
        if (path.includes('/workbench')) {
          breadcrumbItems.push({ title: '工作台', path: '/catalysis/surface/workbench' });
        }
      } else if (path.includes('/electro')) {
        breadcrumbItems.push({ title: '电催化计算', path: '/catalysis/electro' });
        
        if (path.includes('/workbench')) {
          breadcrumbItems.push({ title: '工作台', path: '/catalysis/electro/workbench' });
        }
      } else if (path.includes('/photo')) {
        breadcrumbItems.push({ title: '光催化计算', path: '/catalysis/photo' });
        
        if (path.includes('/workbench')) {
          breadcrumbItems.push({ title: '工作台', path: '/catalysis/photo/workbench' });
        }
      }
    } else if (path.startsWith('/documentation')) {
      breadcrumbItems.push({ title: '文档', path: '/documentation' });
    } else if (path.startsWith('/about')) {
      breadcrumbItems.push({ title: '关于我们', path: '/about' });
    }
    
    dispatch(setBreadcrumbs(breadcrumbItems));
  };
  
  // 处理菜单点击
  const handleMenuClick = ({ key }: { key: string }) => {
    switch (key) {
      case 'home':
        navigate('/');
        break;
      case 'battery':
        navigate('/battery');
        break;
      case 'catalysis':
        navigate('/catalysis');
        break;
      case 'documentation':
        navigate('/documentation');
        break;
      case 'about':
        navigate('/about');
        break;
      default:
        break;
    }
  };
  
  // 处理用户菜单点击
  const handleUserMenuClick = ({ key }: { key: string }) => {
    switch (key) {
      case 'profile':
        navigate('/user/dashboard');
        break;
      case 'tasks':
        navigate('/user/tasks');
        break;
      case 'settings':
        navigate('/user/settings');
        break;
      case 'logout':
        // 清除本地存储
        localStorage.removeItem('token');
        localStorage.removeItem('userInfo');
        // 更新状态
        dispatch(logout());
        // 导航到首页
        navigate('/');
        break;
      default:
        break;
    }
  };
  
  // 用户菜单项
  const userMenuItems = [
    {
      key: 'profile',
      label: '个人控制台',
      icon: <UserOutlined />,
    },
    {
      key: 'tasks',
      label: '我的计算任务',
      icon: <SettingOutlined />,
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
  
  // 切换侧边栏
  const toggleSidebar = () => {
    dispatch(toggleSider());
  };
  
  return (
    <Layout style={{ minHeight: '100vh' }}>
      {/* 头部导航 */}
      <Header style={{ background: '#fff', padding: '0 16px', boxShadow: globalStyles.headerShadow, zIndex: 1000, position: 'sticky', top: 0 }}>
        <Row justify="space-between" align="middle">
          {/* 左侧 Logo 和菜单图标 */}
          <Col xs={8} sm={8} md={6} lg={4}>
            <div style={{ display: 'flex', alignItems: 'center', height: '100%' }}>
              <div style={{ marginRight: 16, display: 'flex', alignItems: 'center', cursor: 'pointer' }} onClick={() => navigate('/')}>
                <img src={logo} alt="Logo" style={{ height: 36 }} />
                <span style={{ fontSize: 18, fontWeight: 600, color: '#1C64F2', marginLeft: 8 }}>计算材料科学平台</span>
              </div>
              <Button 
                type="text" 
                icon={mobileMenuVisible ? <MenuUnfoldOutlined /> : <MenuFoldOutlined />}
                onClick={() => setMobileMenuVisible(!mobileMenuVisible)}
                style={{ display: { xs: 'block', lg: 'none' } }}
              />
            </div>
          </Col>
          
          {/* 中间菜单，只在大屏幕显示 */}
          <Col xs={0} sm={0} md={0} lg={14}>
            <Menu
              mode="horizontal"
              selectedKeys={[activeKey]}
              onClick={handleMenuClick}
              style={{ border: 'none', height: '100%', justifyContent: 'center' }}
            >
              <Menu.Item key="home">首页</Menu.Item>
              <Menu.Item key="battery">电池计算</Menu.Item>
              <Menu.Item key="catalysis">催化计算</Menu.Item>
              <Menu.Item key="documentation">文档</Menu.Item>
              <Menu.Item key="about">关于我们</Menu.Item>
            </Menu>
          </Col>
          
          {/* 右侧用户信息 */}
          <Col xs={16} sm={16} md={18} lg={6} style={{ display: 'flex', justifyContent: 'flex-end', alignItems: 'center' }}>
            {isLoggedIn ? (
              <>
                <Badge count={notifications.filter(n => !n.read).length} dot>
                  <Button icon={<BellOutlined />} type="text" style={{ marginRight: 16 }} />
                </Badge>
                <Dropdown menu={{ items: userMenuItems, onClick: handleUserMenuClick }}>
                  <div style={{ display: 'flex', alignItems: 'center', cursor: 'pointer' }}>
                    <Avatar src={userInfo.avatar} icon={<UserOutlined />} />
                    <span style={{ marginLeft: 8, display: { xs: 'none', sm: 'inline' } }}>{userInfo.username || '用户'}</span>
                  </div>
                </Dropdown>
              </>
            ) : (
              <div>
                <Button type="text" onClick={() => navigate('/login')}>登录</Button>
                <Button type="primary" onClick={() => navigate('/register')}>注册</Button>
              </div>
            )}
          </Col>
        </Row>
        
        {/* 移动端菜单抽屉 */}
        <Drawer
          title="菜单"
          placement="left"
          onClose={() => setMobileMenuVisible(false)}
          open={mobileMenuVisible}
          bodyStyle={{ padding: 0 }}
        >
          <Menu
            mode="inline"
            selectedKeys={[activeKey]}
            onClick={(e) => {
              handleMenuClick(e);
              setMobileMenuVisible(false);
            }}
            style={{ borderRight: 0 }}
          >
            <Menu.Item key="home">首页</Menu.Item>
            <Menu.Item key="battery">电池计算</Menu.Item>
            <Menu.Item key="catalysis">催化计算</Menu.Item>
            <Menu.Item key="documentation">文档</Menu.Item>
            <Menu.Item key="about">关于我们</Menu.Item>
          </Menu>
        </Drawer>
      </Header>
      
      {/* 内容区 */}
      <Content style={{ padding: '0 50px', maxWidth: globalStyles.contentMaxWidth, margin: '16px auto', width: '100%' }}>
        {/* 面包屑导航 */}
        <Breadcrumb style={{ margin: '16px 0' }}>
          {breadcrumbs.map((item, index) => (
            <Breadcrumb.Item key={index}>
              <a onClick={() => navigate(item.path)}>{item.title}</a>
            </Breadcrumb.Item>
          ))}
        </Breadcrumb>
        
        {/* 路由出口 */}
        <div style={{ background: '#fff', padding: 24, minHeight: 280, borderRadius: 4, boxShadow: globalStyles.cardShadow }}>
          <Outlet />
        </div>
      </Content>
      
      {/* 页脚 */}
      <Footer style={{ textAlign: 'center', background: '#f0f2f5' }}>
        计算材料科学平台 ©{new Date().getFullYear()} - 提供先进的材料计算服务
      </Footer>
    </Layout>
  );
};

export default MainLayout; 