import React, { useEffect } from 'react';
import { useSelector } from 'react-redux';
import { useNavigate, Navigate } from 'react-router-dom';
import { RootState } from '@/store';
import { Result, Button, Spin } from 'antd';

interface AuthGuardProps {
  children: React.ReactNode;
}

const AuthGuard: React.FC<AuthGuardProps> = ({ children }) => {
  const { isLoggedIn, userInfo } = useSelector((state: RootState) => state.user);
  const navigate = useNavigate();

  // 如果未登录，重定向到登录页面
  if (!isLoggedIn) {
    return (
      <Result
        status="403"
        title="无权访问"
        subTitle="您需要登录才能访问此页面"
        extra={
          <Button type="primary" onClick={() => navigate('/login', { state: { from: window.location.pathname } })}>
            去登录
          </Button>
        }
      />
    );
  }

  // 如果已登录但没有用户信息，显示加载状态
  if (isLoggedIn && (!userInfo || Object.keys(userInfo).length === 0)) {
    return (
      <div style={{ display: 'flex', justifyContent: 'center', alignItems: 'center', height: '100vh' }}>
        <Spin size="large" tip="加载用户信息..." />
      </div>
    );
  }

  // 已登录且有用户信息，显示子组件
  return <>{children}</>;
};

export default AuthGuard; 