import React from 'react';
import { Navigate, useLocation } from 'react-router-dom';
import { useSelector } from 'react-redux';
import { RootState } from '../../store';

interface AuthGuardProps {
  children: React.ReactNode;
}

const AuthGuard: React.FC<AuthGuardProps> = ({ children }) => {
  const { isLoggedIn } = useSelector((state: RootState) => state.user);
  const location = useLocation();

  // 暂时模拟所有用户都已登录，在实际场景中应该使用以下注释的代码
  // if (!isLoggedIn) {
  //   return <Navigate to="/login" state={{ from: location.pathname }} replace />;
  // }
  
  return <>{children}</>;
};

export default AuthGuard; 