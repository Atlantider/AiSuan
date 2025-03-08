import React, { useEffect } from 'react';
import { useDispatch } from 'react-redux';
import { setLoginState, setUserInfo, setToken } from './store/slices/userSlice';
import AppRoutes from './routes';

function App() {
  const dispatch = useDispatch();

  // 初始化时检查本地存储中的用户信息
  useEffect(() => {
    const token = localStorage.getItem('token');
    const userInfo = localStorage.getItem('userInfo');

    if (token) {
      dispatch(setToken(token));
      dispatch(setLoginState(true));
      
      if (userInfo) {
        try {
          const parsedUserInfo = JSON.parse(userInfo);
          dispatch(setUserInfo(parsedUserInfo));
        } catch (error) {
          console.error('Failed to parse user info:', error);
        }
      }
    }
  }, [dispatch]);

  return <AppRoutes />;
}

export default App;
