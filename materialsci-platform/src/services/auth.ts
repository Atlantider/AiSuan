import axios from 'axios';
import { API_BASE_URL } from './api';

export interface LoginResponse {
  access: string;
  refresh: string;
}

export const login = async (username: string, password: string): Promise<LoginResponse> => {
  try {
    const response = await axios.post(`${API_BASE_URL}/api/token/`, {
      username,
      password,
    });
    
    const { access, refresh } = response.data;
    localStorage.setItem('token', access);
    localStorage.setItem('refreshToken', refresh);
    
    return response.data;
  } catch (error) {
    console.error('登录失败:', error);
    throw error;
  }
};

export const logout = (): void => {
  localStorage.removeItem('token');
  localStorage.removeItem('refreshToken');
};

export const refreshToken = async (): Promise<string> => {
  try {
    const refresh = localStorage.getItem('refreshToken');
    if (!refresh) {
      throw new Error('No refresh token available');
    }
    
    const response = await axios.post(`${API_BASE_URL}/api/token/refresh/`, {
      refresh,
    });
    
    const { access } = response.data;
    localStorage.setItem('token', access);
    
    return access;
  } catch (error) {
    console.error('刷新令牌失败:', error);
    logout(); // 如果刷新失败，注销用户
    throw error;
  }
};

export const isAuthenticated = (): boolean => {
  return !!localStorage.getItem('token');
}; 