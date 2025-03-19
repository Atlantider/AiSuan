import axios from 'axios';

// API基础URL
export const API_BASE_URL = 'http://127.0.0.1:8000';
export const API_PREFIX = '/api/v1';

// 创建axios实例
const api = axios.create({
  baseURL: `${API_BASE_URL}${API_PREFIX}`,
  headers: {
    'Content-Type': 'application/json',
  },
});

// 请求拦截器处理认证令牌
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