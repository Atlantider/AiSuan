import axios from 'axios';

// API基础URL
// export const API_BASE_URL = 'http://192.168.199.211:8000';
export const API_BASE_URL = window.location.hostname === 'localhost' 
  ? 'http://localhost:8000'  // 本地开发环境
  : window.location.protocol + '//' + window.location.hostname + ':8000';  // 使用与前端相同的主机名，保留协议
// export const API_PREFIX = '/api';
export const API_PREFIX = '/api/v1';

console.log('当前API基础URL:', API_BASE_URL);
console.log('当前API前缀:', API_PREFIX);
console.log('当前浏览器URL:', window.location.href);

// 创建axios实例
const api = axios.create({
  baseURL: `${API_BASE_URL}${API_PREFIX}`,
  headers: {
    'Content-Type': 'application/json',
    'X-Requested-With': 'XMLHttpRequest'
  },
  timeout: 60000, // 增加超时时间到60秒
  withCredentials: false  // 修改为false，以避免跨域请求时的凭证问题
});

// 请求拦截器
api.interceptors.request.use(
  (config) => {
    console.log('发送请求:', config.url, config);
    const token = localStorage.getItem('token');
    if (token) {
      config.headers.Authorization = `Bearer ${token}`;
    }
    return config;
  },
  (error) => {
    console.error('请求错误:', error);
    return Promise.reject(error);
  }
);

// 响应拦截器
api.interceptors.response.use(
  (response) => {
    console.log('收到响应:', response.status, response.data);
    return response;
  },
  (error) => {
    console.error('响应错误:', error.message);
    if (error.response) {
      console.error('错误状态:', error.response.status);
      console.error('错误数据:', error.response.data);
    } else if (error.request) {
      console.error('没有收到响应:', error.request);
    }
    return Promise.reject(error);
  }
);

export default api; 