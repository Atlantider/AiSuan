import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'
import path from 'path'

// https://vitejs.dev/config/
export default defineConfig({
  plugins: [react()],
  resolve: {
    alias: {
      '@': path.resolve(__dirname, './src')
    }
  },
  // 根据环境变量判断：开发环境使用'/'，生产环境使用'/AiSuan/'
  base: process.env.NODE_ENV === 'production' ? '/AiSuan/' : '/',
})
