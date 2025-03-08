import { ThemeConfig } from 'antd';

// 全局主题和样式设置
export const globalStyles = {
  contentMaxWidth: 1200,
  contentPadding: '24px',
  headerShadow: '0 1px 4px rgba(0, 21, 41, 0.08)',
  cardShadow: '0 1px 2px rgba(0, 0, 0, 0.05)',
  transitionNormal: 'all 0.3s',
  logoText: {
    color: '#fff',
    fontSize: 18,
    fontWeight: 600,
    letterSpacing: 0.5,
  }
};

// 默认主题设置
export const theme = {
  token: {
    colorPrimary: '#1C64F2',
    colorSuccess: '#10B981',
    colorWarning: '#F59E0B',
    colorError: '#EF4444',
    colorInfo: '#3B82F6',
    borderRadius: 4,
    fontFamily: "'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif",
  },
  components: {
    Layout: {
      headerBg: '#FFFFFF',
      headerHeight: 64,
      siderBg: '#FFFFFF',
    },
    Menu: {
      itemHeight: 48,
      activeBarBorderWidth: 3,
      itemSelectedColor: '#1C64F2',
      itemSelectedBg: '#EBF5FF',
    },
    Button: {
      primaryShadow: 'none',
    },
    Card: {
      headerFontSize: 16,
    },
  },
}; 