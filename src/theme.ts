import { ThemeConfig } from 'antd';

// 主题配置
export const theme: ThemeConfig = {
  token: {
    colorPrimary: '#1C64F2', // 主色调 - 深蓝色
    colorSuccess: '#34D399', // 成功色 - 绿色
    colorWarning: '#FBBF24', // 警告色 - 黄色
    colorError: '#F87171',   // 错误色 - 红色
    colorInfo: '#60A5FA',    // 信息色 - 浅蓝色
    borderRadius: 4,         // 边框圆角
    fontSize: 14,            // 默认字体大小
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

// 全局样式变量
export const globalStyles = {
  contentMaxWidth: 1440,
  contentPadding: '0 24px',
  headerShadow: '0 1px 2px 0 rgba(0, 0, 0, 0.05)',
  cardShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1), 0 1px 2px 0 rgba(0, 0, 0, 0.06)',
  transitionNormal: 'all 0.3s cubic-bezier(0.645, 0.045, 0.355, 1)',
}; 