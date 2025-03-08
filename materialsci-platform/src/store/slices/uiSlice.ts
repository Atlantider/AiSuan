import { createSlice, PayloadAction } from '@reduxjs/toolkit';

export interface UiState {
  siderCollapsed: boolean;
  breadcrumbs: { title: string; path: string }[];
  notifications: {
    id: string;
    type: 'success' | 'info' | 'warning' | 'error';
    message: string;
    description?: string;
    read: boolean;
    createdAt: string;
  }[];
}

const initialState: UiState = {
  siderCollapsed: false,
  breadcrumbs: [],
  notifications: [],
};

export const uiSlice = createSlice({
  name: 'ui',
  initialState,
  reducers: {
    toggleSider: (state) => {
      state.siderCollapsed = !state.siderCollapsed;
    },
    setSiderCollapsed: (state, action: PayloadAction<boolean>) => {
      state.siderCollapsed = action.payload;
    },
    setBreadcrumbs: (state, action: PayloadAction<UiState['breadcrumbs']>) => {
      state.breadcrumbs = action.payload;
    },
    addNotification: (state, action: PayloadAction<Omit<UiState['notifications'][0], 'id' | 'read' | 'createdAt'>>) => {
      const id = Date.now().toString();
      state.notifications = [
        {
          ...action.payload,
          id,
          read: false,
          createdAt: new Date().toISOString(),
        },
        ...state.notifications,
      ];
    },
    markNotificationAsRead: (state, action: PayloadAction<string>) => {
      state.notifications = state.notifications.map((notification) =>
        notification.id === action.payload ? { ...notification, read: true } : notification
      );
    },
    clearAllNotifications: (state) => {
      state.notifications = [];
    },
  },
});

export const {
  toggleSider,
  setSiderCollapsed,
  setBreadcrumbs,
  addNotification,
  markNotificationAsRead,
  clearAllNotifications,
} = uiSlice.actions;

export default uiSlice.reducer; 