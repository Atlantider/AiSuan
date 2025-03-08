import { createSlice, PayloadAction } from '@reduxjs/toolkit';

export interface UserState {
  isLoggedIn: boolean;
  userInfo: {
    id?: string;
    username?: string;
    email?: string;
    role?: 'admin' | 'advanced' | 'regular';
    organization?: string;
    avatar?: string;
  };
  token: string | null;
}

const initialState: UserState = {
  isLoggedIn: false,
  userInfo: {},
  token: null,
};

export const userSlice = createSlice({
  name: 'user',
  initialState,
  reducers: {
    setUserInfo: (state, action: PayloadAction<UserState['userInfo']>) => {
      state.userInfo = action.payload;
    },
    setLoginState: (state, action: PayloadAction<boolean>) => {
      state.isLoggedIn = action.payload;
    },
    setToken: (state, action: PayloadAction<string | null>) => {
      state.token = action.payload;
    },
    logout: (state) => {
      state.isLoggedIn = false;
      state.userInfo = {};
      state.token = null;
    },
  },
});

export const { setUserInfo, setLoginState, setToken, logout } = userSlice.actions;

export default userSlice.reducer; 