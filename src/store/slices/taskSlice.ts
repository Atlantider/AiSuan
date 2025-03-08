import { createSlice, PayloadAction } from '@reduxjs/toolkit';

export interface Task {
  id: string;
  taskType: string;
  module: 'battery' | 'catalysis';
  submodule: string;
  status: 'queued' | 'running' | 'completed' | 'failed' | 'cancelled';
  priority: number;
  createdAt: string;
  startedAt?: string;
  completedAt?: string;
  computeTime?: number;
  progress?: number;
  name?: string;
}

export interface TaskState {
  activeTasks: Task[];
  recentTasks: Task[];
  selectedTask: Task | null;
  loading: boolean;
  error: string | null;
}

const initialState: TaskState = {
  activeTasks: [],
  recentTasks: [],
  selectedTask: null,
  loading: false,
  error: null,
};

export const taskSlice = createSlice({
  name: 'task',
  initialState,
  reducers: {
    setActiveTasks: (state, action: PayloadAction<Task[]>) => {
      state.activeTasks = action.payload;
    },
    setRecentTasks: (state, action: PayloadAction<Task[]>) => {
      state.recentTasks = action.payload;
    },
    setSelectedTask: (state, action: PayloadAction<Task | null>) => {
      state.selectedTask = action.payload;
    },
    addTask: (state, action: PayloadAction<Task>) => {
      state.activeTasks = [action.payload, ...state.activeTasks];
    },
    updateTaskStatus: (state, action: PayloadAction<{ id: string; status: Task['status']; progress?: number }>) => {
      const { id, status, progress } = action.payload;
      
      const updateTask = (task: Task) => {
        if (task.id === id) {
          return { 
            ...task, 
            status, 
            progress: progress !== undefined ? progress : task.progress,
            completedAt: status === 'completed' || status === 'failed' ? new Date().toISOString() : task.completedAt 
          };
        }
        return task;
      };
      
      state.activeTasks = state.activeTasks.map(updateTask);
      state.recentTasks = state.recentTasks.map(updateTask);
      
      if (state.selectedTask && state.selectedTask.id === id) {
        state.selectedTask = updateTask(state.selectedTask);
      }
    },
    setLoading: (state, action: PayloadAction<boolean>) => {
      state.loading = action.payload;
    },
    setError: (state, action: PayloadAction<string | null>) => {
      state.error = action.payload;
    },
  },
});

export const { 
  setActiveTasks, 
  setRecentTasks, 
  setSelectedTask, 
  addTask, 
  updateTaskStatus, 
  setLoading, 
  setError 
} = taskSlice.actions;

export default taskSlice.reducer; 