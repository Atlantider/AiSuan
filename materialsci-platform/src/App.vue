<template>
  <div id="app">
    <!-- 在你的应用顶部添加调试组件链接 -->
    <div style="padding: 10px; background-color: #f5f5f5; text-align: center;">
      <a href="#" @click.prevent="showDebugger = !showDebugger">
        {{ showDebugger ? '隐藏' : '显示' }}任务调试工具
      </a>
    </div>
    
    <!-- 调试组件 -->
    <TaskDetailsViewer v-if="showDebugger" />
    
    <!-- 原始应用内容 -->
    <router-view />
  </div>
</template>

<script>
import TaskDetailsViewer from './components/TaskDetailsViewer.vue';

export default {
  components: {
    TaskDetailsViewer
  },
  data() {
    return {
      showDebugger: false
    };
  }
};
</script>

<style>
#app {
  font-family: Avenir, Helvetica, Arial, sans-serif;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  color: #2c3e50;
}
</style> 