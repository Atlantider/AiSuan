<template>
  <div class="task-details-viewer">
    <h2>任务详情调试器</h2>
    
    <!-- 任务ID输入 -->
    <div class="input-section">
      <a-input-group compact>
        <a-input
          v-model="taskId"
          placeholder="输入任务ID"
          style="width: 200px"
          @pressEnter="loadTaskDetails"
        />
        <a-button type="primary" @click="loadTaskDetails" :loading="loading">
          <a-icon type="search" />查询
        </a-button>
      </a-input-group>
    </div>
    
    <!-- 错误信息 -->
    <a-alert
      v-if="error"
      type="error"
      :message="error"
      show-icon
      style="margin: 16px 0"
    />
    
    <!-- API响应信息 -->
    <div v-if="apiResponse" class="response-section">
      <a-card title="API响应数据" :bordered="false">
        <a-descriptions bordered>
          <a-descriptions-item label="任务ID">{{ apiResponse.id }}</a-descriptions-item>
          <a-descriptions-item label="状态">{{ apiResponse.status }}</a-descriptions-item>
          <a-descriptions-item label="配方ID">{{ apiResponse.formulation_id || '-' }}</a-descriptions-item>
          <a-descriptions-item label="SLURM作业ID">{{ apiResponse.slurm_job_id || '-' }}</a-descriptions-item>
          <a-descriptions-item label="创建时间">{{ apiResponse.created_at ? apiResponse.created_at : '-' }}</a-descriptions-item>
          <a-descriptions-item label="更新时间">{{ apiResponse.updated_at ? apiResponse.updated_at : '-' }}</a-descriptions-item>
        </a-descriptions>
        
        <a-divider orientation="left">详情数据 (details)</a-divider>
        <pre class="details-json">{{ detailsJson }}</pre>
        
        <!-- 调试信息 -->
        <a-collapse v-if="apiResponse.api_version === 'debug'">
          <a-collapse-panel key="1" header="调试信息">
            <p><strong>API版本:</strong> {{ apiResponse.api_version }}</p>
            <p><strong>数据类型:</strong> {{ typeof apiResponse.details }}</p>
            <p><strong>是否为空:</strong> {{ isDetailsEmpty ? '是' : '否' }}</p>
          </a-collapse-panel>
        </a-collapse>
      </a-card>
    </div>
    
    <!-- 结果 -->
    <div v-if="taskResults" class="results-section">
      <a-divider orientation="left">计算结果</a-divider>
      <pre class="details-json">{{ JSON.stringify(taskResults, null, 2) }}</pre>
    </div>
  </div>
</template>

<script>
import axios from 'axios';

export default {
  name: 'TaskDetailsViewer',
  
  data() {
    return {
      taskId: '',
      apiResponse: null,
      taskResults: null,
      loading: false,
      error: null
    };
  },
  
  computed: {
    detailsJson() {
      try {
        if (!this.apiResponse) return '';
        if (!this.apiResponse.details) return '无详情数据';
        
        if (typeof this.apiResponse.details === 'string') {
          try {
            // 尝试解析JSON字符串
            const parsed = JSON.parse(this.apiResponse.details);
            return JSON.stringify(parsed, null, 2);
          } catch (e) {
            // 如果无法解析，直接返回字符串
            return this.apiResponse.details;
          }
        } else {
          // 如果已经是对象，格式化输出
          return JSON.stringify(this.apiResponse.details, null, 2);
        }
      } catch (e) {
        return `解析错误: ${e.message}`;
      }
    },
    
    isDetailsEmpty() {
      if (!this.apiResponse) return true;
      if (!this.apiResponse.details) return true;
      
      if (typeof this.apiResponse.details === 'string') {
        return this.apiResponse.details.trim() === '';
      } else if (typeof this.apiResponse.details === 'object') {
        return Object.keys(this.apiResponse.details).length === 0;
      }
      
      return false;
    }
  },
  
  methods: {
    async loadTaskDetails() {
      if (!this.taskId || isNaN(Number(this.taskId))) {
        this.error = '请输入有效的任务ID';
        return;
      }
      
      this.loading = true;
      this.error = null;
      this.apiResponse = null;
      this.taskResults = null;
      
      try {
        console.log(`获取任务 ${this.taskId} 的详情`);
        const response = await axios.get(`/calculations/${this.taskId}/`);
        console.log('任务详情响应:', response.data);
        this.apiResponse = response.data;
        
        // 如果任务已完成，加载结果
        if (this.apiResponse.status === 'completed') {
          try {
            const resultResponse = await axios.get(`/calculations/${this.taskId}/results/`);
            console.log('任务结果响应:', resultResponse.data);
            this.taskResults = resultResponse.data.results;
          } catch (resultError) {
            console.error('获取结果数据失败:', resultError);
          }
        }
      } catch (error) {
        console.error('获取任务详情失败:', error);
        this.error = '获取任务详情失败，请稍后重试';
      } finally {
        this.loading = false;
      }
    }
  }
};
</script>

<style scoped>
.task-details-viewer {
  max-width: 1000px;
  margin: 0 auto;
  padding: 20px;
}

.input-section {
  margin-bottom: 20px;
}

.response-section, .results-section {
  margin-top: 20px;
}

.details-json {
  background-color: #f5f5f5;
  padding: 12px;
  border-radius: 4px;
  font-family: 'Courier New', monospace;
  white-space: pre-wrap;
  max-height: 400px;
  overflow-y: auto;
}
</style> 