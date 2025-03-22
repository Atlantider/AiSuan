<template>
  <div class="task-manager">
    <h2>任务管理控制台</h2>
    
    <!-- 过滤和搜索 -->
    <div class="filter-controls">
      <a-form layout="inline">
        <a-form-item label="状态">
          <a-select
            v-model="filters.status"
            style="width: 120px"
            @change="loadTasks"
          >
            <a-select-option value="">全部</a-select-option>
            <a-select-option value="pending">等待中</a-select-option>
            <a-select-option value="submitted">已提交</a-select-option>
            <a-select-option value="running">运行中</a-select-option>
            <a-select-option value="completed">已完成</a-select-option>
            <a-select-option value="failed">失败</a-select-option>
            <a-select-option value="cancelled">已取消</a-select-option>
          </a-select>
        </a-form-item>
        <a-form-item label="配方ID">
          <a-input
            v-model="filters.formulation_id"
            placeholder="配方ID"
            style="width: 150px"
            @pressEnter="loadTasks"
          />
        </a-form-item>
        <a-form-item>
          <a-button type="primary" @click="loadTasks">
            <SearchOutlined />搜索
          </a-button>
        </a-form-item>
        <a-form-item>
          <a-button @click="resetFilters">
            <ReloadOutlined />重置
          </a-button>
        </a-form-item>
      </a-form>
    </div>
    
    <!-- 任务列表 -->
    <a-table
      :columns="columns"
      :data-source="tasks"
      :loading="loading"
      :pagination="pagination"
      @change="handleTableChange"
      row-key="id"
    >
      <!-- 任务ID -->
      <template #id="{ text }">
        <span>{{ text }}</span>
      </template>
      
      <!-- 状态 -->
      <template #status="{ text }">
        <a-tag :color="getStatusColor(text)">
          {{ getStatusText(text) }}
        </a-tag>
      </template>
      
      <!-- 输入文件 -->
      <template #input_file="{ record }">
        <span v-if="record.input_file">{{ record.input_file.name }}</span>
        <span v-else>-</span>
      </template>
      
      <!-- 创建时间 -->
      <template #created_at="{ text }">
        <span>{{ formatDate(text) }}</span>
      </template>
      
      <!-- 操作 -->
      <template #action="{ record }">
        <div class="action-buttons">
          <!-- 查看状态 -->
          <a-button
            type="link"
            size="small"
            @click="viewTaskDetails(record)"
            title="查看详情"
          >
            <EyeOutlined />
          </a-button>
          
          <!-- 重新提交 -->
          <a-button
            type="link"
            size="small"
            :disabled="!['failed', 'cancelled'].includes(record.status)"
            @click="resubmitTask(record)"
            title="重新提交"
          >
            <ReloadOutlined />
          </a-button>
          
          <!-- 修改参数 -->
          <a-button
            type="link"
            size="small"
            @click="cloneWithModifications(record)"
            title="修改参数并重新提交"
          >
            <EditOutlined />
          </a-button>
          
          <!-- 取消任务 -->
          <a-button
            type="link"
            size="small"
            :disabled="!['pending', 'submitted', 'running'].includes(record.status)"
            @click="cancelTask(record)"
            title="取消任务"
          >
            <StopOutlined />
          </a-button>
          
          <!-- 删除任务 -->
          <a-button
            type="link"
            size="small"
            danger
            @click="deleteTask(record)"
            title="删除任务"
          >
            <DeleteOutlined />
          </a-button>
        </div>
      </template>
    </a-table>
    
    <!-- 任务详情弹窗 -->
    <a-modal
      v-model="detailsVisible"
      title="任务详情"
      width="800px"
      :footer="null"
    >
      <div v-if="selectedTask">
        <a-descriptions bordered :column="2">
          <a-descriptions-item label="任务ID">
            {{ selectedTask.id }}
          </a-descriptions-item>
          <a-descriptions-item label="状态">
            <a-tag :color="getStatusColor(selectedTask.status)">
              {{ getStatusText(selectedTask.status) }}
            </a-tag>
          </a-descriptions-item>
          <a-descriptions-item label="配方ID">
            {{ typeof selectedTask.formulation_id === 'object' ? JSON.stringify(selectedTask.formulation_id) : (selectedTask.formulation_id || '-') }}
          </a-descriptions-item>
          <a-descriptions-item label="SLURM作业ID">
            {{ selectedTask.slurm_job_id || '-' }}
          </a-descriptions-item>
          <a-descriptions-item label="创建时间">
            {{ formatDate(selectedTask.created_at) }}
          </a-descriptions-item>
          <a-descriptions-item label="更新时间">
            {{ formatDate(selectedTask.updated_at) }}
          </a-descriptions-item>
          <a-descriptions-item label="输入文件" :span="2">
            {{ selectedTask.input_file ? (typeof selectedTask.input_file === 'object' ? JSON.stringify(selectedTask.input_file) : selectedTask.input_file) : '-' }}
          </a-descriptions-item>
        </a-descriptions>
        
        <a-divider />
        
        <!-- 任务详细日志 -->
        <div v-if="taskDetails">
          <h3>任务日志</h3>
          <a-tabs>
            <a-tab-pane key="1" tab="基本信息">
              <a-descriptions bordered :column="1" v-if="taskDetails.output_dir || taskDetails.cmd">
                <a-descriptions-item label="输出目录" v-if="taskDetails.output_dir">
                  {{ typeof taskDetails.output_dir === 'object' ? JSON.stringify(taskDetails.output_dir) : taskDetails.output_dir }}
                </a-descriptions-item>
                <a-descriptions-item label="命令" v-if="taskDetails.cmd">
                  <pre>{{ typeof taskDetails.cmd === 'object' ? JSON.stringify(taskDetails.cmd, null, 2) : taskDetails.cmd }}</pre>
                </a-descriptions-item>
              </a-descriptions>
              <!-- 没有详细信息时显示提示 -->
              <a-empty v-if="isDetailsEmpty" description="暂无详细信息">
                <template #description>
                  <span>
                    该任务暂无详细信息，可能任务尚未完全初始化。<br />
                    若任务状态为"已完成"或"失败"，请联系管理员检查。
                  </span>
                </template>
              </a-empty>
            </a-tab-pane>
            <a-tab-pane key="2" tab="SLURM日志" v-if="slurm_log">
              <pre class="log-content">{{ slurm_log }}</pre>
            </a-tab-pane>
            <a-tab-pane key="3" tab="LAMMPS日志" v-if="lammps_log">
              <pre class="log-content">{{ lammps_log }}</pre>
            </a-tab-pane>
            <a-tab-pane key="4" tab="错误信息" v-if="taskDetails.error">
              <pre class="log-content error-log">{{ typeof taskDetails.error === 'object' ? JSON.stringify(taskDetails.error, null, 2) : taskDetails.error }}</pre>
              <pre class="log-content error-log" v-if="taskDetails.stderr">{{ typeof taskDetails.stderr === 'object' ? JSON.stringify(taskDetails.stderr, null, 2) : taskDetails.stderr }}</pre>
            </a-tab-pane>
          </a-tabs>
        </div>
        
        <!-- 没有详情数据时显示加载提示 -->
        <div v-else>
          <a-empty description="正在加载任务详情...">
            <template #description>
              <span>
                系统正在获取任务详情数据，请稍候...<br />
                如长时间未显示，请刷新页面重试。
              </span>
            </template>
            <a-button type="primary" @click="viewTaskDetails(selectedTask)">
              <ReloadOutlined />重新加载
            </a-button>
          </a-empty>
        </div>
        
        <div style={{ marginTop: '20px', textAlign: 'center' }}>
          <a-button 
            type="primary" 
            @click="handleViewFullDetails(selectedTask)"
          >
            查看完整详情
          </a-button>
        </div>
        
        <!-- 结果分析 -->
        <div v-if="selectedTask.status === 'completed' && taskResults">
          <a-divider />
          <h3>计算结果</h3>
          <a-descriptions bordered :column="2">
            <template v-for="(value, key) in taskResults">
              <a-descriptions-item :label="key" :key="key">
                {{ typeof value === 'object' ? JSON.stringify(value, null, 2) : value }}
              </a-descriptions-item>
            </template>
          </a-descriptions>
        </div>
      </div>
    </a-modal>
    
    <!-- 参数修改弹窗 -->
    <a-modal
      v-model="modifyVisible"
      title="修改参数"
      :confirm-loading="submitting"
      @ok="submitModifiedTask"
    >
      <a-form-model :model="modifiedParams" :label-col="{ span: 6 }" :wrapper-col="{ span: 18 }">
        <a-form-model-item label="配方ID">
          <a-input v-model="modifiedParams.formulation_id" />
        </a-form-model-item>
        
        <!-- TODO: 添加更多可修改的参数 -->
        <a-alert
          type="info"
          show-icon
          message="注意：修改参数将创建一个新的任务"
        />
      </a-form-model>
    </a-modal>
  </div>
</template>

<script>
import axios from 'axios';
import moment from 'moment';
import { 
  EyeOutlined, 
  ReloadOutlined, 
  EditOutlined, 
  StopOutlined, 
  DeleteOutlined,
  SearchOutlined
} from '@ant-design/icons-vue';

export default {
  name: 'TaskManager',
  props: {
    // 任务点击事件处理函数
    onTaskClick: {
      type: Function,
      default: null
    }
  },
  components: {
    EyeOutlined,
    ReloadOutlined,
    EditOutlined,
    StopOutlined,
    DeleteOutlined,
    SearchOutlined
  },
  
  data() {
    return {
      tasks: [],
      loading: false,
      filters: {
        status: '',
        formulation_id: '',
      },
      pagination: {
        current: 1,
        pageSize: 10,
        total: 0,
      },
      columns: [
        {
          title: '任务ID',
          dataIndex: 'id',
          key: 'id',
          scopedSlots: { customRender: 'id' },
        },
        {
          title: '状态',
          dataIndex: 'status',
          key: 'status',
          scopedSlots: { customRender: 'status' },
        },
        {
          title: '配方ID',
          dataIndex: 'formulation_id',
          key: 'formulation_id',
        },
        {
          title: '输入文件',
          dataIndex: 'input_file',
          key: 'input_file',
          scopedSlots: { customRender: 'input_file' },
        },
        {
          title: '创建时间',
          dataIndex: 'created_at',
          key: 'created_at',
          scopedSlots: { customRender: 'created_at' },
        },
        {
          title: '操作',
          key: 'action',
          scopedSlots: { customRender: 'action' },
        },
      ],
      detailsVisible: false,
      modifyVisible: false,
      selectedTask: null,
      taskDetails: null,
      taskResults: null,
      modifiedParams: {
        formulation_id: '',
      },
      submitting: false,
      refreshTimer: null,
      slurm_log: null,
      lammps_log: null,
    };
  },
  
  computed: {
    /**
     * 检查任务详情是否为空
     */
    isDetailsEmpty() {
      if (!this.taskDetails) return true;
      
      // 检查对象是否为空（没有属性或只有空属性）
      const hasValidProperties = Object.keys(this.taskDetails).some(key => {
        const value = this.taskDetails[key];
        // 属性有效的条件：不为null、不为undefined、数组不为空、字符串不为空、对象不为空
        if (value === null || value === undefined) return false;
        if (Array.isArray(value) && value.length === 0) return false;
        if (typeof value === 'string' && value.trim() === '') return false;
        if (typeof value === 'object' && Object.keys(value).length === 0) return false;
        return true;
      });
      
      return !hasValidProperties;
    }
  },
  
  mounted() {
    this.loadTasks();
    // 自动刷新活动任务
    this.refreshTimer = setInterval(() => {
      if (this.tasks.some(task => ['pending', 'submitted', 'running'].includes(task.status))) {
        this.loadTasks(false);
      }
    }, 30000); // 每30秒刷新一次
  },
  
  beforeDestroy() {
    if (this.refreshTimer) {
      clearInterval(this.refreshTimer);
    }
  },
  
  methods: {
    /**
     * 加载任务列表
     */
    async loadTasks(showLoading = true) {
      if (showLoading) {
        this.loading = true;
      }
      
      try {
        const params = {
          page: this.pagination.current,
          page_size: this.pagination.pageSize,
        };
        
        if (this.filters.status) {
          params.status = this.filters.status;
        }
        
        if (this.filters.formulation_id) {
          params.formulation_id = this.filters.formulation_id;
        }
        
        const response = await axios.get('/api/v1/calculations/', { params });
        console.log('获取任务列表响应:', response.data);
        
        // 安全处理API响应
        if (response.data) {
          // 如果响应中有tasks字段并且是数组
          if (response.data.tasks && Array.isArray(response.data.tasks)) {
            this.tasks = response.data.tasks;
            this.pagination.total = response.data.total || response.data.tasks.length;
          }
          // 如果响应本身是数组
          else if (Array.isArray(response.data)) {
            this.tasks = response.data;
            this.pagination.total = response.data.length;
          }
          // 如果响应中有results字段并且是数组
          else if (response.data.results && Array.isArray(response.data.results)) {
            this.tasks = response.data.results;
            this.pagination.total = response.data.count || response.data.results.length;
          }
          // 空数组作为后备
          else {
            console.warn('任务列表数据格式不符合预期:', response.data);
            this.tasks = [];
            this.pagination.total = 0;
          }
        } else {
          this.tasks = [];
          this.pagination.total = 0;
        }
      } catch (error) {
        console.error('加载任务列表出错:', error);
        this.$message.error('加载任务列表失败');
        this.tasks = [];
        this.pagination.total = 0;
      } finally {
        this.loading = false;
      }
    },
    
    /**
     * 处理表格分页、排序、筛选
     */
    handleTableChange(pagination) {
      this.pagination.current = pagination.current;
      this.loadTasks();
    },
    
    /**
     * 重置过滤器
     */
    resetFilters() {
      this.filters.status = '';
      this.filters.formulation_id = '';
      this.loadTasks();
    },
    
    /**
     * 获取状态颜色
     */
    getStatusColor(status) {
      const colors = {
        pending: 'orange',
        submitted: 'blue',
        running: 'cyan',
        completed: 'green',
        failed: 'red',
        cancelled: 'gray',
      };
      
      return colors[status] || 'default';
    },
    
    /**
     * 获取状态文本
     */
    getStatusText(status) {
      const texts = {
        pending: '等待中',
        submitted: '已提交',
        running: '运行中',
        completed: '已完成',
        failed: '失败',
        cancelled: '已取消',
      };
      
      return texts[status] || status;
    },
    
    /**
     * 格式化日期
     */
    formatDate(date) {
      if (!date) return '-';
      return moment(date).format('YYYY-MM-DD HH:mm:ss');
    },
    
    /**
     * 查看任务详情
     */
    async viewTaskDetails(task) {
      this.selectedTask = task;
      this.detailsVisible = true;
      this.taskDetails = null;
      this.taskResults = null;
      this.slurm_log = null;
      this.lammps_log = null;
      
      console.log('查看任务详情:', task);
      
      try {
        // 获取任务详情
        console.log(`获取任务 ${task.id} 的详情`);
        const response = await axios.get(`/api/v1/calculations/${task.id}/`);
        console.log('API响应:', response.data);
        
        // 安全处理details数据
        if (response.data) {
          // 如果details是一个对象，直接使用
          if (response.data.details && typeof response.data.details === 'object') {
            this.taskDetails = response.data.details;
          } 
          // 如果details是整个response或response中的某个字段
          else if (typeof response.data === 'object') {
            // 尝试提取有用的字段作为details
            this.taskDetails = {
              output_dir: response.data.output_dir,
              cmd: response.data.cmd,
              error: response.data.error_message || response.data.error,
              stderr: response.data.stderr
            };
          }
        }
        
        console.log('处理后的任务详情数据:', this.taskDetails);
        
        // 如果任务已完成，获取结果
        if (task.status === 'completed') {
          try {
            const resultResponse = await axios.get(`/api/v1/calculations/${task.id}/results/`);
            console.log('任务结果:', resultResponse.data);
            
            // 确保结果是一个对象
            if (resultResponse.data) {
              if (typeof resultResponse.data === 'object') {
                // 如果结果直接在data中
                if (!resultResponse.data.results) {
                  this.taskResults = resultResponse.data;
                }
                // 如果结果在results字段中且是对象
                else if (typeof resultResponse.data.results === 'object') {
                  this.taskResults = resultResponse.data.results;
                }
                // 如果results不是对象，包装为对象
                else {
                  this.taskResults = { 'result': resultResponse.data.results };
                }
              } else {
                // 如果整个响应不是对象，包装为对象
                this.taskResults = { 'raw': resultResponse.data };
              }
            } else {
              console.warn('未获取到结果数据');
              this.taskResults = { 'message': '未找到结果数据' };
            }
          } catch (error) {
            console.error('获取任务结果出错:', error);
            this.taskResults = { 'error': '获取结果失败' };
          }
        }
        
        // 获取日志文件
        if (task.slurm_job_id) {
          try {
            // SLURM日志
            console.log(`获取SLURM作业 ${task.slurm_job_id} 的日志`);
            const slurmLogResponse = await axios.get(`/api/v1/logs/slurm/${task.slurm_job_id}/`);
            if (slurmLogResponse.data && slurmLogResponse.data.content) {
              this.slurm_log = slurmLogResponse.data.content;
            } else if (typeof slurmLogResponse.data === 'string') {
              this.slurm_log = slurmLogResponse.data;
            }
          } catch (error) {
            console.error('获取SLURM日志出错:', error);
            this.slurm_log = "获取SLURM日志失败: " + (error.message || '未知错误');
          }
        }
        
        // LAMMPS日志
        if (task.status !== 'pending' && task.status !== 'submitted') {
          try {
            console.log(`获取任务 ${task.id} 的LAMMPS日志`);
            const lammpsLogResponse = await axios.get(`/api/v1/logs/lammps/${task.id}/`);
            if (lammpsLogResponse.data && lammpsLogResponse.data.content) {
              this.lammps_log = lammpsLogResponse.data.content;
            } else if (typeof lammpsLogResponse.data === 'string') {
              this.lammps_log = lammpsLogResponse.data;
            }
          } catch (error) {
            console.error('获取LAMMPS日志出错:', error);
            this.lammps_log = "获取LAMMPS日志失败: " + (error.message || '未知错误');
          }
        }
      } catch (error) {
        console.error('获取任务详情出错:', error);
        this.$message.error('获取任务详情失败');
      }
    },
    
    /**
     * 查看完整详情（跳转到详情页面）
     */
    handleViewFullDetails(task) {
      this.detailsVisible = false;
      // 如果有外部传入的处理函数，则调用
      if (this.onTaskClick && typeof this.onTaskClick === 'function') {
        this.onTaskClick(task.id);
      }
    },
    
    /**
     * 重新提交任务
     */
    async resubmitTask(task) {
      if (this.submitting) return;
      
      this.submitting = true;
      
      try {
        // 调用重新提交任务API
        console.log(`重新提交任务: ${task.id}`);
        const response = await axios.post(`/api/v1/calculations/${task.id}/restart/`);
        
        console.log('重新提交响应:', response.data);
        this.$message.success('任务已重新提交');
        this.loadTasks();
      } catch (error) {
        console.error('重新提交任务失败:', error);
        this.$message.error('重新提交任务失败，请稍后重试');
      } finally {
        this.submitting = false;
      }
    },
    
    /**
     * 修改参数并重新提交
     */
    cloneWithModifications(task) {
      this.selectedTask = task;
      this.modifiedParams = {
        formulation_id: task.formulation_id || '',
      };
      this.modifyVisible = true;
    },
    
    /**
     * 提交修改后的任务
     */
    async submitModifiedTask() {
      if (this.submitting) return;
      
      this.submitting = true;
      
      try {
        // 调用克隆任务API
        console.log(`提交修改后的任务: ${this.selectedTask.id}, 参数:`, this.modifiedParams);
        const response = await axios.post(`/api/v1/calculations/${this.selectedTask.id}/clone/`, this.modifiedParams);
        
        console.log('提交修改后的任务响应:', response.data);
        this.$message.success('修改后的任务已提交');
        this.modifyVisible = false;
        this.loadTasks();
      } catch (error) {
        console.error('提交修改后的任务失败:', error);
        this.$message.error('提交修改后的任务失败，请稍后重试');
      } finally {
        this.submitting = false;
      }
    },
    
    /**
     * 取消任务
     */
    async cancelTask(task) {
      try {
        // 调用取消任务API
        console.log(`取消任务: ${task.id}`);
        const response = await axios.post(`/api/v1/calculations/${task.id}/cancel/`);
        
        console.log('取消任务响应:', response.data);
        this.$message.success('任务已取消');
        this.loadTasks();
      } catch (error) {
        console.error('取消任务失败:', error);
        this.$message.error('取消任务失败，请稍后重试');
      }
    },
    
    /**
     * 删除任务
     */
    async deleteTask(task) {
      try {
        // 调用删除任务API
        console.log(`删除任务: ${task.id}`);
        const response = await axios.delete(`/api/v1/calculations/${task.id}/`);
        
        console.log('删除任务响应:', response.data);
        this.$message.success('任务已删除');
        this.loadTasks();
      } catch (error) {
        console.error('删除任务失败:', error);
        this.$message.error('删除任务失败，请稍后重试');
      }
    },
  },
};
</script>

<style scoped>
.task-manager {
  padding: 20px;
  background-color: #ffffff;
  border-radius: 8px;
  box-shadow: 0 1px 3px rgba(0, 0, 0, 0.05);
}

.filter-controls {
  margin-bottom: 20px;
  padding: 16px;
  background-color: #f9f9f9;
  border-radius: 6px;
}

.action-buttons {
  display: flex;
  justify-content: space-between;
}

.log-content {
  max-height: 300px;
  overflow-y: auto;
  padding: 10px;
  background: #f5f5f5;
  border-radius: 4px;
  font-family: monospace;
  white-space: pre-wrap;
}

.error-log {
  background: #fff2f0;
  color: #f5222d;
}
</style> 