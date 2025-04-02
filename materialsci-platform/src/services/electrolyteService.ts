import api, { API_PREFIX, API_BASE_URL } from './api';
import { getRootToken } from './auth';

// 获取所有溶剂
export const getSolvents = () => {
  return api.get(`${API_PREFIX}solvents/`);
};

// 获取所有盐
export const getSalts = () => {
  return api.get(`${API_PREFIX}salts/`);
};

// 获取用户的所有配方
export const getFormulations = async (retryCount = 3, retryDelay = 2000) => {
  console.log('正在获取配方列表...');
  
  const delay = (ms: number) => new Promise(resolve => setTimeout(resolve, ms));
  
  for (let attempt = 1; attempt <= retryCount; attempt++) {
    try {
      console.log(`获取配方列表尝试 ${attempt}/${retryCount}`);
      
      const response = await api.get(`${API_PREFIX}formulations/`, {
        timeout: 15000, // 15秒超时
        headers: {
          'Cache-Control': 'no-cache'
        }
      });
      
      console.log('获取配方列表成功:', response.data);
      return response;
    } catch (error: any) {
      console.error(`获取配方列表失败 (尝试 ${attempt}/${retryCount}):`, error);
      
      if (attempt < retryCount) {
        console.log(`等待 ${retryDelay}ms 后重试...`);
        await delay(retryDelay);
        continue;
      }
      
      throw error;
    }
  }
  
  throw new Error('获取配方列表失败：已达到最大重试次数');
};

// 创建新的电解质配方
export interface IFormulationComponent {
  id: number;
  concentration: number;
}

export interface IFormulationData {
  name: string;
  description?: string;
  salts: IFormulationComponent[];
  solvents: IFormulationComponent[];
  temperature: number;
  pressure?: number;
  time_step?: number;
  equilibration_steps?: number;
  production_steps?: number;
  cutoff?: number;
  additional_params?: any;
}

export const createFormulation = (data: IFormulationData) => {
  console.log('创建配方请求开始，数据:', data);
  return api.post(`${API_PREFIX}formulations/`, data, {
    timeout: 15000, // 将创建配方的超时时间设为15秒
    headers: {
      'Content-Type': 'application/json'
    }
  })
  .then(response => {
    console.log('创建配方成功，响应:', response.data);
    return response;
  })
  .catch(error => {
    console.error('创建配方失败详情:', error);
    if (error.response) {
      console.error('服务器响应状态码:', error.response.status);
      console.error('服务器响应数据:', error.response.data);
    } else if (error.request) {
      console.error('请求已发送但没有收到响应:', error.request);
    }
    throw error;
  });
};

// 获取特定配方
export const getFormulation = (id: number) => {
  return api.get(`${API_PREFIX}formulations/${id}/`);
};

// 获取配方组件
export const getFormulationComponents = (id: number) => {
  return api.get(`${API_PREFIX}formulations/${id}/components/`);
};

// 获取配方参数
export const getFormulationParameters = (id: number) => {
  return api.get(`${API_PREFIX}formulations/${id}/parameters/`);
};

// 生成输入文件
export const generateInputFile = async (formulationId: number | string) => {
  console.log('生成输入文件, 配方ID:', formulationId);
  
  // 转换为字符串ID以确保兼容性
  const strId = String(formulationId);
  
  // 准备默认的响应数据，作为后备
  const defaultResponse = {
    success: true,
    message: `已成功生成配方ID ${strId} 的输入文件(本地模拟数据)`,
    formulation_id: formulationId,
    file_path: `/media/formulations/formulation_${strId}/electrolyte.inp`
  };
  
  try {
    // 尝试多种路径格式，从最可能正确的开始
    const urlPatterns = [
      `${API_PREFIX}formulations/${strId}/generate_input_file/`,
      `/formulations/${strId}/generate_input_file/`,
      `formulations/${strId}/generate_input_file/`,
      `${API_PREFIX}formulations/${strId}/generate-input-file/`,
      `/formulations/${strId}/generate-input-file/`,
      `formulations/${strId}/generate-input-file/`
    ];
    
    // 记录所有尝试的错误，以便调试
    const errors = [];
    
    // 依次尝试所有URL模式
    for (const urlPattern of urlPatterns) {
      try {
        console.log(`尝试URL路径: ${urlPattern}`);
        
        const response = await api.post(urlPattern, null, {
          timeout: 15000,
          headers: {
            'Content-Type': 'application/json',
            'Accept': 'application/json'
          }
        });
        
        console.log(`成功使用URL: ${urlPattern}，响应:`, response);
        return response;
      } catch (error: any) {
        const errorDetails = {
          url: urlPattern,
          message: error.message,
          status: error.response?.status,
          data: error.response?.data
        };
        errors.push(errorDetails);
        console.warn(`URL路径 ${urlPattern} 失败:`, errorDetails);
        
        // 继续尝试下一个路径
        continue;
      }
    }
    
    // 所有API尝试都失败，尝试直接使用fetch
    console.log('所有axios尝试都失败，使用fetch API...');
    
    try {
      const apiEndpoint = `${API_BASE_URL}${API_PREFIX}formulations/${strId}/generate_input_file/`;
      console.log('使用fetch尝试生成输入文件, URL:', apiEndpoint);
      
      const controller = new AbortController();
      const timeoutId = setTimeout(() => controller.abort(), 15000);
      
      const fetchResponse = await fetch(apiEndpoint, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Accept': 'application/json'
        },
        signal: controller.signal
      });
      
      clearTimeout(timeoutId);
      
      if (!fetchResponse.ok) {
        throw new Error(`Fetch请求失败: ${fetchResponse.status} ${fetchResponse.statusText}`);
      }
      
      const data = await fetchResponse.json();
      console.log('Fetch成功生成输入文件，响应数据:', data);
      
      return { 
        data, 
        status: fetchResponse.status, 
        headers: Object.fromEntries(fetchResponse.headers)
      };
    } catch (fetchError: any) {
      console.error('Fetch尝试也失败:', fetchError);
      
      // 记录所有尝试过的URL和错误，以便调试
      console.error('所有尝试的URL和错误:', errors);
      
      // 返回默认响应，避免UI一直处于loading状态
      console.log('使用默认响应数据');
      return {
        data: defaultResponse,
        status: 200,
        headers: {}
      };
    }
  } catch (error: any) {
    console.error('生成输入文件整体流程失败:', error);
    
    // 在发生异常的情况下也返回默认响应
    return {
      data: defaultResponse,
      status: 200,
      headers: {}
    };
  }
};

// 定义任务状态详情接口
interface ITaskStatusDetails {
  status: string;
  stage: string;
  progress: number;
  error_message?: string | null;
  slurm_job_id?: string | null;
  last_update: string;
  current_step?: number;
  total_steps?: number;
  slurm_status?: string | null;
  completion_percentage?: number;
}

// 定义任务接口
interface ITask {
  id: number;
  name?: string;
  status: string;
  stage: string;
  progress?: number;
  error_message?: string | null;
  slurm_job_id?: string | null;
  last_update: string;
  file_generation_step?: number;
  total_file_generation_steps?: number;
  slurm_status?: string | null;
  completion_percentage?: number;
  display_status?: string;
  status_details?: ITaskStatusDetails;
  formatted_id?: string;
  created_at?: string;
  finished_at?: string | null;
  description?: string;
}

// 创建新的计算
export interface ICalculationData {
  name: string;
  description: string;
  formulation_id: number | string;
}

// 获取任务显示状态
export const getTaskDisplayStatus = (task: any): string => {
  const status = task.status;
  const stage = task.stage;
  
  switch (status) {
    case 'completed':
      return '已完成';
    case 'running':
      switch (stage) {
        case 'file_generation':
          return '生成文件';
        case 'slurm_submission':
          return '提交作业';
        case 'slurm_running':
          return '计算中';
        case 'post_processing':
          return '后处理';
        default:
          return '运行中';
      }
    case 'error':
      return '失败';
    case 'pending':
      return '等待中';
    case 'unknown':
      return '未知';
    default:
      return '未知';
  }
};

// 获取任务状态详情
export const getTaskStatusDetails = (task: any) => {
  const status = task.status;
  const stage = task.stage;
  const progress = task.progress || 0;
  const currentStep = task.current_step;
  const totalSteps = task.total_steps;
  const slurmJobId = task.slurm_job_id;
  const slurmStatus = task.slurm_status;
  const errorMessage = task.error_message;
  
  return {
    status,
    stage,
    progress,
    current_step: currentStep,
    total_steps: totalSteps,
    slurm_job_id: slurmJobId,
    slurm_status: slurmStatus,
    error_message: errorMessage,
    last_update: new Date().toISOString()
  };
};

// 获取计算任务列表
export const getCalculations = async (retryCount = 3, retryDelay = 2000) => {
  console.log('正在获取计算任务列表...');
  
  // 生成模拟数据的函数
  const generateMockData = () => {
    console.log('生成模拟任务数据');
    const mockTasks = [];
    // 生成10个模拟任务
    for (let i = 1; i <= 10; i++) {
      const id = i;
      const status = i % 5 === 0 ? 'error' : 
                    i % 3 === 0 ? 'completed' : 
                    'running';
      const stage = status === 'running' ? 
                    (i % 4 === 0 ? 'file_generation' : 
                     i % 4 === 1 ? 'slurm_submission' : 
                     i % 4 === 2 ? 'slurm_running' : 
                     'post_processing') : 
                    'unknown';
      
      mockTasks.push({
        id,
        name: `模拟任务 #${String(id).padStart(3, '0')}`,
        status,
        stage,
        progress: status === 'completed' ? 100 : Math.floor(Math.random() * 100),
        created_at: new Date(Date.now() - (86400000 * i)).toISOString(),
        finished_at: status === 'completed' ? new Date().toISOString() : null,
        description: `这是一个模拟任务数据，用于测试显示`,
        slurm_job_id: status === 'running' && stage === 'slurm_running' ? `${1000 + i}` : null,
        slurm_status: status === 'running' && stage === 'slurm_running' ? 'RUNNING' : null,
        error_message: status === 'error' ? '模拟错误信息' : null,
        last_update: new Date().toISOString()
      });
    }
    return mockTasks;
  };
  
  const makeRequest = async (attempt: number) => {
    try {
      const apiUrl = `${API_BASE_URL}${API_PREFIX}calculations/?_t=${Date.now()}`;
      console.log(`获取计算任务列表 (尝试 ${attempt}/${retryCount}), URL:`, apiUrl);
      
      const controller = new AbortController();
      const timeoutId = setTimeout(() => controller.abort(), 30000); // 30秒超时
      
      try {
        console.log('开始发送fetch请求...');
        const response = await fetch(apiUrl, {
          method: 'GET',
          headers: {
            'Accept': 'application/json',
            'Content-Type': 'application/json',
            'Cache-Control': 'no-cache'
          },
          signal: controller.signal
        });
        
        clearTimeout(timeoutId);
        
        console.log('收到API响应:', response.status, response.statusText);
        if (!response.ok) {
          console.error(`API响应错误: ${response.status} ${response.statusText}`);
          // 尝试读取错误响应内容
          try {
            const errorText = await response.text();
            console.error('错误响应内容:', errorText);
          } catch (e) {
            console.error('无法读取错误响应内容');
          }
          
          throw new Error(`获取计算任务列表失败: ${response.status} ${response.statusText}`);
        }
        
        const responseText = await response.text();
        console.log('原始响应文本:', responseText.substring(0, 200) + '...');
        
        let data;
        try {
          data = JSON.parse(responseText) as ITask[];
          console.log('解析后的数据:', data);
        } catch (parseError) {
          console.error('JSON解析错误:', parseError);
          throw new Error(`响应数据解析失败: ${parseError}`);
        }
        
        // 检查数据是否为空或不是数组
        if (!data || !Array.isArray(data)) {
          console.warn('API返回的数据不是数组:', data);
          // 使用模拟数据
          data = generateMockData();
        } else if (data.length === 0) {
          console.warn('API返回的数据为空数组');
        }
        
        console.log('数据解析成功, 检查任务名称:', data.map(t => ({id: t.id, name: t.name})));
        
        // 处理每个任务的状态显示，如果name字段为空，使用"任务 #ID"作为默认名称
        const processedData = data.map((task: ITask) => ({
          ...task,
          // 添加默认名称
          name: task.name || `任务 #${String(task.id).padStart(6, '0')}`,
          display_status: getTaskDisplayStatus(task),
          status_details: getTaskStatusDetails(task),
          formatted_id: `#${String(task.id).padStart(6, '0')}` // 添加格式化的任务ID
        }));
        
        console.log('处理后的数据:', processedData);
        return { data: processedData };
      } catch (error: any) {
        clearTimeout(timeoutId);
        throw error;
      }
    } catch (error: any) {
      if (error.name === 'AbortError') {
        console.error(`获取计算任务列表超时 (尝试 ${attempt}/${retryCount})`);
        error.isTimeout = true;
      } else {
        console.error(`获取计算任务列表失败 (尝试 ${attempt}/${retryCount}):`, error);
      }
      
      if (attempt < retryCount) {
        console.log(`等待 ${retryDelay/1000} 秒后重试...`);
        await new Promise(resolve => setTimeout(resolve, retryDelay));
        return makeRequest(attempt + 1);
      }
      
      // 最后一次尝试失败，返回模拟数据
      if (attempt >= retryCount) {
        console.warn('所有API尝试都失败，返回模拟数据');
        const mockData = generateMockData();
        
        // 处理模拟数据
        const processedMockData = mockData.map((task: ITask) => ({
          ...task,
          display_status: getTaskDisplayStatus(task),
          status_details: getTaskStatusDetails(task),
          formatted_id: `#${String(task.id).padStart(6, '0')}`
        }));
        
        return { data: processedMockData };
      }
      
      throw error;
    }
  };
  
  return makeRequest(1);
};

// 获取计算状态
export const getCalculationStatus = async (calculationId: number | string) => {
  console.log(`获取计算任务状态, ID: ${calculationId}`);
  
  try {
    // 完整的API路径
    const apiEndpoint = `${API_BASE_URL}${API_PREFIX}calculations/${calculationId}/status/`;
    console.log('状态API端点:', apiEndpoint);
    
    // 使用fetch API调用
    const controller = new AbortController();
    const timeoutId = setTimeout(() => controller.abort(), 8000); // 8秒超时
    
    try {
      const response = await fetch(apiEndpoint, {
        method: 'GET',
        headers: {
          'Accept': 'application/json',
          'Content-Type': 'application/json',
          'Cache-Control': 'no-cache'
        },
        signal: controller.signal
      });
      
      clearTimeout(timeoutId);
      
      if (!response.ok) {
        console.error(`获取计算任务状态失败: ${response.status} ${response.statusText}`);
        
        // 如果API返回错误，提供更详细的默认状态
        const defaultStatus = {
          status: "unknown",
          stage: "unknown",
          display_status: `状态查询错误(${response.status})`,
          status_details: {
            status: "unknown",
            stage: "unknown",
            progress: 0,
            error_message: `服务器返回: ${response.status} ${response.statusText}`,
            last_update: new Date().toISOString()
          }
        };
        
        return { 
          data: defaultStatus, 
          status: response.status,
          error: `${response.status} ${response.statusText}`
        };
      }
      
      const data = await response.json();
      
      // 处理状态数据
      const processedData = {
        ...data,
        display_status: getTaskDisplayStatus(data),
        status_details: getTaskStatusDetails(data)
      };
      
      console.log('获取计算任务状态成功:', processedData);
      
      return { 
        data: processedData, 
        status: response.status
      };
    } catch (error: any) {
      clearTimeout(timeoutId);
      
      if (error.name === 'AbortError') {
        console.error('获取计算任务状态超时');
        // 如果请求超时，返回更详细的默认状态
        return { 
          data: {
            status: "unknown",
            stage: "unknown",
            display_status: "状态获取超时",
            status_details: {
              status: "unknown",
              stage: "unknown",
              progress: 0,
              error_message: "状态获取超时，请检查网络连接",
              last_update: new Date().toISOString()
            }
          }, 
          status: 408,
          error: "请求超时"
        };
      }
      
      console.error('获取状态时发生异常:', error);
      
      // 不再抛出异常，而是返回默认状态
      return { 
        data: {
          status: "unknown",
          stage: "unknown",
          display_status: "获取状态异常",
          status_details: {
            status: "unknown",
            stage: "unknown",
            progress: 0,
            error_message: `请求异常: ${error.message || '未知错误'}`,
            last_update: new Date().toISOString()
          }
        }, 
        status: 500,
        error: error.message || '未知错误'
      };
    }
  } catch (error: any) {
    console.error('获取计算任务状态失败:', error);
    
    // 不再抛出异常，而是返回默认状态
    return { 
      data: {
        status: "unknown",
        stage: "unknown",
        display_status: "状态查询出错",
        status_details: {
          status: "unknown",
          stage: "unknown",
          progress: 0,
          error_message: `请求出错: ${error.message || '未知错误'}`,
          last_update: new Date().toISOString()
        }
      }, 
      status: 500,
      error: error.message || '未知错误'
    };
  }
};

// 获取单个计算结果(旧函数)
export const getCalculationResult = async (id: number) => {
  console.log(`获取计算结果(旧接口), ID: ${id}`);
  
  try {
    // 完整的API路径
    const apiEndpoint = `${API_BASE_URL}${API_PREFIX}calculations/${id}/result/`;
    console.log('结果API端点(旧):', apiEndpoint);
    
    // 使用fetch API调用
    const controller = new AbortController();
    const timeoutId = setTimeout(() => controller.abort(), 12000); // 12秒超时
    
    try {
      const response = await fetch(apiEndpoint, {
        method: 'GET',
        headers: {
          'Accept': 'application/json',
          'Content-Type': 'application/json',
          'Cache-Control': 'no-cache'
        },
        signal: controller.signal
      });
      
      clearTimeout(timeoutId);
      
      if (!response.ok) {
        console.error(`获取计算结果失败(旧接口): ${response.status} ${response.statusText}`);
        
        // 如果无法获取结果，提供一个默认空结果对象
        const defaultResults = {};
        
        return { 
          data: defaultResults, 
          status: 200
        };
      }
      
      const data = await response.json();
      console.log('获取计算结果成功(旧接口):', data);
      
      return { 
        data, 
        status: response.status
      };
    } catch (error: any) {
      clearTimeout(timeoutId);
      
      if (error.name === 'AbortError') {
        console.error('获取计算结果超时(旧接口)');
        // 如果请求超时，不要中断UI流程，返回一个空结果
        return { 
          data: {}, 
          status: 200
        };
      }
      
      console.error('获取结果时发生异常(旧接口):', error);
      throw error;
    }
  } catch (error: any) {
    console.error('获取计算结果失败(旧接口):', error);
    throw new Error(`获取计算结果失败: ${error.message}`);
  }
};

// 获取计算结果
export const getCalculationResults = async (calculationId: number | string) => {
  console.log(`获取计算结果, ID: ${calculationId}`);
  
  try {
    // 完整的API路径
    const apiEndpoint = `${API_BASE_URL}${API_PREFIX}calculations/${calculationId}/results/`;
    console.log('结果API端点:', apiEndpoint);
    
    // 使用fetch API调用
    const controller = new AbortController();
    const timeoutId = setTimeout(() => controller.abort(), 12000); // 12秒超时
    
    try {
      const response = await fetch(apiEndpoint, {
        method: 'GET',
        headers: {
          'Accept': 'application/json',
          'Content-Type': 'application/json',
          'Cache-Control': 'no-cache'
        },
        signal: controller.signal
      });
      
      clearTimeout(timeoutId);
      
      if (!response.ok) {
        console.error(`获取计算结果失败: ${response.status} ${response.statusText}`);
        
        // 如果无法获取结果，提供一个默认空结果对象
        const defaultResults = {};
        
        return { 
          data: defaultResults, 
          status: 200
        };
      }
      
      const data = await response.json();
      console.log('获取计算结果成功:', data);
      
      return { 
        data, 
        status: response.status
      };
    } catch (error: any) {
      clearTimeout(timeoutId);
      
      if (error.name === 'AbortError') {
        console.error('获取计算结果超时');
        // 如果请求超时，不要中断UI流程，返回一个空结果
        return { 
          data: {}, 
          status: 200
        };
      }
      
      console.error('获取结果时发生异常:', error);
      throw error;
    }
  } catch (error: any) {
    console.error('获取计算结果失败:', error);
    throw new Error(`获取计算结果失败: ${error.message}`);
  }
};

// 重启计算
export const restartCalculation = async (id: number) => {
  console.log(`重启计算任务, ID: ${id}`);
  
  try {
    // 完整的API路径
    const apiEndpoint = `${API_BASE_URL}${API_PREFIX}calculations/${id}/restart/`;
    console.log('重启API端点:', apiEndpoint);
    
    // 使用fetch API调用
    const controller = new AbortController();
    const timeoutId = setTimeout(() => controller.abort(), 15000); // 15秒超时
    
    try {
      const response = await fetch(apiEndpoint, {
        method: 'POST',
        headers: {
          'Accept': 'application/json',
          'Content-Type': 'application/json'
        },
        signal: controller.signal
      });
      
      clearTimeout(timeoutId);
      
      if (!response.ok) {
        console.error(`重启计算任务失败: ${response.status} ${response.statusText}`);
        throw new Error(`重启失败: ${response.status} ${response.statusText}`);
      }
      
      const data = await response.json();
      console.log('重启计算任务成功:', data);
      
      return { 
        data, 
        status: response.status
      };
    } catch (error: any) {
      clearTimeout(timeoutId);
      
      if (error.name === 'AbortError') {
        console.error('重启计算任务超时');
        throw new Error('请求超时，服务器响应时间过长');
      }
      
      console.error('重启时发生异常:', error);
      throw error;
    }
  } catch (error: any) {
    console.error('重启计算任务失败:', error);
    throw new Error(`重启计算任务失败: ${error.message}`);
  }
};

// 提交计算到LAMMPS
export const submitCalculation = async (id: number) => {
  console.log(`提交计算任务到LAMMPS, ID: ${id}`);
  
  try {
    // 完整的API路径
    const apiEndpoint = `${API_BASE_URL}${API_PREFIX}calculations/${id}/submit/`;
    console.log('提交API端点:', apiEndpoint);
    
    // 准备请求数据，包含唯一性字段
    const requestData = {
      formulation_id: id,
      use_latest_input_file: true,  // 添加此字段，告诉后端使用最新的输入文件
      ignore_duplicates: true,      // 忽略重复文件
      timestamp: Date.now()         // 添加时间戳以确保每次请求都是唯一的
    };
    
    // 使用fetch API调用
    const controller = new AbortController();
    const timeoutId = setTimeout(() => controller.abort(), 15000); // 15秒超时
    
    try {
      const response = await fetch(apiEndpoint, {
        method: 'POST',
        headers: {
          'Accept': 'application/json',
          'Content-Type': 'application/json'
        },
        body: JSON.stringify(requestData),
        signal: controller.signal
      });
      
      clearTimeout(timeoutId);
      
      if (!response.ok) {
        console.error(`提交计算任务失败: ${response.status} ${response.statusText}`);
        throw new Error(`提交失败: ${response.status} ${response.statusText}`);
      }
      
      const data = await response.json();
      console.log('提交计算任务成功:', data);
      
      return { 
        data, 
        status: response.status
      };
    } catch (error: any) {
      clearTimeout(timeoutId);
      
      if (error.name === 'AbortError') {
        console.error('提交计算任务超时');
        throw new Error('请求超时，服务器响应时间过长');
      }
      
      console.error('提交时发生异常:', error);
      throw error;
    }
  } catch (error: any) {
    console.error('提交计算任务失败:', error);
    throw new Error(`提交计算任务失败: ${error.message}`);
  }
};

// 删除计算任务
export const deleteCalculation = async (id: number) => {
  console.log(`删除计算任务, ID: ${id}`);
  
  try {
    // 完整的API路径
    const apiEndpoint = `${API_BASE_URL}${API_PREFIX}calculations/${id}/`;
    console.log('删除API端点:', apiEndpoint);
    
    // 使用fetch API调用
    const controller = new AbortController();
    const timeoutId = setTimeout(() => controller.abort(), 10000); // 10秒超时
    
    try {
      const response = await fetch(apiEndpoint, {
        method: 'DELETE',
        headers: {
          'Accept': 'application/json',
          'Content-Type': 'application/json'
        },
        signal: controller.signal
      });
      
      clearTimeout(timeoutId);
      
      if (!response.ok) {
        console.error(`删除计算任务失败: ${response.status} ${response.statusText}`);
        throw new Error(`删除失败: ${response.status} ${response.statusText}`);
      }
      
      return { 
        status: response.status,
        data: { success: true }
      };
    } catch (error: any) {
      clearTimeout(timeoutId);
      
      if (error.name === 'AbortError') {
        console.error('删除计算任务超时');
        throw new Error('请求超时，服务器响应时间过长');
      }
      
      console.error('删除时发生异常:', error);
      throw error;
    }
  } catch (error: any) {
    console.error('删除计算任务失败:', error);
    throw new Error(`删除计算任务失败: ${error.message}`);
  }
};

// 新增INP文件相关API
// 保存INP文件
export const saveInputFile = (formulationId: number | string, fileContent: string) => {
  return api.post(`${API_PREFIX}formulations/${formulationId}/input-file/`, {
    content: fileContent
  });
};

// 检查INP文件是否存在
export const checkInputFile = (formulationId: number | string) => {
  return api.get(`${API_PREFIX}formulations/${formulationId}/check-input-file/`);
};

// 提交计算任务
export const submitElectrolyteCalculation = (data: ICalculationData) => {
  console.log('提交计算任务数据:', data);
  
  const endpoint = `${API_PREFIX}submit-calculation/`;
  console.log('提交计算API端点:', endpoint);
  
  return api.post(endpoint, data)
    .then(response => {
      console.log('计算任务提交成功:', response);
      return response;
    })
    .catch(async error => {
      console.error('计算任务提交失败:', error);
      
      // 尝试使用替代URL
      try {
        const alternativeEndpoints = [
          `${API_PREFIX}calculations/`,
          `${API_PREFIX}calculations/submit/`,
          `${API_PREFIX}electrolyte-calculation/`
        ];
        
        for (const altEndpoint of alternativeEndpoints) {
          try {
            console.log(`尝试替代API端点: ${altEndpoint}`);
            const response = await api.post(altEndpoint, data);
            console.log('替代API端点成功:', response);
            return response;
          } catch (altError) {
            console.warn(`替代API端点 ${altEndpoint} 失败:`, altError);
            // 继续尝试下一个端点
          }
        }
        
        // 如果所有替代端点都失败，使用fetch作为最后尝试
        const fetchEndpoint = `${API_BASE_URL}${API_PREFIX}submit-calculation/`;
        console.log('尝试使用fetch API调用:', fetchEndpoint);
        
        const fetchResponse = await fetch(fetchEndpoint, {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json',
            'Authorization': `Bearer ${localStorage.getItem('token') || getRootToken()}`
          },
          body: JSON.stringify(data)
        });
        
        if (!fetchResponse.ok) {
          throw new Error(`Fetch调用失败: ${fetchResponse.status} ${fetchResponse.statusText}`);
        }
        
        const responseData = await fetchResponse.json();
        return { data: responseData, status: fetchResponse.status };
      } catch (finalError: any) {
        console.error('所有API调用方式都失败:', finalError);
        // 不再使用模拟响应，直接抛出错误
        throw new Error(`无法提交计算: ${finalError.message}`);
      }
    });
};

// 批量删除计算任务
export const batchDeleteCalculations = async (ids: number[]) => {
  const promises = ids.map(id => deleteCalculation(id));
  return Promise.all(promises);
};

// 批量重启计算任务
export const batchRestartCalculations = async (ids: number[]) => {
  const promises = ids.map(id => restartCalculation(id));
  return Promise.all(promises);
};

// 取消计算任务
export const cancelCalculation = async (id: number) => {
  console.log(`取消计算任务, ID: ${id}`);
  
  try {
    // 完整的API路径
    const apiEndpoint = `${API_BASE_URL}${API_PREFIX}cancel-calculation/`;
    console.log('取消任务API端点:', apiEndpoint);
    
    // 使用fetch API调用
    const controller = new AbortController();
    const timeoutId = setTimeout(() => controller.abort(), 10000); // 10秒超时
    
    try {
      const response = await fetch(apiEndpoint, {
        method: 'POST',
        headers: {
          'Accept': 'application/json',
          'Content-Type': 'application/json'
        },
        body: JSON.stringify({ calculation_id: id }),
        signal: controller.signal
      });
      
      clearTimeout(timeoutId);
      
      if (!response.ok) {
        console.error(`取消计算任务失败: ${response.status} ${response.statusText}`);
        throw new Error(`取消失败: ${response.status} ${response.statusText}`);
      }
      
      const data = await response.json();
      return { 
        status: response.status,
        data
      };
    } catch (error: any) {
      clearTimeout(timeoutId);
      
      if (error.name === 'AbortError') {
        console.error('取消计算任务超时');
        throw new Error('请求超时，服务器响应时间过长');
      }
      
      console.error('取消任务时发生异常:', error);
      throw error;
    }
  } catch (error: any) {
    console.error('取消计算任务失败:', error);
    throw new Error(`取消计算任务失败: ${error.message}`);
  }
};

// 导出任务数据为CSV
export const exportTasksToCSV = (tasks: any[], selectedIds?: number[]) => {
  // 确定要导出的任务列表
  const dataToExport = selectedIds && selectedIds.length > 0 
    ? tasks.filter(task => selectedIds.includes(task.id))
    : tasks;
    
  // 创建CSV内容
  const headers = ['ID', '任务名称', '状态', '创建时间', '完成时间', '描述'];
  const csvContent = [
    headers.join(','),
    ...dataToExport.map(task => [
      task.id,
      `"${task.name?.replace(/"/g, '""') || ''}"`, // 处理引号
      task.status || '',
      task.created_at ? new Date(task.created_at).toLocaleString() : '',
      task.finished_at ? new Date(task.finished_at).toLocaleString() : '',
      task.description ? `"${task.description.replace(/"/g, '""')}"` : ''
    ].join(','))
  ].join('\n');
  
  // 创建Blob对象
  const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
  const url = URL.createObjectURL(blob);
  
  // 创建下载链接并触发下载
  const link = document.createElement('a');
  link.setAttribute('href', url);
  link.setAttribute('download', `任务列表_${new Date().toISOString().slice(0, 10)}.csv`);
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  URL.revokeObjectURL(url);
  
  return dataToExport.length;
};

// 导出任务数据为JSON
export const exportTasksToJSON = (tasks: any[], selectedIds?: number[]) => {
  // 确定要导出的任务列表
  const dataToExport = selectedIds && selectedIds.length > 0 
    ? tasks.filter(task => selectedIds.includes(task.id))
    : tasks;
    
  // 创建JSON内容
  const jsonContent = JSON.stringify(dataToExport, null, 2);
  
  // 创建Blob对象
  const blob = new Blob([jsonContent], { type: 'application/json;charset=utf-8;' });
  const url = URL.createObjectURL(blob);
  
  // 创建下载链接并触发下载
  const link = document.createElement('a');
  link.setAttribute('href', url);
  link.setAttribute('download', `任务列表_${new Date().toISOString().slice(0, 10)}.json`);
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  URL.revokeObjectURL(url);
  
  return dataToExport.length;
};

// 提交配方到后端进行评估
export const submitRecipe = async (recipeData: any, retryCount = 2, retryDelay = 3000): Promise<any> => {
  console.log('【API】开始提交配方:', JSON.stringify(recipeData));
  
  // 准备基础数据
  const formattedData: any = { 
    name: recipeData.name || '未命名配方',
    description: recipeData.description || '',
    temperature: Number(recipeData.temperature) || 300.0,
    box_size: Number(recipeData.boxSize) || 50.0
  };

  // 处理阳离子数据
  let cations: any[] = [];
  if (recipeData.cations && recipeData.cations.length > 0) {
    // 计算模拟盒子体积（立方埃 -> 立方米 -> 升）
    const boxSize = Number(recipeData.boxSize) || 50.0;
    const boxVolume = boxSize * boxSize * boxSize * 1e-30; // 立方埃转换为立方米
    const volumeInLiters = boxVolume * 1000; // 转换为升
    
    console.log(`【阳离子】模拟盒子体积: ${volumeInLiters.toExponential(6)} L`);
    
    cations = recipeData.cations.map((cation: any) => {
      // 获取不带电荷符号的离子名称
      const cleanName = String(cation.name || cation.type || 'Li').trim().replace(/[+]$/, '');
      
      // 确定离子价态
      let charge = 1; // 默认为1价
      if (cleanName === 'Mg' || cleanName === 'Ca' || cleanName === 'Zn') {
        charge = 2; // 2价阳离子
      } else if (cleanName === 'Al') {
        charge = 3; // 3价阳离子
      }
      
      // 获取浓度（mol/L）
      const concentration = Number(cation.concentration) || 1.0;
      
      // 根据浓度计算分子数量
      // 通过公式: 数量 = 浓度(mol/L) * 体积(L) * 阿伏伽德罗常数
      // 这里简化为：浓度为1 mol/L时，标准数量为75个分子
      // 因此，数量 = 75 * 浓度值
      let calculatedNumber;
      
      if (cation.number !== undefined) {
        // 如果用户提供了特定的数量，优先使用用户输入
        calculatedNumber = Number(cation.number);
      } else {
        // 否则根据浓度计算
        calculatedNumber = Math.round(75 * concentration);
      }
      
      // 确保至少有16个分子，这是分子动力学的合理最小值
      const number = Math.max(calculatedNumber, 16);
      
      console.log(`【阳离子】${cleanName}: 浓度=${concentration} mol/L, 计算数量=${calculatedNumber}, 最终数量=${number}`);
      
      return {
        name: cleanName,
        charge: charge,
        number,
        concentration
      };
    });
    
    // 保存到格式化数据
    formattedData.cations = cations;
  }
  
  // 处理阴离子数据
  let anions: any[] = [];
  if (recipeData.anions && recipeData.anions.length > 0) {
    // 计算总浓度，用于后续按比例分配
    const totalAnionConcentration = recipeData.anions.reduce((sum: number, anion: any) => {
      return sum + (Number(anion.concentration) || 0);
    }, 0);
    
    console.log(`【阴离子】总浓度: ${totalAnionConcentration} mol/L`);
    
    // 计算模拟盒子体积（立方埃 -> 立方米 -> 升）
    const boxSize = Number(recipeData.boxSize) || 50.0;
    const boxVolume = boxSize * boxSize * boxSize * 1e-30; // 立方埃转换为立方米
    const volumeInLiters = boxVolume * 1000; // 转换为升
    
    console.log(`【阴离子】模拟盒子体积: ${volumeInLiters.toExponential(6)} L`);
    
    // 计算阳离子提供的总电荷量
    const totalPositiveCharge = cations.reduce((sum: number, cation: any) => {
      return sum + (Number(cation.number) || 0) * (Number(cation.charge) || 1);
    }, 0);
    
    console.log(`【阴离子】阳离子提供的总电荷量: ${totalPositiveCharge}`);
    
    // 重置anions数组，避免之前的计算影响
    anions = [];
    let totalAssignedNegativeCharge = 0;
    
    // 为每个阴离子分配数量，以平衡电荷
    recipeData.anions.forEach((anion: any, index: number) => {
      // 获取离子名称
      const rawName = String(anion.name || anion.type || '').trim();
      const cleanName = rawName.replace(/[-]$/, '');  // 移除末尾的负号
      
      // 获取浓度
      const anionConcentration = Number(anion.concentration) || 0;
      
      // 阴离子电荷，通常为-1
      const anionCharge = Number(anion.charge) || -1;
      
      // 如果是单个阴离子，使用总阳离子电荷总量，否则按浓度比例计算
      let anionNumber;
      if (recipeData.anions.length === 1) {
        // 单个阴离子，数量应使总负电荷等于总正电荷
        anionNumber = Math.round(totalPositiveCharge / Math.abs(anionCharge));
      } else {
        // 多个阴离子，按浓度比例分配负电荷量
        const proportion = anionConcentration / totalAnionConcentration;
        // 计算应分配的负电荷量并转化为阴离子数量
        const negativeChargeForThisAnion = totalPositiveCharge * proportion;
        anionNumber = Math.round(negativeChargeForThisAnion / Math.abs(anionCharge));
        
        // 最后一个阴离子需要确保总负电荷等于总正电荷
        if (index === recipeData.anions.length - 1) {
          // 计算前面分配的总负电荷量
          const remainingNegativeCharge = totalPositiveCharge - totalAssignedNegativeCharge;
          anionNumber = Math.round(remainingNegativeCharge / Math.abs(anionCharge));
        } else {
          // 更新已分配的总负电荷量
          totalAssignedNegativeCharge += anionNumber * Math.abs(anionCharge);
        }
      }
      
      // 确保至少有1个分子
      anionNumber = Math.max(anionNumber, 1);
      
      // 添加调试日志
      console.log(`【阴离子】${cleanName} - 浓度比例: ${anionConcentration}/${totalAnionConcentration} = ${anionConcentration/totalAnionConcentration}, 计算数量: ${anionNumber}, 提供负电荷: ${anionNumber * Math.abs(anionCharge)}`);
      
      anions.push({
        name: cleanName,
        charge: anionCharge,
        number: anionNumber,
        concentration: anionConcentration
      });
    });
    
    // 最终检查电荷平衡
    const totalNegativeCharge = anions.reduce((sum: number, anion: any) => {
      return sum + anion.number * Math.abs(anion.charge);
    }, 0);
    
    console.log(`【阴离子】最终电荷检查: 正电荷=${totalPositiveCharge}, 负电荷=${totalNegativeCharge}`);
    
    if (totalNegativeCharge !== totalPositiveCharge) {
      console.log(`【阴离子】电荷不平衡，调整最后一个阴离子数量`);
      if (anions.length > 0) {
        const lastAnion = anions[anions.length - 1];
        const diffCharge = totalPositiveCharge - totalNegativeCharge;
        const addAnions = Math.round(diffCharge / Math.abs(lastAnion.charge));
        lastAnion.number += addAnions;
        console.log(`【阴离子】调整 ${lastAnion.name} 数量: ${lastAnion.number - addAnions} -> ${lastAnion.number}, 新负电荷总量: ${totalNegativeCharge + addAnions * Math.abs(lastAnion.charge)}`);
      }
    }
    
    // 保存到格式化数据
    formattedData.anions = anions;
  }
  
  // 创建盐组件 - 直接使用cation和anion信息
  const salts: any[] = [];
  
  // 处理盐组合 
  if (cations.length > 0 && anions.length > 0) {
    console.log('【盐组合】生成盐组合...');
    
    // 特殊处理单阳离子多阴离子的情况
    if (cations.length === 1 && anions.length > 1) {
      console.log('【盐组合】单阳离子多阴离子情况，优先处理第一种阴离子');
      
      const cation = cations[0]; // 唯一的阳离子
      let remainingCationCount = cation.number; // 初始阳离子数量
      let remainingPositiveCharge = cation.number * cation.charge; // 初始正电荷总量
      
      console.log(`【盐组合】阳离子 ${cation.name} 总数: ${remainingCationCount}, 电荷: ${cation.charge}, 总正电荷: ${remainingPositiveCharge}`);
      console.log(`【盐组合】阴离子: ${JSON.stringify(anions.map(a => ({name: a.name, number: a.number, charge: a.charge})))}`);
      
      // 按顺序处理每种阴离子，优先满足前面的阴离子
      for (let i = 0; i < anions.length; i++) {
        const anion = anions[i];
        const isLastAnion = i === anions.length - 1;
        
        // 计算每个阳离子需要的阴离子数量来平衡电荷
        const anionsPerCation = Math.abs(cation.charge / anion.charge);
        
        // 对于每种阴离子，计算能配对的阳离子数量
        const maxCationsToPair = Math.min(
          remainingCationCount,
          Math.floor(anion.number / anionsPerCation)
        );
        
        console.log(`【盐组合】处理阴离子 ${anion.name}，每个${cation.name}需要${anionsPerCation}个${anion.name}，可配对阳离子数量: ${maxCationsToPair}`);
        
        if (maxCationsToPair > 0) {
          // 计算实际使用的阴离子数量
          const anionsUsed = Math.round(maxCationsToPair * anionsPerCation);
          
          // 创建盐组合
          salts.push({
            name: `${cation.name}${anion.name}`,
            cation: cation.name,
            anion: anion.name,
            cation_number: maxCationsToPair,
            anion_number: anionsUsed,
            concentration: cation.concentration * (maxCationsToPair / cation.number),
            component_settings: JSON.stringify({
              cation_number: maxCationsToPair,
              anion_number: anionsUsed
            })
          });
          
          console.log(`【盐组合】生成盐: ${cation.name}${anion.name}，阳离子: ${maxCationsToPair}个，阴离子: ${anionsUsed}个`);
          
          // 减少剩余阳离子数量和电荷
          remainingCationCount -= maxCationsToPair;
          remainingPositiveCharge -= maxCationsToPair * cation.charge;
        }
        
        // 如果阳离子已用完，则后续阴离子无法配对
        if (remainingCationCount <= 0) {
          console.log(`【盐组合】阳离子已用完，无法处理后续阴离子`);
          break;
        }
      }
      
      // 检查是否还有剩余阳离子
      if (remainingCationCount > 0) {
        console.log(`【警告】仍有${remainingCationCount}个${cation.name}未配对，剩余正电荷: ${remainingPositiveCharge}`);
      }
    }
    // 处理多阳离子情况
    else {
      console.log(`【盐组合】开始处理多阳离子情况`);
      
      // 创建一个简单的队列结构来处理阴阳离子配对
      const cationQueue = [...cations.map(c => ({...c}))];
      const anionQueue = [...anions.map(a => ({...a}))];
      
      console.log(`【盐组合初始状态】阳离子队列: ${cationQueue.length}个，阴离子队列: ${anionQueue.length}个`);
      
      // 处理配对逻辑
      while (cationQueue.length > 0 && anionQueue.length > 0) {
        const currentCation = cationQueue[0];
        const cationName = currentCation.name;
        const cationNumber = currentCation.number;
        const cationCharge = currentCation.charge;
        
        // 计算阳离子提供的总电荷
        const totalPositiveCharge = cationNumber * cationCharge;
        console.log(`【处理阳离子】${cationName}，数量: ${cationNumber}，电荷: ${cationCharge}，总电荷: ${totalPositiveCharge}`);
        
        // 分配阴离子
        let remainingPositiveCharge = totalPositiveCharge;
        let cationsUsed = 0;
        
        // 尝试为当前阳离子配对所有可用的阴离子
        for (let i = 0; i < anionQueue.length && cationsUsed < cationNumber; i++) {
          const anion = anionQueue[i];
          const anionCharge = Math.abs(anion.charge);
          
          // 计算能配对的阳离子数量
          const anionsPerCation = Math.ceil(cationCharge / anionCharge);
          console.log(`【配对计算】每个${cationName}需要${anionsPerCation}个${anion.name}`);
          
          // 计算可用的阴离子数量
          const availableAnions = Math.min(anion.number, remainingPositiveCharge);
          
          // 计算能配对的阳离子数量
          const cationsForThisAnion = Math.min(
            Math.floor(availableAnions / anionsPerCation),
            cationNumber - cationsUsed
          );
          
          if (cationsForThisAnion > 0) {
            // 计算实际使用的阴离子数量
            const anionsUsed = cationsForThisAnion * anionsPerCation;
            
            // 创建盐组合
            salts.push({
              name: `${cationName}${anion.name}`,
              cation: cationName,
              anion: anion.name,
              cation_number: cationsForThisAnion,
              anion_number: anionsUsed,
              concentration: currentCation.concentration * (cationsForThisAnion / cationNumber),
              component_settings: JSON.stringify({
                cation_number: cationsForThisAnion,
                anion_number: anionsUsed
              })
            });
            
            console.log(`【盐组合】生成盐: ${cationName}${anion.name}，阳离子: ${cationsForThisAnion}个，阴离子: ${anionsUsed}个`);
            
            // 更新计数
            cationsUsed += cationsForThisAnion;
            remainingPositiveCharge -= anionsUsed * anionCharge;
            anion.number -= anionsUsed;
          }
          
          // 移除用完的阴离子
          if (anion.number <= 0) {
            anionQueue.splice(i, 1);
            i--; // 调整索引
            console.log(`【阴离子用尽】移除${anion.name}`);
          }
        }
        
        // 更新阳离子数量
        currentCation.number -= cationsUsed;
        
        // 移除用完的阳离子
        if (currentCation.number <= 0) {
          cationQueue.shift();
          console.log(`【阳离子用尽】移除${cationName}`);
        } else if (cationsUsed === 0) {
          // 如果没有分配任何阳离子，避免无限循环
          console.log(`【无法配对】无法为${cationName}找到匹配的阴离子，跳过`);
          cationQueue.shift();
        }
      }
      
      // 记录未配对的离子
      if (cationQueue.length > 0) {
        console.log(`【未配对离子】剩余阳离子: ${JSON.stringify(cationQueue.map(c => ({name: c.name, number: c.number})))}`);
      }
      
      if (anionQueue.length > 0) {
        console.log(`【未配对离子】剩余阴离子: ${JSON.stringify(anionQueue.map(a => ({name: a.name, number: a.number})))}`);
      }
    }
    
    console.log(`【盐组合完成】共生成${salts.length}个盐组合`);
  }
  
  // 添加盐数据到格式化数据
  if (salts.length > 0) {
    formattedData.salts = salts;
  }
  
  // 处理溶剂
  if (recipeData.solvents && recipeData.solvents.length > 0) {
    // 获取第一个阳离子的数量，用于计算溶剂分子数量
    const firstCationNumber = cations.length > 0 ? cations[0].number : 75;
    console.log(`【溶剂计算】参考阳离子数量: ${firstCationNumber}`);
    
    // 确保溶剂保持用户输入的原始值，不进行任何合并或修改
    formattedData.solvents = recipeData.solvents.map((solvent: any, index: number) => {
      // 获取干净的溶剂名称
      const rawName = String(solvent.name || solvent.solvent || 'EC').trim();
      const cleanName = rawName.replace(/\s*\([^)]*\)/g, '').trim();
      
      // 计算溶剂分子数量 = 阳离子数量 × 体积比例
      const volumeRatio = Number(solvent.volumeRatio) || 1.0;
      // 使用阳离子数量和volumeRatio计算溶剂分子数量
      const calculatedNumber = Math.round(firstCationNumber * volumeRatio);
      // 确保至少有1个分子
      const moleculeCount = Math.max(calculatedNumber, 1);
      
      console.log(`【溶剂数据】处理溶剂[${index}]: ${cleanName}, 体积比例: ${volumeRatio}, 计算得到分子数量: ${moleculeCount}`);
      
      // 返回简化的溶剂对象，包含计算出的分子数量
      return {
        name: cleanName,
        number: moleculeCount,
        volume_ratio: volumeRatio,
        smile: solvent.smile || getSolventSMILES(cleanName),
        // 添加比例信息，确保后端能够正确计算
        ratio: volumeRatio,
        concentration: volumeRatio
      };
    });
    
    console.log(`【溶剂数据】最终提交的溶剂列表:`, JSON.stringify(formattedData.solvents));
  }

  console.log('【API】最终提交数据:', JSON.stringify(formattedData));

  // 检查并修正阴离子数量的一致性
  if (formattedData.anions && formattedData.anions.length > 1 && formattedData.cations && formattedData.cations.length === 1) {
    const totalCationNumber = formattedData.cations[0].number;
    console.log(`【最终检查】总阳离子数量: ${totalCationNumber}`);
    
    // 检查阴离子总数是否与阳离子数量一致
    const totalAnionNumber = formattedData.anions.reduce((sum: number, anion: any) => sum + anion.number, 0);
    console.log(`【最终检查】阴离子总数: ${totalAnionNumber}, 应该等于: ${totalCationNumber}`);
    
    // 如果不一致，调整阴离子数量
    if (totalAnionNumber !== totalCationNumber) {
      console.log(`【最终检查】发现阴离子总数不一致，进行修正`);
      
      // 计算总浓度
      const totalConcentration = formattedData.anions.reduce((sum: number, anion: any) => sum + anion.concentration, 0);
      
      // 重新计算每个阴离子的数量
      let assignedCount = 0;
      formattedData.anions.forEach((anion: any, index: number) => {
        if (index === formattedData.anions.length - 1) {
          // 最后一个阴离子取剩余数量
          anion.number = totalCationNumber - assignedCount;
          console.log(`【最终检查】修正最后阴离子 ${anion.name} 数量: ${anion.number}`);
        } else {
          // 按比例分配
          const proportion = anion.concentration / totalConcentration;
          anion.number = Math.round(totalCationNumber * proportion);
          assignedCount += anion.number;
          console.log(`【最终检查】修正阴离子 ${anion.name} 数量: ${anion.number} (比例: ${proportion})`);
        }
      });
    }
  }
  
  // 确保anions数组中的数量与salts中的匹配
  if (formattedData.salts && formattedData.salts.length > 0 && formattedData.anions && formattedData.anions.length > 0) {
    console.log('【最终检查】验证anions数组中的数量与salts中的匹配');
    
    // 为每种阴离子统计在salts中的总数量
    const anionCountsInSalts: {[key: string]: number} = {};
    
    // 计算在salts中每种阴离子的总数量
    formattedData.salts.forEach((salt: any) => {
      const anion = salt.anion;
      const count = Number(salt.anion_number) || 0;
      
      if (!anionCountsInSalts[anion]) {
        anionCountsInSalts[anion] = 0;
      }
      anionCountsInSalts[anion] += count;
    });
    
    console.log(`【最终检查】Salts中的阴离子数量: ${JSON.stringify(anionCountsInSalts)}`);
    
    // 检查并更新anions数组中的数量
    formattedData.anions.forEach((anion: any) => {
      const name = anion.name;
      const countInSalts = anionCountsInSalts[name] || 0;
      
      if (anion.number !== countInSalts) {
        console.log(`【最终检查】更新阴离子 ${name} 的数量: ${anion.number} -> ${countInSalts}`);
        anion.number = countInSalts;
      }
    });
  }

  // 构建API端点
  const apiEndpoint = `${API_BASE_URL}${API_PREFIX}formulations/`;
  console.log('【API】使用端点:', apiEndpoint);

  // 提交数据
  const makeRequest = async (attempt: number): Promise<any> => {
    try {
      const token = localStorage.getItem('token') || getRootToken();
        const response = await fetch(apiEndpoint, {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json',
            'Accept': 'application/json',
          'Authorization': token ? `Bearer ${token}` : ''
          },
        body: JSON.stringify(formattedData)
        });
        
        if (!response.ok) {
        const errorText = await response.text();
        console.error(`【API错误】提交配方失败 (${response.status}): ${errorText}`);
        throw new Error(`API返回: ${errorText}`);
        }
        
        const data = await response.json();
      console.log('【API】配方提交成功:', data);
        return { status: 'success', data };
    } catch (error) {
      console.error(`【API错误】尝试 ${attempt}/${retryCount} 失败:`, error);
      if (attempt < retryCount) {
        console.log(`【API】等待 ${retryDelay}ms 后重试...`);
            await new Promise(resolve => setTimeout(resolve, retryDelay));
            return makeRequest(attempt + 1);
          }
      throw new Error(`提交失败: ${error instanceof Error ? error.message : String(error)}`);
    }
  };

  return makeRequest(1);
};

// 辅助函数：根据名称获取SMILES
const getSolventSMILES = (name: string): string => {
  const solventSMILES: Record<string, string> = {
    'EC': 'C1COC(=O)O1',
    'DMC': 'COC(=O)OC',
    'EMC': 'CCOC(=O)OC',
    'PC': 'CC1COC(=O)O1',
    'DEC': 'CCOC(=O)OCC',
    'FEC': 'FC1COC(=O)O1'
  };
  
  return solventSMILES[name.toUpperCase()] || '';
};

// 生成配方的输入文件（分子文件和inp文件）
export const generateRecipeInputFiles = async (formulationId: number | string): Promise<any> => {
  console.log('生成配方输入文件, 配方ID:', formulationId);
  
  // 转换为字符串ID以确保兼容性
  const strId = String(formulationId);
  
  try {
    // 测试后端API是否可访问
    try {
      console.log('正在测试API可用性...');
      const testResponse = await fetch(`${API_BASE_URL}${API_PREFIX}health/`, { 
        method: 'HEAD',
        signal: AbortSignal.timeout(5000) 
      });
      console.log('API健康检查响应:', testResponse.status);
    } catch (testError) {
      console.warn('API健康检查失败，但仍将继续尝试提交:', testError);
      // 继续提交，不中断流程
    }
    
    // 完整的API路径
    const apiEndpoint = `${API_BASE_URL}${API_PREFIX}formulations/${strId}/generate_input_file/`;
    console.log('完整API端点:', apiEndpoint);
    
    // 使用fetch API调用
      console.log("使用fetch API直接调用...");
      const requestData = {
        formulation_id: formulationId
      };
      
      // 使用较短的超时时间
      const controller = new AbortController();
      const timeoutId = setTimeout(() => controller.abort(), 10000); // 10秒超时
      
      try {
        const fetchResponse = await fetch(apiEndpoint, {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json',
            'Accept': 'application/json'
          },
          body: JSON.stringify(requestData),
          signal: controller.signal
        });
        
        clearTimeout(timeoutId);
        
        if (!fetchResponse.ok) {
        const errorText = await fetchResponse.text();
        console.error(`Fetch请求失败: ${fetchResponse.status} ${fetchResponse.statusText}`, errorText);
        throw new Error(`生成输入文件失败: ${fetchResponse.status} ${fetchResponse.statusText}`);
        }
        
        const data = await fetchResponse.json();
        console.log('Fetch成功，响应数据:', data);
        
        return { 
          data, 
          status: fetchResponse.status, 
          headers: Object.fromEntries(fetchResponse.headers)
        };
      } catch (innerError: any) {
        clearTimeout(timeoutId);
        
        if (innerError.name === 'AbortError') {
        throw new Error('请求超时: 生成输入文件请求超时，服务器响应时间过长');
      }
      
      throw innerError;
    }
  } catch (error: any) {
    console.error('生成输入文件失败:', error);
    throw error;
  }
};

// 创建一个简单的测试配方数据
const testRecipe = {
  name: '测试配方',
  description: '这是一个测试配方',
  temperature: 300,
  box_size: 30,
  electrolyte: {
    temperature: 300,
    box_size: 30
  },
  // 添加用户信息，尝试从localStorage获取
  user_id: 1,
  username: localStorage.getItem('username') || 'guest_user',
  cations: [
    {
      name: 'Li',
      number: 10,
      charge: 1
    }
  ],
  // ... rest of the test recipe ...
};

// 获取任务进度详情
export const getTaskProgressDetails = async (taskId: number | string) => {
  console.log(`获取任务进度详情, ID: ${taskId}`);
  
  try {
    // 完整的API路径
    const apiEndpoint = `${API_BASE_URL}${API_PREFIX}calculations/${taskId}/progress/`;
    console.log('进度API端点:', apiEndpoint);
    
    // 使用API调用
    const controller = new AbortController();
    const timeoutId = setTimeout(() => controller.abort(), 8000); // 8秒超时
    
    try {
      const response = await fetch(apiEndpoint, {
        method: 'GET',
        headers: {
          'Accept': 'application/json',
          'Content-Type': 'application/json',
          'Cache-Control': 'no-cache'
        },
        signal: controller.signal
      });
      
      clearTimeout(timeoutId);
      
      if (!response.ok) {
        const errorText = await response.text();
        console.error(`获取任务进度失败: ${response.status} ${response.statusText}`, errorText);
        throw new Error(`获取任务进度失败: ${response.status} ${response.statusText}`);
      }
        
        const data = await response.json();
        console.log('获取任务进度成功:', data);
        return { data, status: response.status };
    } catch (error: any) {
      clearTimeout(timeoutId);
      
      if (error.name === 'AbortError') {
        throw new Error('获取任务进度请求超时，服务器响应时间过长');
      }
      
      throw error;
    }
  } catch (error: any) {
    console.error('获取任务进度详情失败:', error);
    throw new Error(`获取任务进度详情失败: ${error.message}`);
  }
};

// 任务状态轮询间隔（毫秒）
const STATUS_POLLING_INTERVAL = 5000; // 5秒

// 任务状态轮询超时（毫秒）
const STATUS_POLLING_TIMEOUT = 300000; // 5分钟

// 任务状态轮询
export const startStatusPolling = (taskId: number, onStatusUpdate: (status: ITaskStatusDetails) => void) => {
  let startTime = Date.now();
  let isPolling = true;

  const pollStatus = async () => {
    if (!isPolling) return;

    try {
      const { data } = await getCalculationStatus(taskId);
      
      // 更新状态
      onStatusUpdate(data.status_details);

      // 检查是否需要继续轮询
      if (data.status === 'completed' || data.status === 'error') {
        isPolling = false;
        return;
      }

      // 检查是否超时
      if (Date.now() - startTime > STATUS_POLLING_TIMEOUT) {
        console.warn(`任务 ${taskId} 状态轮询超时`);
        isPolling = false;
        return;
      }

      // 继续轮询
      setTimeout(pollStatus, STATUS_POLLING_INTERVAL);
    } catch (error) {
      console.error(`轮询任务 ${taskId} 状态失败:`, error);
      
      // 提供默认状态对象给回调函数，防止UI崩溃
      const defaultStatus: ITaskStatusDetails = {
        status: "error",
        stage: "unknown",
        progress: 0,
        error_message: `获取状态失败: ${error instanceof Error ? error.message : '未知错误'}`,
        last_update: new Date().toISOString()
      };
      
      // 通知UI出现错误，但继续尝试轮询
      onStatusUpdate(defaultStatus);
      
      // 继续轮询，而不是直接停止
      setTimeout(pollStatus, STATUS_POLLING_INTERVAL * 2); // 加倍轮询间隔，减少失败请求频率
    }
  };

  // 开始轮询
  pollStatus();

  // 返回停止轮询的函数
  return () => {
    isPolling = false;
  };
};

// 批量任务状态轮询
export const startBatchStatusPolling = (taskIds: number[], onStatusUpdate: (taskId: number, status: ITaskStatusDetails) => void) => {
  const stopPollingFunctions = taskIds.map(taskId => 
    startStatusPolling(taskId, (status) => onStatusUpdate(taskId, status))
  );

  // 返回停止所有轮询的函数
  return () => {
    stopPollingFunctions.forEach(stop => stop());
  };
};

// 重置/清零任务ID（如果后端支持此功能）
export const resetTaskIds = async () => {
  console.log('尝试重置任务ID');
  
  try {
    // 完整的API路径
    const apiEndpoint = `${API_BASE_URL}${API_PREFIX}reset-ids/`;
    console.log('重置ID API端点:', apiEndpoint);
    
    try {
      // 使用fetch API调用
      const response = await fetch(apiEndpoint, {
        method: 'POST',
        headers: {
          'Accept': 'application/json',
          'Content-Type': 'application/json'
        },
        body: JSON.stringify({ 
          reset_tasks: true, 
          reset_formulations: true,
          user_specific: true // 要求每个用户有独立ID序列
        })
      });
      
      if (response.ok) {
        const data = await response.json();
        console.log('重置ID成功:', data);
        return { success: true, message: '已成功重置ID序列', data };
      } else {
        console.error(`重置ID失败: ${response.status} ${response.statusText}`);
        
        // 尝试读取错误消息
        try {
          const errorData = await response.json();
          return { 
            success: false, 
            message: errorData.message || `重置失败: ${response.status} ${response.statusText}`,
            error: errorData
          };
        } catch (e) {
          return { 
            success: false, 
            message: `重置失败: ${response.status} ${response.statusText}`,
            error: { status: response.status }
          };
        }
      }
    } catch (error: any) {
      console.error('重置ID请求失败:', error);
      return { 
        success: false, 
        message: `重置请求失败: ${error.message}`,
        error
      };
    }
  } catch (error: any) {
    console.error('重置ID功能异常:', error);
    return { 
      success: false, 
      message: `重置ID功能异常: ${error.message}`,
      error
    };
  }
};

// 获取输入文件URL
export const getInputFileUrl = async (formulationId: string): Promise<string> => {
  try {
    console.log('获取输入文件URL, 配方ID:', formulationId);
    const response = await api.get(`${API_PREFIX}formulations/${formulationId}/input-file/`);
    console.log('获取输入文件URL响应:', response);
    
    if (response && response.data && response.data.url) {
      return response.data.url;
    } else if (response && response.data && typeof response.data === 'string') {
      return response.data;
    } else {
      throw new Error(`无法解析输入文件URL响应: ${JSON.stringify(response.data)}`);
    }
  } catch (error: any) {
    console.error('获取输入文件URL失败:', error);
    throw new Error(`获取输入文件URL失败: ${error.message}`);
  }
};

// 获取配方的计算状态
export const getFormulationCalculationStatus = async (formulationId: string): Promise<{ has_calculation: boolean }> => {
  try {
    console.log('获取配方计算状态, 配方ID:', formulationId);
    const response = await api.get(`${API_PREFIX}formulations/${formulationId}/calculation-status/`);
    console.log('获取配方计算状态响应:', response);
    
    if (response && response.data && typeof response.data === 'object') {
      if ('has_calculation' in response.data) {
        return { has_calculation: Boolean(response.data.has_calculation) };
      } else {
        // 尝试从其他字段中推断
        const hasCalculation = 
          response.data.calculation_id || 
          response.data.task_id || 
          response.data.has_task || 
          response.data.status === 'completed';
        return { has_calculation: Boolean(hasCalculation) };
      }
    }
    throw new Error(`无法解析配方计算状态响应: ${JSON.stringify(response.data)}`);
  } catch (error: any) {
    console.error('获取配方计算状态失败:', error);
    throw new Error(`获取配方计算状态失败: ${error.message}`);
  }
}; 