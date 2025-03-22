import api from './api';

// 获取所有溶剂
export const getSolvents = () => {
  return api.get('/solvents/');
};

// 获取所有盐
export const getSalts = () => {
  return api.get('/salts/');
};

// 获取用户的所有配方
export const getFormulations = () => {
  return api.get('/formulations/');
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
  console.log('创建配方:', data);
  return api.post('/formulations/', data);
};

// 获取特定配方
export const getFormulation = (id: number) => {
  return api.get(`/formulations/${id}/`);
};

// 获取配方组件
export const getFormulationComponents = (id: number) => {
  return api.get(`/formulations/${id}/components/`);
};

// 获取配方参数
export const getFormulationParameters = (id: number) => {
  return api.get(`/formulations/${id}/parameters/`);
};

// 生成输入文件
export const generateInputFile = (id: number) => {
  console.log(`调用生成输入文件API: /formulations/${id}/generate_input_file/`);
  console.log(`参数类型: ${typeof id}, 参数值: ${id}`);
  
  // 确保id是数字
  const formattedId = typeof id === 'string' ? parseInt(id) : id;
  
  console.log(`格式化后ID: ${formattedId}, 类型: ${typeof formattedId}`);
  console.log(`API基础URL: ${api.defaults.baseURL}`);
  console.log(`完整请求URL: ${api.defaults.baseURL}/formulations/${formattedId}/generate_input_file/`);
  
  // 添加详细错误处理和超时设置
  return api.post(`/formulations/${formattedId}/generate_input_file/`, {}, {
    timeout: 30000, // 设置30秒超时
    headers: {
      'Content-Type': 'application/json',
      'X-Requested-With': 'XMLHttpRequest'
    }
  })
    .then(response => {
      console.log('生成输入文件API调用成功:', response);
      console.log('响应状态码:', response.status);
      console.log('响应数据:', response.data);
      return response;
    })
    .catch(error => {
      console.error('生成输入文件API调用失败:', error);
      console.error('错误类型:', error.constructor.name);
      
      if (error.response) {
        // 服务器响应了错误状态码
        console.error('错误状态码:', error.response.status);
        console.error('错误响应数据:', error.response.data);
        console.error('错误响应头:', error.response.headers);
      } else if (error.request) {
        // 请求已发送但没有收到响应
        console.error('请求已发送但没有收到响应:', error.request);
        console.error('请求配置:', error.config);
      } else {
        // 请求设置时发生错误
        console.error('请求设置错误:', error.message);
      }
      
      console.error('完整错误对象:', JSON.stringify(error, null, 2));
      throw error;
    });
};

// 获取计算任务列表
export const getCalculations = () => {
  return api.get('/calculations/');
};

// 创建新的计算
export interface ICalculationData {
  name: string;
  description?: string;
  formulation_id: number;
}

export const createCalculation = (data: ICalculationData) => {
  return api.post('/calculations/', data);
};

// 获取计算详情
export const getCalculation = (id: number) => {
  return api.get(`/calculations/${id}/`);
};

// 获取计算结果
export const getCalculationResult = (id: number) => {
  return api.get(`/calculations/${id}/result/`);
};

// 重启计算
export const restartCalculation = (id: number) => {
  return api.post(`/calculations/${id}/restart/`);
};

// 提交计算到LAMMPS
export const submitCalculation = (id: number) => {
  return api.post(`/calculations/${id}/submit/`);
};

// 新增INP文件相关API
// 保存INP文件
export const saveInputFile = (formulationId: number | string, fileContent: string) => {
  return api.post(`/formulations/${formulationId}/input-file/`, {
    content: fileContent
  });
};

// 检查INP文件是否存在
export const checkInputFile = (formulationId: number | string) => {
  return api.get(`/formulations/${formulationId}/check-input-file/`);
};

// 提交计算任务
export const submitElectrolyteCalculation = (data: ICalculationData) => {
  return api.post('/submit-calculation/', data);
};

// 获取计算状态
export const getCalculationStatus = (calculationId: number | string) => {
  return api.get(`/calculations/${calculationId}/status/`);
};

// 获取计算结果
export const getCalculationResults = (calculationId: number | string) => {
  return api.get(`/calculations/${calculationId}/results/`);
};

// 提交配方到后端处理
export const submitRecipe = async (recipeData: any) => {
  console.log('提交配方数据:', recipeData);
  try {
    const response = await api.post('/process-recipe/', recipeData);
    console.log('配方提交响应:', response);
    return response;  // 返回完整的response对象
  } catch (error: any) {  // 明确指定error类型为any
    console.error('配方提交失败:', error);
    if (error.response) {
      console.error('错误响应:', error.response.data);
      console.error('状态码:', error.response.status);
      console.error('请求URL:', error.config.url);
    }
    throw error;
  }
}; 