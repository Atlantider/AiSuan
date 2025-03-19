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
  return api.post(`/formulations/${id}/generate_input_file/`);
};

// 获取所有计算
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