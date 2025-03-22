import React, { useState, useEffect, useRef } from 'react';
import { Typography, Card, Row, Col, Button, Divider, Space, Tabs, Steps, Form, Input, InputNumber, Select, Upload, message, Checkbox, theme, notification, Modal, Alert, Tag, Drawer } from 'antd';
import { useNavigate } from 'react-router-dom';
import { 
  NodeIndexOutlined,
  BranchesOutlined,
  InteractionOutlined,
  AppstoreOutlined,
  BarChartOutlined,
  SettingOutlined,
  UploadOutlined,
  DownloadOutlined,
  PlusOutlined,
  LeftOutlined,
  RightOutlined,
  DeleteOutlined,
  BugOutlined,
  LoginOutlined
} from '@ant-design/icons';
import { isAuthenticated, login } from '../../services/auth';
import * as electrolyteService from '../../services/electrolyteService';
import { ICalculationData, IFormulationData } from '../../services/electrolyteService';
import api from '../../services/api';

// 解构Antd组件
const { Title, Paragraph } = Typography;
const { Step } = Steps;
const { useToken } = theme;

// 修改样式对象的类型定义
interface StylesType {
  [key: string]: React.CSSProperties;
}

// 定义样式
const styles: StylesType = {
  container: {
    padding: 24
  },
  stepContent: {
    marginBottom: 24
  },
  contentCard: {
    marginBottom: 16
  },
  actionButtons: {
    display: 'flex',
    justifyContent: 'center',
    marginTop: 24
  },
  componentRow: {
    marginBottom: 16
  },
  pageContainer: {
    padding: '24px',
    backgroundColor: '#f5f5f5'
  },
  containerCard: {
    padding: 24
  },
  innerCard: {
    marginTop: 24,
    boxShadow: '0 2px 8px rgba(0,0,0,0.08)',
    borderRadius: 8
  },
  formItem: {
    marginBottom: 16
  },
  stepButtons: {
    display: 'flex',
    justifyContent: 'space-between',
    marginTop: 24
  }
};

const ElectrolyteCalculation: React.FC = () => {
  const navigate = useNavigate();
  const [activeTab, setActiveTab] = useState('1');
  const [form] = Form.useForm();
  const [step, setStep] = useState(0);
  const { token } = useToken();
  
  // 新增状态
  const [solvents, setSolvents] = useState<any[]>([]);
  const [salts, setSalts] = useState<any[]>([]);
  const [loading, setLoading] = useState(false);
  const [creating, setCreating] = useState(false);
  const [formSubmitted, setFormSubmitted] = useState(false);
  const [currentFormulationId, setCurrentFormulationId] = useState<number | null>(null);
  const [inputFileUrl, setInputFileUrl] = useState<string | null>(null);
  const [hasInputFile, setHasInputFile] = useState(false);

  // 在组件内的状态中添加登录相关状态
  const [isLoginModalVisible, setIsLoginModalVisible] = useState(false);
  const [loginForm] = Form.useForm();
  const [loginLoading, setLoginLoading] = useState(false);

  // 创建状态
  const [generating, setGenerating] = useState(false);
  const [downloading, setDownloading] = useState(false);
  const [submitting, setSubmitting] = useState(false);
  const [loginVisible, setLoginVisible] = useState(false);

  // 添加配方列表状态
  const [formulations, setFormulations] = useState<any[]>([]);
  const [formulationsLoading, setFormulationsLoading] = useState<boolean>(false);

  // 添加日志状态和函数
  const [logMessages, setLogMessages] = useState<string[]>([]);
  const [showLogDrawer, setShowLogDrawer] = useState(false);
  const logEndRef = useRef<HTMLDivElement>(null);
  
  // 添加日志函数
  const addLog = (msg: string) => {
    const timestamp = new Date().toLocaleTimeString();
    const logMsg = `[${timestamp}] ${msg}`;
    setLogMessages(prev => [...prev, logMsg]);
    
    // 滚动到日志底部
    setTimeout(() => {
      if (logEndRef.current) {
        logEndRef.current.scrollIntoView({ behavior: 'smooth' });
      }
    }, 100);
  };

  // 渲染状态日志组件
  const renderLogDrawer = () => {
    return (
      <Drawer
        title="处理状态日志"
        placement="right"
        width={500}
        open={showLogDrawer}
        onClose={() => setShowLogDrawer(false)}
        maskClosable={true}
      >
        <div style={{ 
          height: 'calc(100vh - 150px)', 
          overflowY: 'auto',
          backgroundColor: '#f0f0f0',
          padding: '10px',
          fontFamily: 'monospace',
          borderRadius: '4px'
        }}>
          {logMessages.map((msg, index) => (
            <div key={index} style={{ 
              marginBottom: '5px',
              whiteSpace: 'pre-wrap',
              wordBreak: 'break-word'
            }}>
              {msg}
            </div>
          ))}
          <div ref={logEndRef} />
        </div>
      </Drawer>
    );
  };

  // 从后端获取溶剂和盐数据
  const fetchSolventsAndSalts = () => {
    if (!isAuthenticated()) {
      // 如果未登录，设置默认数据而不是调用API
      setSolvents([
        { id: 1, name: 'EC (ethylene carbonate)', smile: 'C1OC(=O)O1', description: '碳酸乙烯酯' },
        { id: 2, name: 'DMC (dimethyl carbonate)', smile: 'COC(=O)OC', description: '碳酸二甲酯' },
        { id: 3, name: 'EMC (ethyl methyl carbonate)', smile: 'CCOC(=O)OC', description: '碳酸甲乙酯' }
      ]);
      setSalts([
        { id: '1', name: 'LiPF6', cation: 'Li', anion: 'PF6', description: '六氟磷酸锂' },
        { id: '2', name: 'LiBF4', cation: 'Li', anion: 'BF4', description: '四氟硼酸锂' },
        { id: '3', name: 'LiTFSI', cation: 'Li', anion: 'TFSI', description: '双(三氟甲磺酰)亚胺锂' },
        { id: '4', name: 'LiFSI', cation: 'Li', anion: 'FSI', description: '双(氟磺酰)亚胺锂' },
        { id: '5', name: 'NaPF6', cation: 'Na', anion: 'PF6', description: '六氟磷酸钠' }
      ]);
      setLoading(false);
      return;
    }

    setLoading(true);
    Promise.all([electrolyteService.getSolvents(), electrolyteService.getSalts()])
      .then(([solventsRes, saltsRes]) => {
        setSolvents(solventsRes.data);
        setSalts(saltsRes.data);
      })
      .catch(error => {
        console.error('获取数据失败:', error);
        message.error('获取溶剂和盐数据失败，将使用默认数据');
        
        // 设置默认数据
        setSolvents([
          { id: 1, name: 'EC (ethylene carbonate)', smile: 'C1OC(=O)O1', description: '碳酸乙烯酯' },
          { id: 2, name: 'DMC (dimethyl carbonate)', smile: 'COC(=O)OC', description: '碳酸二甲酯' },
          { id: 3, name: 'EMC (ethyl methyl carbonate)', smile: 'CCOC(=O)OC', description: '碳酸甲乙酯' }
        ]);
        setSalts([
          { id: '1', name: 'LiPF6', cation: 'Li', anion: 'PF6', description: '六氟磷酸锂' },
          { id: '2', name: 'LiBF4', cation: 'Li', anion: 'BF4', description: '四氟硼酸锂' },
          { id: '3', name: 'LiTFSI', cation: 'Li', anion: 'TFSI', description: '双(三氟甲磺酰)亚胺锂' },
          { id: '4', name: 'LiFSI', cation: 'Li', anion: 'FSI', description: '双(氟磺酰)亚胺锂' },
          { id: '5', name: 'NaPF6', cation: 'Na', anion: 'PF6', description: '六氟磷酸钠' }
        ]);
        
        message.error('获取数据失败，请检查网络连接和API权限');
      })
      .finally(() => {
        setLoading(false);
      });
  };

  // 阳离子选项
  const cationOptions = [
    { label: 'Li+', value: 'Li' },
    { label: 'Na+', value: 'Na' },
    { label: 'K+', value: 'K' },
    { label: 'Mg2+', value: 'Mg' },
    { label: 'Ca2+', value: 'Ca' }
  ];

  // 阴离子选项
  const anionOptions = [
    { label: 'PF6-', value: 'PF6' },
    { label: 'BF4-', value: 'BF4' },
    { label: 'TFSI-', value: 'TFSI' },
    { label: 'FSI-', value: 'FSI' },
    { label: 'ClO4-', value: 'ClO4' }
  ];

  // 默认溶剂选项（当API加载失败时使用）
  const defaultSolventOptions = [
    { 
      id: '1',
      name: 'EC (碳酸乙烯酯)', 
      smile: 'C1COC(=O)O1',
      description: '环状碳酸酯'
    },
    { 
      id: '2',
      name: 'DMC (碳酸二甲酯)', 
      smile: 'COC(=O)OC',
      description: '链状碳酸酯'
    },
    { 
      id: '3',
      name: 'EMC (碳酸甲乙酯)', 
      smile: 'CCOC(=O)OC',
      description: '链状碳酸酯'
    },
    { 
      id: '4',
      name: 'PC (碳酸丙烯酯)', 
      smile: 'CC1COC(=O)O1',
      description: '环状碳酸酯'
    },
    { 
      id: '5',
      name: 'DEC (碳酸二乙酯)', 
      smile: 'CCOC(=O)OCC',
      description: '链状碳酸酯'
    },
    { 
      id: '6',
      name: 'FEC (氟代碳酸乙烯酯)', 
      smile: 'C1OC(=O)OC1F',
      description: '氟代环状碳酸酯'
    }
  ];
  
  // 默认盐选项（当API加载失败时使用）
  const defaultSaltOptions = [
    {
      id: '1',
      name: 'LiPF6',
      cation: 'Li',
      anion: 'PF6',
      description: '六氟磷酸锂'
    },
    {
      id: '2',
      name: 'LiBF4',
      cation: 'Li',
      anion: 'BF4',
      description: '四氟硼酸锂'
    },
    {
      id: '3',
      name: 'LiTFSI',
      cation: 'Li',
      anion: 'TFSI',
      description: '双(三氟甲磺酰)亚胺锂'
    },
    {
      id: '4',
      name: 'LiFSI',
      cation: 'Li',
      anion: 'FSI',
      description: '双(氟磺酰)亚胺锂'
    },
    {
      id: '5',
      name: 'NaPF6',
      cation: 'Na',
      anion: 'PF6',
      description: '六氟磷酸钠'
    }
  ];
  
  // 如果API加载失败，确保使用默认选项
  useEffect(() => {
    if (!loading && solvents.length === 0) {
      setSolvents(defaultSolventOptions);
    }
  }, [loading, solvents]);
  
  useEffect(() => {
    if (!loading && salts.length === 0) {
      setSalts(defaultSaltOptions);
    }
  }, [loading, salts]);

  // 获取动态生成的溶剂选项
  const getSolventOptions = () => {
    return solvents.map(solvent => ({
      label: `${solvent.name || solvent.id} ${solvent.smile ? `(${solvent.smile})` : ''}`,
      value: solvent.id.toString(),
      smile: solvent.smile,
      description: solvent.description || '无描述'
    }));
  };

  // 动态生成盐选项
  const saltOptions = salts.map(salt => ({
    label: `${salt.name || salt.id} ${salt.cation && salt.anion ? `(${salt.cation}/${salt.anion})` : ''}`,
    value: salt.id.toString(),
    cation: salt.cation,
    anion: salt.anion,
    description: salt.description || '无描述'
  }));

  // 计算类型选项
  const calculationOptions = [
    {
      label: '离子电导率计算',
      value: 'conductivity',
      icon: <NodeIndexOutlined style={{ fontSize: 24, color: token.colorPrimary }} />,
      description: '计算电解液中的离子电导率，评估电解液的导电能力。'
    },
    {
      label: '扩散系数计算',
      value: 'diffusion',
      icon: <BranchesOutlined style={{ fontSize: 24, color: token.colorPrimary }} />,
      description: '计算电解液中离子的扩散系数，评估离子在电解液中的迁移速率。'
    },
    {
      label: '密度与黏度计算',
      value: 'density_viscosity',
      icon: <InteractionOutlined style={{ fontSize: 24, color: token.colorPrimary }} />,
      description: '计算电解液的密度和黏度，评估电解液的基本物理性质。'
    },
    {
      label: '径向分布函数计算',
      value: 'rdf',
      icon: <BarChartOutlined style={{ fontSize: 24, color: token.colorPrimary }} />,
      description: '计算电解液中原子间的径向分布函数，分析溶液结构。'
    }
  ];

  // 定义计算类型值的类型
  type CalculationType = 'conductivity' | 'diffusion' | 'density_viscosity' | 'rdf';

  // Excel上传配置
  const uploadProps = {
    name: 'file',
    accept: '.xlsx,.xls',
    action: '/api/upload-electrolyte-formula',
    onChange(info: any) {
      if (info.file.status === 'done') {
        message.success(`${info.file.name} 上传成功`);
        handleExcelData(info.file.response);
      } else if (info.file.status === 'error') {
        message.error(`${info.file.name} 上传失败`);
      }
    },
  };

  // 处理Excel数据
  const handleExcelData = (data: any) => {
    form.setFieldsValue(data);
  };

  // 下载Excel模板
  const downloadTemplate = () => {
    window.open('/api/download-template', '_blank');
  };

  // 表单提交处理
  const handleSubmit = async (values: any) => {
    if (step === 2) {
      try {
        // 清空之前的日志
        setLogMessages([]);
        // 显示日志抽屉
        setShowLogDrawer(true);
        
        setSubmitting(true);
        addLog('开始提交配方数据...');
        addLog(`表单数据: ${JSON.stringify(values, null, 2)}`);
        
        // 检查表单数据结构
        addLog(`发现阳离子数据: ${values.cations?.length || 0}个`);
        addLog(`发现阴离子数据: ${values.anions?.length || 0}个`);
        addLog(`发现溶剂数据: ${values.solvents?.length || 0}个`);
        
        // 准备配方数据 - 完全匹配后端期望的格式
        const recipeData = {
          name: values.formulaName || "未命名配方",
          temperature: values.temperature || 300,
          box_size: values.boxSize || 30,
          // 阳离子数据
          cations: values.cations?.map((cation: any) => ({
            name: cation.type || "Li",
            number: Math.round(cation.concentration * 16),
            charge: 1
          })) || [{
            name: "Li", 
            number: 16, 
            charge: 1
          }],
          // 阴离子数据
          anions: values.anions?.map((anion: any) => ({
            name: anion.type || "PF6",
            number: Math.round(anion.concentration * 16),
            charge: -1
          })) || [{
            name: "PF6", 
            number: 16, 
            charge: -1
          }],
          // 溶剂数据
          solvents: values.solvents?.map((s: any) => {
            const solventOption = getSolventOptions().find(opt => opt.value === s.solvent);
            return {
              name: solventOption?.label.split(' ')[0] || 'EC',
              smile: solventOption?.smile || 'C1COC(=O)O1',
              number: Math.round(s.volumeRatio * 54)
            };
          }) || [{
            name: "EC", 
            smile: "C1COC(=O)O1", 
            number: 54
          }],
          // 可选：传递计算类型
          calculation_types: values.calculationTypes || ['conductivity']
        };

        addLog(`准备提交的配方数据: ${JSON.stringify(recipeData, null, 2)}`);
        
        try {
          // 显示加载提示
          message.loading({ content: '正在提交配方到后端...', key: 'submitRecipe', duration: 0 });
          
          // 调用API提交配方前的日志
          addLog('开始调用API提交配方...');
          
          // 直接使用原始axios发送请求，添加更多调试信息
          addLog(`提交到URL: ${api.defaults.baseURL}/process-recipe/`);
          
          // 调用API提交配方
          addLog('发送API请求...');
          try {
            const response = await electrolyteService.submitRecipe(recipeData);
            
            // 调用API后的日志
            addLog('API请求完成!');
            addLog(`响应状态: ${response.status}`);
            addLog(`响应数据: ${JSON.stringify(response.data, null, 2)}`);
            
            if (response?.data?.success) {
              message.success({ content: '配方提交成功！', key: 'submitRecipe' });
              addLog('配方提交成功!');
              
              // 保存创建的配方ID
              if (response.data.formulation_id) {
                setCurrentFormulationId(response.data.formulation_id);
                addLog(`设置当前配方ID: ${response.data.formulation_id}`);
                
                // 显示弹窗表明成功
                Modal.success({
                  title: '配方提交成功',
                  content: (
                    <div>
                      <p>配方ID: {response.data.formulation_id}</p>
                      {response.data.calculation_id && (
                        <p>计算任务ID: {response.data.calculation_id}</p>
                      )}
                      <p>配方已成功提交到后端处理</p>
                    </div>
                  ),
                });
              }
              
              // 如果返回了计算ID，说明计算任务已经创建
              if (response.data.calculation_id) {
                addLog(`收到计算任务ID: ${response.data.calculation_id}`);
                notification.success({
                  message: '计算任务已创建',
                  description: `计算任务ID: ${response.data.calculation_id}，正在跳转到计算任务列表...`,
                  duration: 5
                });
                
                // 跳转到计算任务列表页面
                addLog('延迟3秒后将跳转到任务管理页面...');
                setTimeout(() => {
                  navigate('/user/tasks', { 
                    state: { 
                      refresh: true,
                      highlightId: response.data.calculation_id 
                    } 
                  });
                }, 3000);
              }
            } else {
              // 如果返回的success为false
              addLog(`请求返回success=false: ${response.data.error || '未知错误'}`);
              throw new Error(response?.data?.error || '提交失败，未收到成功响应');
            }
          } catch (apiError: any) {
            addLog(`API调用出错: ${apiError.message}`);
            
            if (apiError.response) {
              addLog(`HTTP状态码: ${apiError.response.status}`);
              addLog(`响应数据: ${JSON.stringify(apiError.response.data || {}, null, 2)}`);
            } else if (apiError.request) {
              addLog('请求已发送但未收到响应，可能是网络问题或后端未启动');
            } else {
              addLog(`请求配置问题: ${apiError.message}`);
            }
            
            message.error({
              content: apiError.response?.data?.error || apiError.message || '配方提交失败，请检查网络连接',
              key: 'submitRecipe'
            });
            
            // 显示详细错误信息
            notification.error({
              message: '提交失败',
              description: (
                <div>
                  <p>错误详情：{apiError.message}</p>
                  <p>请求URL：{apiError.config?.url || '未知'}</p>
                  <p>状态码：{apiError.response?.status || '未知'}</p>
                  <p>响应数据：{JSON.stringify(apiError.response?.data || '无数据')}</p>
                </div>
              ),
              duration: 0
            });
            
            // 显示错误弹窗
            Modal.error({
              title: '提交失败',
              content: (
                <div>
                  <p>提交配方时出现错误：</p>
                  <p>{apiError.message || '未知错误'}</p>
                  <p>请检查网络连接和服务器状态，然后重试。</p>
                  <p>详细日志请查看右侧状态日志面板。</p>
                </div>
              ),
            });
          }
        } catch (error: any) {
          addLog(`处理过程中出现异常: ${error.message}`);
          addLog(`错误堆栈: ${error.stack}`);
          
          message.error({
            content: '提交失败，请检查数据格式',
            key: 'submitRecipe'
          });
        } finally {
          setSubmitting(false);
          addLog('提交流程结束');
        }
      } catch (error: any) {
        console.error('提交处理失败:', error);
        console.error('错误堆栈:', error.stack);
        message.error({
          content: '提交失败，请检查数据格式',
          key: 'submitRecipe'
        });
      }
    } else {
      // 如果不是最后一步，继续下一步
      setStep(step + 1);
    }
  };

  // 根据阴阳离子获取盐的ID
  const getSaltIdByCationAnion = (cation: string, anion: string): number | null => {
    // 这里可以实现一个根据阳离子和阴离子查找盐ID的逻辑
    // 简单起见，这里返回null，表示需要后端创建新的盐
    return null;
  };

  // 组件挂载时执行
  useEffect(() => {
    // 加载数据
    fetchSolventsAndSalts();
    loadFormulations(); // 添加此行加载配方列表
    
    // 初始化Form默认值
    form.setFieldsValue({
      cations: [{ type: 'Li', concentration: 1.0 }],
      anions: [{ type: 'PF6', concentration: 1.0 }],
      solvents: [{ solvent: '1', volumeRatio: 1.0 }]
    });
  }, []);

  // 生成LAMMPS输入文件
  const generateLammpsInputs = async () => {
    // 添加UI反馈
    message.loading({ content: '正在生成输入文件...', key: 'generateFile', duration: 0 });
    console.log('========开始执行generateLammpsInputs函数=========');
    console.log('生成输入文件按钮被点击');
    
    try {
      console.log('设置generating状态为true');
      setGenerating(true);
      
      // 获取表单数据
      console.log('开始获取表单数据');
      const values = form.getFieldsValue();
      console.log('表单数据:', values);
      
      console.log('开始生成molyte-cursor输入文件');
      
      // 创建INP文件内容
      let fileContent = '';
      
      // 添加开始标识符
      fileContent += `START\n`;
      
      // 添加基本信息
      fileContent += `name ${values.formulaName || '未命名配方'}\n`;
      fileContent += `T ${values.temperature || 298}\n`;
      fileContent += `Box_size ${values.boxSize || 50}\n`;
      fileContent += `concentration ${values.cations && values.cations[0] ? values.cations[0].concentration || 1.0 : 1.0}\n`;
      fileContent += `time_step ${values.timeStep || 1.0}\n`;
      fileContent += `equilibration_steps ${values.equilibrationSteps || 1000000}\n`;
      fileContent += `production_steps ${values.productionSteps || 2000000}\n`;
      fileContent += `cutoff ${values.cutoff || 12.0}\n`;
      fileContent += `pressure ${values.pressure || 1.0}\n`;
      
      // 检查阳离子数据并添加
      console.log('开始处理阳离子数据...');
      const cations = values.cations || [];
      console.log('阳离子数据:', JSON.stringify(cations, null, 2));
      
      if (cations.length > 0) {
        cations.forEach((cation: any, index: number) => {
          const i = index + 1;
          const ratio = index === 0 ? 1.0 : (cation.concentration / (cations[0]?.concentration || 1.0));
          
          // 只添加名称和比例，去掉null和smile
          fileContent += `Cation${i}_name ${cation.type || 'Li'}\n`;
          fileContent += `Cation${i}_ratio ${parseFloat(ratio.toFixed(2))}\n`;
          
          console.log(`已添加Cation${i}相关信息`);
        });
      } else {
        // 如果没有阳离子数据，添加默认值
        fileContent += `Cation1_name Li\n`;
        fileContent += `Cation1_ratio 1.0\n`;
        
        console.log('未找到阳离子数据，已添加默认值');
      }
      
      // 检查阴离子数据并添加
      console.log('开始处理阴离子数据...');
      const anions = values.anions || [];
      console.log('阴离子数据:', JSON.stringify(anions, null, 2));
      
      if (anions.length > 0) {
        anions.forEach((anion: any, index: number) => {
          const i = index + 1;
          const cationConcentration = cations[0]?.concentration || 1.0;
          const ratio = anion.concentration / cationConcentration;
          
          // 只添加名称和比例，去掉null和smile
          fileContent += `Anion${i}_name ${anion.type || 'PF6'}\n`;
          fileContent += `Anion${i}_ratio ${parseFloat(ratio.toFixed(2))}\n`;
          
          console.log(`已添加Anion${i}相关信息`);
        });
      } else {
        // 如果没有阴离子数据，添加默认值
        fileContent += `Anion1_name PF6\n`;
        fileContent += `Anion1_ratio 1.0\n`;
        
        console.log('未找到阴离子数据，已添加默认值');
      }
      
      // 检查溶剂数据并添加
      console.log('开始处理溶剂数据...');
      const solvents = values.solvents || [];
      console.log('溶剂数据:', JSON.stringify(solvents, null, 2));
      console.log('溶剂选项:', JSON.stringify(getSolventOptions(), null, 2));
      
      if (solvents.length > 0) {
        solvents.forEach((solvent: any, index: number) => {
          const i = index + 1;
          const solventOption = getSolventOptions().find(opt => opt.value === solvent.solvent);
          console.log(`查找溶剂${i}:`, solvent.solvent, '结果:', solventOption);
          
          // 从溶剂选项中提取名称和SMILE
          const solventName = solventOption ? 
            (solventOption.label.split(' ')[0]) : `Solvent${i}`;
          const solventSmile = solventOption?.smile || '';
          const ratio = (solvent.volumeRatio || 1.0) * 5; // 溶剂比例通常更高
          
          // 保留溶剂的完整信息
          fileContent += `Sol${i}_name ${solventName}\n`;
          fileContent += `Sol${i}_smile ${solventSmile}\n`;
          fileContent += `Sol${i}_ratio ${parseFloat(ratio.toFixed(2))}\n`;
          
          console.log(`已添加Sol${i}相关信息:`, { name: solventName, smile: solventSmile, ratio });
        });
      } else {
        // 如果没有溶剂数据，添加默认值
        fileContent += `Sol1_name EC\n`;
        fileContent += `Sol1_smile C1COC(=O)O1\n`;
        fileContent += `Sol1_ratio 5.0\n`;
        
        console.log('未找到溶剂数据，已添加默认值');
      }
      
      // 添加计算类型
      fileContent += `calculation_types ${((values.calculationTypes || ['conductivity']) as CalculationType[]).join(',')}\n`;
      
      // 添加结束标识符
      fileContent += `END\n`;
      
      console.log('生成的inp文件内容:');
      console.log(fileContent);
      
      try {
        let targetFormulationId = currentFormulationId;
        
        // 如果没有当前配方ID，则先创建新配方
        if (!targetFormulationId) {
          console.log('没有选择现有配方，将创建新配方...');
          message.loading({ content: '正在创建新配方...', key: 'generateFile', duration: 0 });
          
          // 准备创建配方的数据
          const formData = {
            name: values.formulaName || '临时电解液配方',
            description: values.description || '通过前端自动创建的临时配方',
            salts: [
              {
                id: 1,
                concentration: 1.0
              }
            ],
            solvents: (values.solvents || []).map((s: any) => ({
              id: parseInt(s.solvent) || 1,
              concentration: s.volumeRatio || 1.0
            })),
            temperature: values.temperature || 298,
            pressure: values.pressure || 1.0,
            time_step: values.timeStep || 1.0,
            equilibration_steps: values.equilibrationSteps || 1000000,
            production_steps: values.productionSteps || 2000000,
            cutoff: values.cutoff || 12.0
          };
          
          try {
            const response = await electrolyteService.createFormulation(formData);
            console.log('创建配方响应:', response);
            
            if (response && response.data && response.data.id) {
              // 获取真实配方ID
              targetFormulationId = response.data.id;
              console.log('获取到新配方ID:', targetFormulationId);
              
              // 保存配方ID到状态
              setCurrentFormulationId(targetFormulationId);
              message.success(`已创建新配方(ID=${targetFormulationId})`);
            } else {
              throw new Error('创建配方成功但未返回配方ID');
            }
          } catch (createError: any) {
            console.error('创建配方失败:', createError);
            message.error({
              content: '创建配方失败，请检查数据格式或网络连接',
              key: 'generateFile'
            });
            if (createError.response) {
              console.error('错误响应:', createError.response.data);
              notification.error({
                message: '创建配方失败',
                description: `详细错误: ${JSON.stringify(createError.response.data)}`
              });
            }
            setGenerating(false);
            return;
          }
        }
        
        // 确保此时targetFormulationId一定有值
        if (!targetFormulationId) {
          message.error({
            content: '无法获取有效的配方ID',
            key: 'generateFile'
          });
          setGenerating(false);
          return;
        }
        
        // 使用获取到的配方ID保存INP文件
        console.log('使用配方ID:', targetFormulationId);
        message.loading({ content: '正在保存输入文件...', key: 'generateFile', duration: 0 });
        
        try {
          // 更新INP文件
          const saveResponse = await electrolyteService.saveInputFile(targetFormulationId, fileContent);
          console.log('保存INP文件响应:', saveResponse);
          
          // 调用后端生成输入文件API
          console.log('调用后端生成输入文件API...');
          message.loading({ content: '正在生成LAMMPS输入文件...', key: 'generateFile', duration: 0 });
          
          try {
            console.log(`准备调用generateInputFile API，配方ID=${targetFormulationId}`);
            const generateResponse = await electrolyteService.generateInputFile(targetFormulationId);
            console.log('生成输入文件响应:', generateResponse);
            
            if (generateResponse.status === 200) {
              setHasInputFile(true);
              message.success({
                content: '成功生成输入文件!',
                key: 'generateFile'
              });
              notification.success({
                message: '成功生成输入文件',
                description: `输入文件已成功生成并保存到服务器。配方ID: ${targetFormulationId}`
              });
            } else {
              throw new Error('生成输入文件API调用返回非200状态码');
            }
          } catch (genError: any) {
            console.error('调用生成输入文件API失败:', genError);
            console.error('详细错误:', genError.response?.data || '无详细信息');
            
            // 显示错误提示
            message.error({
              content: '生成输入文件失败，但已保存基本内容',
              key: 'generateFile'
            });
            
            // 显示详细错误信息
            notification.error({
              message: '生成输入文件失败',
              description: (
                <div>
                  <p>API调用失败但已保存基本输入文件内容。</p>
                  <p>错误详情: {genError.response?.data?.error || genError.message || '未知错误'}</p>
                  <p>配方ID: {targetFormulationId}</p>
                </div>
              ),
              duration: 0
            });
            
            // 即使生成API失败，也设置为true以允许用户继续
            setHasInputFile(true);
          }
        } catch (apiError: any) {
          console.error('保存INP文件失败:', apiError);
          
          // 显示错误提示
          message.error({
            content: '保存INP文件失败',
            key: 'generateFile'
          });
          
          if (apiError.response) {
            console.error('错误响应:', apiError.response.data);
            notification.error({
              message: '保存INP文件失败',
              description: `错误详情: ${JSON.stringify(apiError.response.data)}`,
              duration: 0
            });
          }
        }
      } catch (error: any) {
        console.error('生成输入文件过程中出错:', error);
        message.error({
          content: '生成过程中发生错误',
          key: 'generateFile'
        });
        notification.error({
          message: '生成输入文件错误',
          description: error.message || '未知错误',
          duration: 0
        });
      } finally {
        setGenerating(false);
      }
    } catch (error: any) {
      console.error('生成输入文件出错:', error);
      message.error({
        content: '生成输入文件失败',
        key: 'generateFile'
      });
      setGenerating(false);
    }
  };
  
  // 提交计算
  const submitCalculation = async () => {
    try {
      setSubmitting(true);
      
      const values = form.getFieldsValue();
      
      if (!currentFormulationId) {
        setSubmitting(false);
        return;
      }
      
      // 在提交计算前先保存最新的INP文件内容
      try {
        // 创建INP文件内容
        let fileContent = '';
        
        // 添加开始标识符
        fileContent += `START\n`;
        
        // 添加基本信息
        fileContent += `name ${values.formulaName || '未命名配方'}\n`;
        fileContent += `T ${values.temperature || 298}\n`;
        fileContent += `Box_size ${values.boxSize || 50}\n`;
        fileContent += `concentration ${values.cations && values.cations[0] ? values.cations[0].concentration || 1.0 : 1.0}\n`;
        fileContent += `time_step ${values.timeStep || 1.0}\n`;
        fileContent += `equilibration_steps ${values.equilibrationSteps || 1000000}\n`;
        fileContent += `production_steps ${values.productionSteps || 2000000}\n`;
        fileContent += `cutoff ${values.cutoff || 12.0}\n`;
        fileContent += `pressure ${values.pressure || 1.0}\n`;
        
        // 添加阳离子数据
        const cations = values.cations || [];
        if (cations.length > 0) {
          cations.forEach((cation: any, index: number) => {
            const i = index + 1;
            const ratio = index === 0 ? 1.0 : (cation.concentration / (cations[0]?.concentration || 1.0));
            fileContent += `Cation${i}_name ${cation.type || 'Li'}\n`;
            fileContent += `Cation${i}_ratio ${parseFloat(ratio.toFixed(2))}\n`;
          });
        } else {
          fileContent += `Cation1_name Li\n`;
          fileContent += `Cation1_ratio 1.0\n`;
        }
        
        // 添加阴离子数据
        const anions = values.anions || [];
        if (anions.length > 0) {
          anions.forEach((anion: any, index: number) => {
            const i = index + 1;
            const cationConcentration = cations[0]?.concentration || 1.0;
            const ratio = anion.concentration / cationConcentration;
            fileContent += `Anion${i}_name ${anion.type || 'PF6'}\n`;
            fileContent += `Anion${i}_ratio ${parseFloat(ratio.toFixed(2))}\n`;
          });
        } else {
          fileContent += `Anion1_name PF6\n`;
          fileContent += `Anion1_ratio 1.0\n`;
        }
        
        // 添加溶剂数据
        const solvents = values.solvents || [];
        if (solvents.length > 0) {
          solvents.forEach((solvent: any, index: number) => {
            const i = index + 1;
            const solventOption = getSolventOptions().find(opt => opt.value === solvent.solvent);
            const solventName = solventOption ? 
              (solventOption.label.split(' ')[0]) : `Solvent${i}`;
            const solventSmile = solventOption?.smile || '';
            const ratio = (solvent.volumeRatio || 1.0) * 5; // 溶剂比例通常更高
            fileContent += `Sol${i}_name ${solventName}\n`;
            fileContent += `Sol${i}_smile ${solventSmile}\n`;
            fileContent += `Sol${i}_ratio ${parseFloat(ratio.toFixed(2))}\n`;
          });
        } else {
          fileContent += `Sol1_name EC\n`;
          fileContent += `Sol1_smile C1COC(=O)O1\n`;
          fileContent += `Sol1_ratio 5.0\n`;
        }
        
        // 添加计算类型
        fileContent += `calculation_types ${((values.calculationTypes || ['conductivity']) as CalculationType[]).join(',')}\n`;
        
        // 添加结束标识符
        fileContent += `END\n`;
        
        // 保存INP文件到后端
        await electrolyteService.saveInputFile(currentFormulationId, fileContent);
        console.log('已更新INP文件到后端');
      } catch (fileError) {
        console.error('更新INP文件失败:', fileError);
      }
      
      const calculationData: ICalculationData = {
        name: values.calculationName || `电解液计算${new Date().toLocaleDateString('zh-CN')}`,
        description: values.calculationDescription || '',
        formulation_id: currentFormulationId
      };
      
      console.log('提交计算数据:', calculationData);
      
      try {
        // 调用API提交计算任务
        const response = await electrolyteService.submitElectrolyteCalculation(calculationData);
        console.log('计算任务提交响应:', response);
        
        // 立即跳转到任务管理页面
        navigate('/user/tasks', { state: { refresh: true } }); // 跳转到任务管理页面并带上刷新标志
      } catch (apiError) {
        console.error('API调用失败:', apiError);
      }
      
      setSubmitting(false);
    } catch (error) {
      console.error('提交计算失败:', error);
      setSubmitting(false);
    }
  };

  // 添加下载输入文件的函数
  const downloadInputFile = async () => {
    console.log('下载输入文件按钮被点击');
    
    if (!hasInputFile) {
      console.log('没有可用的输入文件，需要先生成');
      return;
    }
    
    try {
      // 获取表单数据
      const values = form.getFieldsValue(true); // 使用true获取所有字段值
      console.log('下载文件时的表单数据:', JSON.stringify(values, null, 2));
      
      // 创建文本内容字符串，每行一个参数
      let fileContent = '';
      
      // 添加开始标识符
      fileContent += `START\n`;
      
      // 添加基本信息
      fileContent += `name ${values.formulaName || '未命名配方'}\n`;
      fileContent += `T ${values.temperature || 298}\n`;
      fileContent += `Box_size ${values.boxSize || 50}\n`;
      fileContent += `concentration ${values.cations && values.cations[0] ? values.cations[0].concentration || 1.0 : 1.0}\n`;
      fileContent += `time_step ${values.timeStep || 1.0}\n`;
      fileContent += `equilibration_steps ${values.equilibrationSteps || 1000000}\n`;
      fileContent += `production_steps ${values.productionSteps || 2000000}\n`;
      fileContent += `cutoff ${values.cutoff || 12.0}\n`;
      fileContent += `pressure ${values.pressure || 1.0}\n`;
      
      // 检查阳离子数据并添加（不包含null和smile）
      console.log('开始处理阳离子数据...');
      const cations = values.cations || [];
      console.log('阳离子数据:', JSON.stringify(cations, null, 2));
      
      if (cations.length > 0) {
        cations.forEach((cation: any, index: number) => {
          const i = index + 1;
          const ratio = index === 0 ? 1.0 : (cation.concentration / (cations[0]?.concentration || 1.0));
          
          // 只添加名称和比例，去掉null和smile
          fileContent += `Cation${i}_name ${cation.type || 'Li'}\n`;
          fileContent += `Cation${i}_ratio ${parseFloat(ratio.toFixed(2))}\n`;
          
          console.log(`已添加Cation${i}相关信息（已去除null和smile）`);
        });
      } else {
        // 如果没有阳离子数据，添加默认值（不包含null和smile）
        fileContent += `Cation1_name Li\n`;
        fileContent += `Cation1_ratio 1.0\n`;
        
        console.log('未找到阳离子数据，已添加默认值（已去除null和smile）');
      }
      
      // 检查阴离子数据并添加（不包含null和smile）
      console.log('开始处理阴离子数据...');
      const anions = values.anions || [];
      console.log('阴离子数据:', JSON.stringify(anions, null, 2));
      
      if (anions.length > 0) {
        anions.forEach((anion: any, index: number) => {
          const i = index + 1;
          const cationConcentration = cations[0]?.concentration || 1.0;
          const ratio = anion.concentration / cationConcentration;
          
          // 只添加名称和比例，去掉null和smile
          fileContent += `Anion${i}_name ${anion.type || 'PF6'}\n`;
          fileContent += `Anion${i}_ratio ${parseFloat(ratio.toFixed(2))}\n`;
          
          console.log(`已添加Anion${i}相关信息（已去除null和smile）`);
        });
      } else {
        // 如果没有阴离子数据，添加默认值（不包含null和smile）
        fileContent += `Anion1_name PF6\n`;
        fileContent += `Anion1_ratio 1.0\n`;
        
        console.log('未找到阴离子数据，已添加默认值（已去除null和smile）');
      }
      
      // 检查溶剂数据并添加
      console.log('开始处理溶剂数据...');
      const solvents = values.solvents || [];
      console.log('溶剂数据:', JSON.stringify(solvents, null, 2));
      console.log('溶剂选项:', JSON.stringify(getSolventOptions(), null, 2));
      
      if (solvents.length > 0) {
        solvents.forEach((solvent: any, index: number) => {
          const i = index + 1;
          const solventOption = getSolventOptions().find(opt => opt.value === solvent.solvent);
          console.log(`查找溶剂${i}:`, solvent.solvent, '结果:', solventOption);
          
          // 从溶剂选项中提取名称和SMILE
          const solventName = solventOption ? 
            (solventOption.label.split(' ')[0]) : `Solvent${i}`;
          const solventSmile = solventOption?.smile || '';
          const ratio = (solvent.volumeRatio || 1.0) * 5; // 溶剂比例通常更高
          
          // 保留溶剂的完整信息（但去掉null）
          fileContent += `Sol${i}_name ${solventName}\n`;
          fileContent += `Sol${i}_smile ${solventSmile}\n`;
          fileContent += `Sol${i}_ratio ${parseFloat(ratio.toFixed(2))}\n`;
          
          console.log(`已添加Sol${i}相关信息（已去除null）:`, { name: solventName, smile: solventSmile, ratio });
        });
      } else {
        // 如果没有溶剂数据，添加默认值（不包含null）
        fileContent += `Sol1_name EC\n`;
        fileContent += `Sol1_smile C1COC(=O)O1\n`;
        fileContent += `Sol1_ratio 5.0\n`;
        
        console.log('未找到溶剂数据，已添加默认值（已去除null）');
      }
      
      // 添加计算类型
      fileContent += `calculation_types ${((values.calculationTypes || ['conductivity']) as CalculationType[]).join(',')}\n`;
      
      // 添加结束标识符
      fileContent += `END\n`;
      
      console.log('生成的inp文件内容:');
      console.log(fileContent);
      
      // 创建Blob对象
      const blob = new Blob([fileContent], { type: 'text/plain' });
      
      // 创建下载链接
      const url = window.URL.createObjectURL(blob);
      console.log('创建下载URL:', url);
      
      // 创建下载链接并触发点击
      try {
        const a = document.createElement('a');
        a.href = url;
        a.download = `electrolyte_input_${currentFormulationId || 'new'}.inp`;
        a.style.display = 'none';
        
        // 添加到文档并触发点击
        document.body.appendChild(a);
        console.log('添加下载链接到文档并触发点击');
        a.click();
        
        // 清理
        setTimeout(() => {
          document.body.removeChild(a);
          window.URL.revokeObjectURL(url);
          console.log('清理下载链接和URL');
        }, 100);
      } catch (error) {
        console.error('文件下载失败:', error);
      }
      
    } catch (error) {
      console.error('下载输入文件过程中出错:', error);
      message.error('下载输入文件失败，请稍后重试');
    }
  };

  // 渲染步骤内容
  const renderStepContent = () => {
    switch (step) {
      case 0:
        return renderFormulaDesign();
      case 1:
        return renderConditionsSettings();
      case 2:
  return (
          <>
            {renderCalculationSelection()}
            
            <div style={{ display: 'flex', justifyContent: 'center', marginTop: 24 }}>
              <Space size="large">
          <Button 
            type="primary" 
                  onClick={() => {
                    console.log('生成输入文件按钮被点击');
                    generateLammpsInputs();
                  }} 
                  loading={generating}
                  icon={<SettingOutlined />}
            size="large" 
          >
                  生成输入文件
          </Button>
                
          <Button 
                  type="default" 
                  onClick={downloadInputFile}
                  disabled={!hasInputFile}
                  icon={<DownloadOutlined />}
            size="large"
                >
                  下载输入文件
                </Button>
                
                <Button 
                  type="primary" 
                  onClick={submitCalculation} 
                  loading={submitting}
                  disabled={!hasInputFile}
                  icon={<AppstoreOutlined />}
                  size="large"
                >
                  提交计算
          </Button>
        </Space>
      </div>
          </>
        );
      default:
        return null;
    }
  };

  // 渲染配方设计（包含阴阳离子和溶剂）
  const renderFormulaDesign = () => {
    return (
      <>
        <Card 
          title={<Title level={4}>电解液配方设计</Title>} 
          style={styles.contentCard}
          extra={
            <Space>
              <Upload {...uploadProps}>
                <Button icon={<UploadOutlined />}>导入Excel</Button>
              </Upload>
              <Button 
                icon={<DownloadOutlined />}
                onClick={downloadTemplate}
              >
                下载模板
              </Button>
            </Space>
          }
        >
          <Form.Item
            name="formulaName"
            label="配方名称"
            rules={[{ required: true, message: '请输入配方名称' }]}
          >
            <Input placeholder="例如：LiPF6 in EC/DMC" style={{ width: '100%' }} />
          </Form.Item>
          
          <Form.Item
            name="description"
            label="配方描述"
          >
            <Input.TextArea placeholder="请输入该电解液配方的描述信息（可选）" rows={2} />
          </Form.Item>
        </Card>
        
        <Divider />
        
        {/* 使用新的阳离子和阴离子选择组件 */}
        {renderCationSelection()}
        
        <Divider />
        
        {renderAnionSelection()}
        
        <Divider />
        
        {renderSolventSelection()}
      </>
    );
  };

  // 渲染阳离子选择
  const renderCationSelection = () => {
    return (
      <Card title="阳离子选择" style={styles.contentCard}>
        <Form.List name="cations">
          {(fields, { add, remove }) => (
            <>
              <Typography.Paragraph>
                请选择阳离子并设置浓度(mol/L)
              </Typography.Paragraph>
              
              {fields.map(field => (
                <Row key={field.key} gutter={16} align="middle" style={styles.componentRow}>
                  <Col span={12}>
                    <Form.Item
                      {...field}
                      name={[field.name, 'type']}
                      label="阳离子"
                      rules={[{ required: true, message: '请选择阳离子' }]}
                    >
                      <Select options={cationOptions} placeholder="选择阳离子" />
                    </Form.Item>
                  </Col>
                  
                  <Col span={8}>
                    <Form.Item
                      {...field}
                      name={[field.name, 'concentration']}
                      label="浓度 (mol/L)"
                      rules={[{ required: true, message: '请输入浓度' }]}
                    >
                      <InputNumber min={0.01} max={10} step={0.01} precision={2} style={{ width: '100%' }} />
                    </Form.Item>
                  </Col>
                  
                  <Col span={4} style={{ textAlign: 'center', paddingTop: '30px' }}>
                    {fields.length > 1 ? (
                      <Button 
                        type="text" 
                        danger 
                        icon={<DeleteOutlined />} 
                        onClick={() => remove(field.name)}
                      />
                    ) : null}
                  </Col>
                </Row>
              ))}
              
              <Form.Item>
                <Button 
                  type="dashed" 
                  onClick={() => add({ type: undefined, concentration: 1.0 })} 
                  block 
                  icon={<PlusOutlined />}
                >
                  添加阳离子
                </Button>
              </Form.Item>
            </>
          )}
        </Form.List>
      </Card>
    );
  };

  // 渲染阴离子选择
  const renderAnionSelection = () => {
    return (
      <Card title="阴离子选择" style={styles.contentCard}>
        <Form.List name="anions">
          {(fields, { add, remove }) => (
            <>
              <Typography.Paragraph>
                请选择阴离子并设置浓度(mol/L)
              </Typography.Paragraph>
              
              {fields.map(field => (
                <Row key={field.key} gutter={16} align="middle" style={styles.componentRow}>
                  <Col span={12}>
                    <Form.Item
                      {...field}
                      name={[field.name, 'type']}
                      label="阴离子"
                      rules={[{ required: true, message: '请选择阴离子' }]}
                    >
                      <Select options={anionOptions} placeholder="选择阴离子" />
                    </Form.Item>
                  </Col>
                  
                  <Col span={8}>
                    <Form.Item
                      {...field}
                      name={[field.name, 'concentration']}
                      label="浓度 (mol/L)"
                      rules={[{ required: true, message: '请输入浓度' }]}
                    >
                      <InputNumber min={0.01} max={10} step={0.01} precision={2} style={{ width: '100%' }} />
                    </Form.Item>
                  </Col>
                  
                  <Col span={4} style={{ textAlign: 'center', paddingTop: '30px' }}>
                    {fields.length > 1 ? (
                      <Button 
                        type="text" 
                        danger 
                        icon={<DeleteOutlined />} 
                        onClick={() => remove(field.name)}
                      />
                    ) : null}
                  </Col>
                </Row>
              ))}
              
              <Form.Item>
                <Button 
                  type="dashed" 
                  onClick={() => add({ type: undefined, concentration: 1.0 })} 
                  block 
                  icon={<PlusOutlined />}
                >
                  添加阴离子
                </Button>
              </Form.Item>
            </>
          )}
        </Form.List>
      </Card>
    );
  };

  // 渲染溶剂选择
  const renderSolventSelection = () => {
    return (
      <Card title="溶剂选择" style={styles.contentCard}>
        <Form.List name="solvents">
          {(fields, { add, remove }) => (
            <>
              <Typography.Paragraph>
                请选择溶剂并设置体积比
              </Typography.Paragraph>
              
              {fields.map(field => (
                <Row key={field.key} gutter={16} align="middle" style={styles.componentRow}>
                  <Col span={12}>
                    <Form.Item
                      {...field}
                      name={[field.name, 'solvent']}
                      label="溶剂"
                      rules={[{ required: true, message: '请选择溶剂' }]}
                    >
                      <Select 
                        options={getSolventOptions()} 
                        placeholder="选择溶剂" 
                        loading={loading}
                      />
                    </Form.Item>
                  </Col>
                  
                  <Col span={8}>
                    <Form.Item
                      {...field}
                      name={[field.name, 'volumeRatio']}
                      label="体积比"
                      rules={[{ required: true, message: '请输入体积比' }]}
                    >
                      <InputNumber min={0.01} max={10} step={0.01} precision={2} style={{ width: '100%' }} />
                    </Form.Item>
                  </Col>
                  
                  <Col span={4} style={{ textAlign: 'center', paddingTop: '30px' }}>
                    {fields.length > 1 ? (
                      <Button 
                        type="text" 
                        danger 
                        icon={<DeleteOutlined />} 
                        onClick={() => remove(field.name)}
                      />
                    ) : null}
                  </Col>
                </Row>
              ))}
              
              <Form.Item>
                <Button 
                  type="dashed" 
                  onClick={() => add({ solvent: undefined, volumeRatio: 1.0 })} 
                  block 
                  icon={<PlusOutlined />}
                >
                  添加溶剂
                </Button>
              </Form.Item>
            </>
          )}
        </Form.List>
      </Card>
    );
  };

  // 渲染条件设置（温度和模拟参数）
  const renderConditionsSettings = () => {
    return (
      <Card 
        title={<Title level={4}>模拟条件设置</Title>} 
        style={styles.contentCard}
      >
        <Row gutter={[24, 16]}>
          <Col span={12}>
            <Form.Item
              name="temperature"
              label="温度"
              rules={[{ required: true, message: '请输入温度' }]}
              initialValue={298}
            >
              <InputNumber 
                min={100} 
                max={500} 
                placeholder="温度" 
                style={{ width: '100%' }}
                addonAfter="K"
              />
            </Form.Item>
          </Col>
          <Col span={12}>
            <Form.Item
              name="pressure"
              label="压力"
              rules={[{ required: true, message: '请输入压力' }]}
              initialValue={1}
            >
              <InputNumber 
                min={0.1} 
                max={100} 
                step={0.1} 
                placeholder="压力" 
                style={{ width: '100%' }}
                addonAfter="atm"
              />
            </Form.Item>
          </Col>
          
          {/* 新增的Box Size设置 */}
          <Col span={12}>
            <Form.Item
              name="boxSize"
              label="模拟盒子尺寸"
              rules={[{ required: true, message: '请输入模拟盒子尺寸' }]}
              initialValue={50}
              tooltip="模拟系统的盒子边长"
            >
              <InputNumber 
                min={20} 
                max={200} 
                step={5} 
                placeholder="盒子尺寸" 
                style={{ width: '100%' }}
                addonAfter="Å"
              />
            </Form.Item>
          </Col>
        </Row>

        <Card 
          type="inner" 
          title={<span style={styles.cardTitle}>分子动力学参数</span>}
          style={styles.innerCard}
          bordered={false}
          className="inner-card"
        >
          <Row gutter={[24, 16]}>
            <Col span={8}>
              <Form.Item
                name="timeStep"
                label="时间步长"
                rules={[{ required: true, message: '请输入时间步长' }]}
                initialValue={1.0}
              >
                <InputNumber 
                  min={0.1} 
                  max={2.0} 
                  step={0.1} 
                  placeholder="时间步长" 
                  style={{ width: '100%' }}
                  addonAfter="fs"
                />
              </Form.Item>
            </Col>
            <Col span={8}>
              <Form.Item
                name="equilibrationSteps"
                label="平衡步数"
                rules={[{ required: true, message: '请输入平衡步数' }]}
                initialValue={1000000}
              >
                <InputNumber 
                  min={10000} 
                  max={10000000} 
                  step={10000} 
                  placeholder="平衡步数" 
                  style={{ width: '100%' }}
                />
              </Form.Item>
            </Col>
            <Col span={8}>
              <Form.Item
                name="productionSteps"
                label="生产步数"
                rules={[{ required: true, message: '请输入生产步数' }]}
                initialValue={2000000}
              >
                <InputNumber 
                  min={10000} 
                  max={10000000} 
                  step={10000} 
                  placeholder="生产步数" 
                  style={{ width: '100%' }}
                />
              </Form.Item>
            </Col>
          </Row>
          <Row gutter={[24, 16]}>
            <Col span={8}>
              <Form.Item
                name="cutoff"
                label="截断半径"
                rules={[{ required: true, message: '请输入截断半径' }]}
                initialValue={12.0}
              >
                <InputNumber 
                  min={8.0} 
                  max={20.0} 
                  step={0.5} 
                  placeholder="截断半径" 
                  style={{ width: '100%' }}
                  addonAfter="Å"
                />
              </Form.Item>
            </Col>
            <Col span={16}>
              <Form.Item
                name="additionalParams"
                label="额外参数 (JSON格式，可选)"
              >
                <Input.TextArea 
                  placeholder='{"pair_style": "lj/cut/coul/long", ...}' 
                  rows={2}
                />
              </Form.Item>
            </Col>
          </Row>
        </Card>
      </Card>
    );
  };

  // 渲染计算类型选择
  const renderCalculationSelection = () => {
    return (
      <Card 
        title={<Title level={4}>计算类型选择</Title>} 
        style={styles.contentCard}
      >
        <Paragraph style={{ marginBottom: 16 }}>
          请选择需要进行的计算类型，系统将根据选择生成相应的LAMMPS输入文件。
        </Paragraph>
        
        <Form.Item 
          name="calculationTypes" 
          label="计算类型"
          rules={[{ required: true, message: '请选择至少一种计算类型' }]}
          initialValue={['conductivity', 'diffusion', 'density_viscosity']}
        >
          <Checkbox.Group style={{ width: '100%' }}>
      <Row gutter={[24, 24]}>
              {calculationOptions.map(option => (
                <Col span={12} key={option.value}>
                  <Card 
                    hoverable 
                    style={{ height: '100%', transition: 'all 0.3s', borderRadius: 6 }}
                  >
              <div style={{ display: 'flex', alignItems: 'flex-start' }}>
                      <div style={{ ...styles.iconContainer, color: token.colorPrimary }}>
                        {option.icon}
                </div>
                <div>
                        <Checkbox value={option.value}>
                          <Title level={5} style={{ margin: 0 }}>{option.label}</Title>
                        </Checkbox>
                        <Paragraph style={{ marginTop: 8, color: token.colorTextSecondary }}>
                          {option.description}
                        </Paragraph>
                </div>
              </div>
            </Card>
          </Col>
        ))}
      </Row>
          </Checkbox.Group>
        </Form.Item>
        
        <Form.Item
          name="calculationName"
          label="计算名称"
          rules={[{ required: true, message: '请输入计算名称' }]}
          initialValue={`电解液计算${new Date().toLocaleDateString('zh-CN')}`}
        >
          <Input placeholder="请输入计算名称，用于区分不同的计算任务" />
        </Form.Item>
        
        <Form.Item
          name="calculationDescription"
          label="计算描述（可选）"
        >
          <Input.TextArea placeholder="请输入计算任务的描述信息" rows={2} />
        </Form.Item>
            </Card>
    );
  };

  // 添加登录处理函数
  const handleLogin = async (values: any) => {
    try {
      setLoginLoading(true);
      await login(values.username, values.password);
      message.success('登录成功');
      setIsLoginModalVisible(false);
      // 重新加载数据
      fetchSolventsAndSalts();
    } catch (error: any) {
      console.error('登录失败:', error);
      if (error.response && error.response.data) {
        message.error(`登录失败: ${JSON.stringify(error.response.data)}`);
      } else {
        message.error('登录失败，请检查用户名和密码');
      }
    } finally {
      setLoginLoading(false);
    }
  };

  // 添加加载配方列表的函数
  const loadFormulations = async () => {
    setFormulationsLoading(true);
    try {
      const response = await electrolyteService.getFormulations();
      console.log('获取配方列表成功:', response.data);
      setFormulations(response.data || []);
    } catch (error) {
      console.error('获取配方列表失败:', error);
      message.error('获取配方列表失败');
    } finally {
      setFormulationsLoading(false);
    }
  };

  // 主页面布局渲染
  return (
    <div className={styles.pageContainer as string}>
      {renderLogDrawer()}
      <Card className={styles.containerCard as string}>
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 16 }}>
          <Typography.Title level={2}>电解液配方计算工具</Typography.Title>
          <Button 
            type="default"
            icon={<BugOutlined />}
            onClick={() => setShowLogDrawer(true)}
          >
            显示处理日志
          </Button>
        </div>
        <Typography.Paragraph>
          设计电解液配方，进行分子动力学模拟，并获取电解液相关性能参数。
          {!isAuthenticated() && (
            <Button 
              type="link" 
              onClick={() => setIsLoginModalVisible(true)}
              icon={<LoginOutlined />}
            >
              点击登录
            </Button>
          )}
        </Typography.Paragraph>
        
        <Steps current={step} onChange={setStep}>
          <Step title="配方设计" description="设计电解液配方" icon={<SettingOutlined />} />
          <Step title="模拟条件" description="设置模拟参数" icon={<AppstoreOutlined />} />
          <Step title="计算提交" description="设置计算类型并提交" icon={<DownloadOutlined />} />
        </Steps>
        
        <Card className={styles.innerCard as string} style={{ marginTop: 24 }}>
          <Form
            form={form}
            layout="vertical"
            onFinish={handleSubmit}
            initialValues={{
              formulaName: `锂盐电解液${new Date().toLocaleDateString('zh-CN')}`,
              description: '',
              cations: [{ type: 'Li', concentration: 1.0 }],
              anions: [{ type: 'PF6', concentration: 1.0 }],
              solvents: [{ solvent: '1', volumeRatio: 1.0 }],
              temperature: 298,
              pressure: 1,
              timeStep: 1.0,
              equilibrationSteps: 1000000,
              productionSteps: 2000000,
              cutoff: 12.0,
              calculationTypes: ['conductivity', 'diffusion', 'density_viscosity'],
              calculationName: `电解液计算${new Date().toLocaleDateString('zh-CN')}`,
              calculationDescription: '',
            }}
          >
            {renderStepContent()}
            
            <div className={styles.stepButtons as string}>
              {step > 0 && (
                <Button 
                  icon={<LeftOutlined />} 
                  onClick={() => setStep(step - 1)}
                >
                  上一步
                </Button>
              )}
              
              {step < 2 && (
        <Button 
          type="primary" 
                  icon={<RightOutlined />} 
                  onClick={() => setStep(step + 1)}
        >
                  下一步
        </Button>
              )}
              
              {step === 2 && (
                <Button 
                  type="primary" 
                  onClick={() => form.submit()}
                  loading={submitting}
                >
                  保存配方
                </Button>
              )}
      </div>
          </Form>
        </Card>
      
        {/* 登录模态框 */}
        <Modal
          title="用户登录"
          open={isLoginModalVisible}
          onCancel={() => setIsLoginModalVisible(false)}
          footer={null}
        >
          <Form form={loginForm} onFinish={handleLogin} layout="vertical">
            <Form.Item
              name="username"
              label="用户名"
              rules={[{ required: true, message: '请输入用户名' }]}
            >
              <Input placeholder="请输入用户名" />
            </Form.Item>
            
            <Form.Item
              name="password"
              label="密码"
              rules={[{ required: true, message: '请输入密码' }]}
            >
              <Input.Password placeholder="请输入密码" />
            </Form.Item>
            
            <Form.Item>
              <Button type="primary" htmlType="submit" loading={loginLoading} block>
                登录
              </Button>
            </Form.Item>
          </Form>
        </Modal>
      </Card>
    </div>
  );
};

export default ElectrolyteCalculation; 