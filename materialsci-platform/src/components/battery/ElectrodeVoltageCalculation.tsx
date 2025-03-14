import React, { useState, useEffect } from 'react';
import { 
  Typography, 
  Card, 
  Steps, 
  Form, 
  Input, 
  Select, 
  Button, 
  Upload, 
  Radio, 
  Slider, 
  InputNumber, 
  Row, 
  Col, 
  Divider, 
  Space, 
  Alert,
  Spin,
  message,
  Checkbox,
  theme
} from 'antd';
import { 
  UploadOutlined, 
  InboxOutlined, 
  ThunderboltOutlined, 
  SettingOutlined, 
  LineChartOutlined, 
  CheckCircleOutlined,
  ExperimentOutlined,
  DatabaseOutlined,
  BuildOutlined
} from '@ant-design/icons';
import type { UploadProps } from 'antd';
import StructureViewer from '../StructureViewer';
import CrystalViewer from '../visualization/CrystalViewer';

const { Title, Paragraph, Text } = Typography;
const { Step } = Steps;
const { Option } = Select;
const { Dragger } = Upload;

const ElectrodeVoltageCalculation: React.FC = () => {
  const { token } = theme.useToken();
  const [currentStep, setCurrentStep] = useState(0);
  const [form] = Form.useForm();
  const [fileList, setFileList] = useState<any[]>([]);
  const [structurePreview, setStructurePreview] = useState<boolean>(false);
  const [loading, setLoading] = useState<boolean>(false);
  const [concentrationRange, setConcentrationRange] = useState<[number, number]>([0, 1]);
  const [concentrationSteps, setConcentrationSteps] = useState<number>(5);
  const [supercellMethod, setSupercellMethod] = useState<string>('auto');
  const [targetAtomCount, setTargetAtomCount] = useState<number>(100);
  const [supercellDimensions, setSupercellDimensions] = useState<{a: number, b: number, c: number}>({a: 1, b: 1, c: 1});
  const [modelsGenerated, setModelsGenerated] = useState<boolean>(false);
  const [width, setWidth] = useState<number>(window.innerWidth);
  
  // 监听窗口大小变化
  useEffect(() => {
    const handleResize = () => {
      setWidth(window.innerWidth);
    };
    
    window.addEventListener('resize', handleResize);
    return () => {
      window.removeEventListener('resize', handleResize);
    };
  }, []);

  // 文件上传配置
  const uploadProps: UploadProps = {
    name: 'file',
    multiple: false,
    fileList: fileList,
    action: 'https://run.mocky.io/v3/435e224c-44fb-4773-9faf-380c5e6a2188', // 使用一个模拟API
    beforeUpload: (file) => {
      const isCIF = file.name.endsWith('.cif') || 
                    file.name.endsWith('.poscar') || 
                    file.name.endsWith('.pdb') || 
                    file.name.endsWith('.xyz') || 
                    file.name.endsWith('.cml');
      if (!isCIF) {
        message.error('只能上传CIF、POSCAR、PDB、XYZ或CML格式的文件!');
      }
      
      // 清空旧的文件列表，确保只有一个文件被处理
      setFileList([]);
      setStructurePreview(false);
      form.resetFields(['structureInfo']);
      
      return isCIF || Upload.LIST_IGNORE;
    },
    onChange(info) {
      setFileList(info.fileList.slice(-1)); // 只保留最后上传的文件
      
      // 当文件上传成功时，解析结构
      if (info.file.status === 'done') {
        setLoading(true);
        // 实际项目中应该调用后端API解析文件内容
        // 这里模拟解析过程
        setTimeout(() => {
          // 基于文件类型和名称解析结构信息
          const fileName = info.file.name.toLowerCase();
          let structureData = {
            formula: '',
            spaceGroup: '',
            latticeParams: {
              a: 0,
              b: 0,
              c: 0,
              alpha: 0,
              beta: 0,
              gamma: 0
            },
            volume: 0,
            atomCount: 0,
            density: 0,
            elements: ['']
          };
          
          let atomPositions: any[] = [];
          
          if (fileName.includes('si') || fileName.includes('silicon')) {
            // 硅结构
            structureData = {
              formula: 'Si',
              spaceGroup: 'Fd-3m (227)',
              latticeParams: {
                a: 5.4307,
                b: 5.4307,
                c: 5.4307,
                alpha: 90,
                beta: 90,
                gamma: 90
              },
              volume: 160.1,
              atomCount: 8,
              density: 2.33,
              elements: ['Si (100%)']
            };
            
            // 硅原子位置 - 完整的金刚石结构
            atomPositions = [
              { element: 'Si', x: 0, y: 0, z: 0 },
              { element: 'Si', x: 2.7154, y: 0, z: 2.7154 },
              { element: 'Si', x: 0, y: 2.7154, z: 2.7154 },
              { element: 'Si', x: 2.7154, y: 2.7154, z: 0 },
              { element: 'Si', x: 1.3577, y: 1.3577, z: 1.3577 },
              { element: 'Si', x: 4.0731, y: 1.3577, z: 4.0731 },
              { element: 'Si', x: 1.3577, y: 4.0731, z: 4.0731 },
              { element: 'Si', x: 4.0731, y: 4.0731, z: 1.3577 }
            ];
          } else if (fileName.includes('licoo2') || fileName.includes('lico')) {
            // LiCoO2结构
            structureData = {
              formula: 'LiCoO₂',
              spaceGroup: 'R-3m (166)',
              latticeParams: {
                a: 2.8156,
                b: 2.8156,
                c: 14.0542,
                alpha: 90,
                beta: 90, 
                gamma: 120
              },
              volume: 96.48,
              atomCount: 12,
              density: 5.06,
              elements: ['Li (25%)', 'Co (25%)', 'O (50%)']
            };
            
            // LiCoO2原子位置 - 更完整的结构数据
            atomPositions = [
              { element: 'Li', x: 0, y: 0, z: 0, charge: 'Li+' },
              { element: 'Li', x: 1.4078, y: 1.4078, z: 0, charge: 'Li+' },
              { element: 'Li', x: 0, y: 1.4078, z: 7.0271, charge: 'Li+' },
              { element: 'Li', x: 1.4078, y: 0, z: 7.0271, charge: 'Li+' },
              
              { element: 'Co', x: 0, y: 0, z: 3.5135, charge: 'Co3+' },
              { element: 'Co', x: 1.4078, y: 1.4078, z: 3.5135, charge: 'Co3+' },
              { element: 'Co', x: 0, y: 1.4078, z: 10.5406, charge: 'Co3+' },
              { element: 'Co', x: 1.4078, y: 0, z: 10.5406, charge: 'Co3+' },
              
              { element: 'O', x: 0, y: 0, z: 1.6, charge: 'O2-' },
              { element: 'O', x: 1.4078, y: 1.4078, z: 1.6, charge: 'O2-' },
              { element: 'O', x: 0, y: 1.4078, z: 5.4, charge: 'O2-' },
              { element: 'O', x: 1.4078, y: 0, z: 5.4, charge: 'O2-' },
              { element: 'O', x: 0, y: 0, z: 8.6271, charge: 'O2-' },
              { element: 'O', x: 1.4078, y: 1.4078, z: 8.6271, charge: 'O2-' },
              { element: 'O', x: 0, y: 1.4078, z: 12.4542, charge: 'O2-' },
              { element: 'O', x: 1.4078, y: 0, z: 12.4542, charge: 'O2-' }
            ];
          } else {
            // 默认使用通用结构，实际情况应该解析文件内容
            structureData = {
              formula: fileName.split('.')[0],
              spaceGroup: '未知',
              latticeParams: {
                a: 5.0,
                b: 5.0,
                c: 5.0,
                alpha: 90,
                beta: 90,
                gamma: 90
              },
              volume: 125.0,
              atomCount: 8,
              density: 3.0,
              elements: ['未知']
            };
            
            // 默认原子位置（创建一个2x2x2的网格）
            atomPositions = [];
            for (let i = 0; i < 2; i++) {
              for (let j = 0; j < 2; j++) {
                for (let k = 0; k < 2; k++) {
                  atomPositions.push({
                    element: 'C', // 默认元素
                    x: i * 2.5,
                    y: j * 2.5,
                    z: k * 2.5
                  });
                }
              }
            }
          }
          
          // 保存解析结果
          form.setFieldsValue({
            structureInfo: structureData,
            atomPositions: atomPositions
          });
          
          setStructurePreview(true);
          setLoading(false);
          message.success(`${info.file.name} 文件解析成功`);
        }, 1500);
      } else if (info.file.status === 'error') {
        message.error(`${info.file.name} 文件上传失败`);
      }
    },
    onRemove: () => {
      setStructurePreview(false);
      form.resetFields(['structureInfo', 'atomPositions']);
      return true;
    },
    // 完全模拟上传过程，不发送实际请求
    customRequest: ({ onSuccess, onError, file }) => {
      setTimeout(() => {
        if (onSuccess) {
          onSuccess("ok" as any);
        }
      }, 1000);
    }
  };
  
  // 下一步
  const nextStep = () => {
    form.validateFields()
      .then(() => {
        if (currentStep === 2) {
          // 模拟模型生成过程
          setLoading(true);
          setTimeout(() => {
            setModelsGenerated(true);
            setLoading(false);
            setCurrentStep(currentStep + 1);
          }, 2000);
        } else {
          setCurrentStep(currentStep + 1);
        }
      })
      .catch(errorInfo => {
        console.log('表单验证失败:', errorInfo);
      });
  };
  
  // 上一步
  const prevStep = () => {
    setCurrentStep(currentStep - 1);
  };
  
  // 提交计算
  const submitCalculation = () => {
    setLoading(true);
    // 模拟提交过程
    setTimeout(() => {
      setLoading(false);
      setCurrentStep(4);
    }, 2000);
  };
  
  // 获取步骤图标
  const getStepIcon = (step: number) => {
    switch (step) {
      case 0: return <DatabaseOutlined />;
      case 1: return <ExperimentOutlined />;
      case 2: return <BuildOutlined />;
      case 3: return <ThunderboltOutlined />;
      case 4: return <SettingOutlined />;
      case 5: return <CheckCircleOutlined />;
      default: return null;
    }
  };

  // 渲染步骤内容
  const renderStepContent = () => {
    switch (currentStep) {
      case 0:
        return renderStructureUpload();
      case 1:
        return renderIonConfiguration();
      case 2:
        return renderSupercellSetup();
      case 3:
        return renderModelPreview();
      case 4:
        return renderCalculationSetup();
      case 5:
        return renderSubmission();
      default:
        return null;
    }
  };
  
  // 步骤1: 结构上传
  const renderStructureUpload = () => {
    // 获取结构信息
    const structureInfo = form.getFieldValue('structureInfo') || {
      formula: '',
      spaceGroup: '',
      latticeParams: {
        a: 0,
        b: 0,
        c: 0,
        alpha: 0,
        beta: 0,
        gamma: 0
      },
      volume: 0,
      atomCount: 0,
      density: 0,
      elements: ['']
    };
    
    // 获取原子位置数据
    const atomPositions = form.getFieldValue('atomPositions') || [];
    
    return (
      <div>
        <Card 
          title={
            <div style={{ display: 'flex', alignItems: 'center' }}>
              <DatabaseOutlined style={{ fontSize: 18, marginRight: 8, color: token.colorPrimary }} />
              <span>上传结构文件</span>
            </div>
          }
          className="step-card"
          bordered={false}
          style={{ 
            borderRadius: '8px', 
            boxShadow: '0 4px 12px rgba(0,0,0,0.05)',
            marginBottom: '24px'
          }}
        >
          <Paragraph style={{ fontSize: '16px', marginBottom: '20px' }}>
            请上传电极材料的晶体结构文件。支持的格式包括CIF、POSCAR、PDB、XYZ和CML。
          </Paragraph>
          
          <Form form={form} layout="vertical">
            <Form.Item 
              name="structureFile" 
              rules={[{ required: true, message: '请上传结构文件' }]}
            >
              <Dragger 
                {...uploadProps}
                style={{ 
                  background: token.colorBgContainer,
                  borderRadius: '8px',
                  borderColor: token.colorBorder,
                  padding: '10px'
                }}
              >
                <p className="ant-upload-drag-icon">
                  <InboxOutlined style={{ color: token.colorPrimary, fontSize: '48px' }} />
                </p>
                <p className="ant-upload-text" style={{ fontSize: '16px', fontWeight: 500 }}>点击或拖拽文件到此区域上传</p>
                <p className="ant-upload-hint" style={{ color: token.colorTextSecondary }}>
                  支持格式: .cif, .poscar, .pdb, .xyz, .cml
                </p>
              </Dragger>
            </Form.Item>
            
            <Form.Item name="structureInfo" hidden>
              <Input />
            </Form.Item>
            
            <Form.Item name="atomPositions" hidden>
              <Input />
            </Form.Item>
          </Form>
          
          {loading && (
            <div style={{ 
              textAlign: 'center', 
              margin: '24px 0', 
              padding: '30px', 
              background: token.colorBgContainerDisabled,
              borderRadius: '8px' 
            }}>
              <Spin size="large" tip="正在解析结构文件..." />
            </div>
          )}
        </Card>
        
        {structurePreview && (
          <>
            <Card 
              title={
                <div style={{ display: 'flex', alignItems: 'center' }}>
                  <ThunderboltOutlined style={{ fontSize: 18, marginRight: 8, color: '#52c41a' }} />
                  <span>结构信息</span>
                </div>
              }
              className="info-card"
              style={{ 
                borderRadius: '8px', 
                boxShadow: '0 4px 12px rgba(0,0,0,0.05)',
                marginBottom: '24px',
                borderLeft: `4px solid ${token.colorSuccess}`
              }}
            >
              <Row gutter={24}>
                <Col xs={24} md={12}>
                  <div style={{ 
                    background: token.colorBgContainer, 
                    padding: '16px',
                    borderRadius: '8px',
                    border: `1px solid ${token.colorBorderSecondary}`,
                    height: '100%'
                  }}>
                    <p style={{ fontSize: '16px', margin: '8px 0' }}>
                      <Text strong style={{ color: token.colorTextHeading }}>化学式:</Text> 
                      <Text style={{ marginLeft: '8px', fontSize: '16px' }}>{structureInfo.formula}</Text>
                    </p>
                    <p style={{ fontSize: '16px', margin: '8px 0' }}>
                      <Text strong style={{ color: token.colorTextHeading }}>空间群:</Text> 
                      <Text style={{ marginLeft: '8px' }}>{structureInfo.spaceGroup}</Text>
                    </p>
                    <p style={{ fontSize: '16px', margin: '8px 0' }}>
                      <Text strong style={{ color: token.colorTextHeading }}>晶胞参数:</Text>
                    </p>
                    <div style={{ 
                      padding: '8px 16px', 
                      background: token.colorFillQuaternary,
                      borderRadius: '6px',
                      margin: '8px 0 12px'
                    }}>
                      <Row gutter={[16, 8]}>
                        <Col span={8}>a = {structureInfo.latticeParams.a} Å</Col>
                        <Col span={8}>b = {structureInfo.latticeParams.b} Å</Col>
                        <Col span={8}>c = {structureInfo.latticeParams.c} Å</Col>
                        <Col span={8}>α = {structureInfo.latticeParams.alpha}°</Col>
                        <Col span={8}>β = {structureInfo.latticeParams.beta}°</Col>
                        <Col span={8}>γ = {structureInfo.latticeParams.gamma}°</Col>
                      </Row>
                    </div>
                    <p style={{ fontSize: '16px', margin: '8px 0' }}>
                      <Text strong style={{ color: token.colorTextHeading }}>体积:</Text> 
                      <Text style={{ marginLeft: '8px' }}>{structureInfo.volume} Å³</Text>
                    </p>
                    <p style={{ fontSize: '16px', margin: '8px 0' }}>
                      <Text strong style={{ color: token.colorTextHeading }}>原子数:</Text> 
                      <Text style={{ marginLeft: '8px' }}>{structureInfo.atomCount}</Text>
                    </p>
                    <p style={{ fontSize: '16px', margin: '8px 0' }}>
                      <Text strong style={{ color: token.colorTextHeading }}>密度:</Text> 
                      <Text style={{ marginLeft: '8px' }}>{structureInfo.density} g/cm³</Text>
                    </p>
                    <p style={{ fontSize: '16px', margin: '8px 0' }}>
                      <Text strong style={{ color: token.colorTextHeading }}>元素组成:</Text> 
                      <Text style={{ marginLeft: '8px' }}>{structureInfo.elements.join(', ')}</Text>
                    </p>
                  </div>
                </Col>
                <Col xs={24} md={12}>
                  <div style={{ 
                    background: token.colorBgContainer,
                    padding: '16px',
                    borderRadius: '8px',
                    height: '100%',
                    display: 'flex',
                    flexDirection: 'column',
                    justifyContent: 'center',
                    alignItems: 'center',
                    border: `1px solid ${token.colorBorderSecondary}`
                  }}>
                    <p style={{ 
                      fontSize: '16px', 
                      fontWeight: 500, 
                      color: token.colorTextHeading,
                      marginBottom: '16px'
                    }}>
                      <BuildOutlined style={{ marginRight: '8px' }} />
                      三维结构预览
                    </p>
                    {/* 使用新的晶体结构查看器 */}
                    <CrystalViewer 
                      atoms={atomPositions} 
                      width={width > 1200 ? 450 : 350} 
                      height={400} 
                      style="polyhedra"
                      showUnitCell={true}
                      showLabels={true}
                      backgroundColor="#f9f9f9"
                      rotationSpeed={0.5}
                      latticeParams={structureInfo.latticeParams}
                    />
                  </div>
                </Col>
              </Row>
            </Card>
            
            <Card
              title={
                <div style={{ display: 'flex', alignItems: 'center' }}>
                  <BuildOutlined style={{ fontSize: 18, marginRight: 8, color: token.colorPrimary }} />
                  <span>高级结构可视化</span>
                </div>
              }
              bordered={false}
              style={{ 
                borderRadius: '8px', 
                boxShadow: '0 4px 12px rgba(0,0,0,0.05)',
                marginBottom: '24px'
              }}
            >
              <div style={{ width: '100%', display: 'flex', justifyContent: 'center' }}>
                <CrystalViewer 
                  atoms={atomPositions} 
                  width={width > 1200 ? 1100 : Math.min(900, width - 50)} 
                  height={700} 
                  style="polyhedra"
                  colorScheme="element"
                  showUnitCell={true}
                  showLabels={true}
                  backgroundColor="#f5f5f5"
                  rotationSpeed={0.3}
                  latticeParams={structureInfo.latticeParams}
                />
              </div>
            </Card>
          </>
        )}
        
        <div style={{ 
          display: 'flex', 
          justifyContent: 'flex-end', 
          marginTop: '16px' 
        }}>
          <Button 
            type="primary" 
            onClick={nextStep} 
            disabled={!structurePreview}
            style={{ 
              height: '40px', 
              borderRadius: '6px',
              padding: '0 24px',
              fontSize: '16px'
            }}
          >
            下一步
          </Button>
        </div>
      </div>
    );
  };
  
  // 步骤2: 离子配置
  const renderIonConfiguration = () => {
    return (
      <div>
        <Card 
          title={
            <div style={{ display: 'flex', alignItems: 'center' }}>
              <ExperimentOutlined style={{ fontSize: 18, marginRight: 8, color: token.colorPrimary }} />
              <span>离子与浓度设置</span>
            </div>
          }
          className="step-card"
          bordered={false}
          style={{ 
            borderRadius: '8px', 
            boxShadow: '0 4px 12px rgba(0,0,0,0.05)',
            marginBottom: '24px'
          }}
        >
          <Paragraph style={{ fontSize: '16px', marginBottom: '20px' }}>
            选择工作离子类型并设置浓度范围，系统将在指定范围内生成多个不同浓度的计算模型。
          </Paragraph>
          
          <Form form={form} layout="vertical">
            <Form.Item 
              name="ionType" 
              label={<Text strong style={{ fontSize: '16px' }}>工作离子类型</Text>}
              initialValue="Li" 
              rules={[{ required: true, message: '请选择工作离子类型' }]}
            >
              <Radio.Group size="large" style={{ marginBottom: '16px' }}>
                <Radio.Button value="Li" style={{ borderRadius: '4px', marginRight: '8px', marginBottom: '8px' }}>
                  <span style={{ padding: '0 8px' }}>Li⁺</span>
                </Radio.Button>
                <Radio.Button value="Na" style={{ borderRadius: '4px', marginRight: '8px', marginBottom: '8px' }}>
                  <span style={{ padding: '0 8px' }}>Na⁺</span>
                </Radio.Button>
                <Radio.Button value="K" style={{ borderRadius: '4px', marginRight: '8px', marginBottom: '8px' }}>
                  <span style={{ padding: '0 8px' }}>K⁺</span>
                </Radio.Button>
                <Radio.Button value="Mg" style={{ borderRadius: '4px', marginRight: '8px', marginBottom: '8px' }}>
                  <span style={{ padding: '0 8px' }}>Mg²⁺</span>
                </Radio.Button>
                <Radio.Button value="Ca" style={{ borderRadius: '4px', marginRight: '8px', marginBottom: '8px' }}>
                  <span style={{ padding: '0 8px' }}>Ca²⁺</span>
                </Radio.Button>
                <Radio.Button value="Al" style={{ borderRadius: '4px', marginRight: '8px', marginBottom: '8px' }}>
                  <span style={{ padding: '0 8px' }}>Al³⁺</span>
                </Radio.Button>
              </Radio.Group>
            </Form.Item>
            
            <Form.Item 
              label={<Text strong style={{ fontSize: '16px' }}>离子插入比例范围</Text>}
              required
            >
              <div style={{ 
                background: token.colorBgContainer, 
                padding: '16px',
                borderRadius: '8px',
                border: `1px solid ${token.colorBorderSecondary}`,
                marginBottom: '16px'
              }}>
                <Row gutter={16} align="middle">
                  <Col span={16}>
                    <Slider
                      range
                      min={0}
                      max={1}
                      step={0.01}
                      value={concentrationRange}
                      onChange={(value) => setConcentrationRange(value as [number, number])}
                      tooltip={{ formatter: (value) => `${(value as number * 100).toFixed(0)}%` }}
                      styles={{ track: { background: token.colorPrimary } }}
                    />
                  </Col>
                  <Col span={8}>
                    <Space align="center">
                      <InputNumber
                        min={0}
                        max={1}
                        step={0.01}
                        value={concentrationRange[0]}
                        onChange={(value) => setConcentrationRange([value || 0, concentrationRange[1]])}
                        style={{ width: '80px' }}
                        addonAfter="%"
                        formatter={(value) => `${(Number(value) * 100).toFixed(0)}`}
                        parser={(value) => (Number(value) || 0) / 100}
                      />
                      <Text>至</Text>
                      <InputNumber
                        min={0}
                        max={1}
                        step={0.01}
                        value={concentrationRange[1]}
                        onChange={(value) => setConcentrationRange([concentrationRange[0], value || 1])}
                        style={{ width: '80px' }}
                        addonAfter="%"
                        formatter={(value) => `${(Number(value) * 100).toFixed(0)}`}
                        parser={(value) => (Number(value) || 0) / 100}
                      />
                    </Space>
                  </Col>
                </Row>
              </div>
            </Form.Item>
            
            <Form.Item 
              label={<Text strong style={{ fontSize: '16px' }}>浓度步数</Text>}
              required
            >
              <div style={{ 
                background: token.colorBgContainer, 
                padding: '16px',
                borderRadius: '8px',
                border: `1px solid ${token.colorBorderSecondary}`
              }}>
                <Row gutter={16} align="middle">
                  <Col span={16}>
                    <Slider
                      min={2}
                      max={10}
                      value={concentrationSteps}
                      onChange={(value) => setConcentrationSteps(value)}
                      marks={{
                        2: '2',
                        5: '5',
                        10: '10'
                      }}
                      styles={{ track: { background: token.colorPrimary } }}
                    />
                  </Col>
                  <Col span={8}>
                    <InputNumber
                      min={2}
                      max={10}
                      value={concentrationSteps}
                      onChange={(value) => setConcentrationSteps(value || 5)}
                      style={{ width: '80px' }}
                      addonAfter="个"
                    />
                  </Col>
                </Row>
              </div>
            </Form.Item>
          </Form>
          
          <Alert
            message="浓度设置说明"
            description={
              <div style={{ fontSize: '14px', lineHeight: '1.6' }}>
                系统将在 <Text strong>{(concentrationRange[0] * 100).toFixed(0)}%</Text> 到 <Text strong>{(concentrationRange[1] * 100).toFixed(0)}%</Text> 的范围内生成 <Text strong>{concentrationSteps}</Text> 个不同浓度的模型进行计算。
                <div style={{ marginTop: '8px' }}>
                  这将产生更全面的电极电位曲线，帮助您了解材料在不同离子浓度下的性能变化。
                </div>
              </div>
            }
            type="info"
            showIcon
            style={{ 
              marginBottom: 16,
              borderRadius: '8px',
              borderLeft: `4px solid ${token.colorInfo}`
            }}
          />
        </Card>
        
        <div style={{ 
          display: 'flex', 
          justifyContent: 'space-between', 
          marginTop: '16px' 
        }}>
          <Button 
            onClick={prevStep} 
            style={{ 
              height: '40px', 
              borderRadius: '6px',
              padding: '0 24px',
              fontSize: '16px'
            }}
          >
            上一步
          </Button>
          <Button 
            type="primary" 
            onClick={nextStep}
            style={{ 
              height: '40px', 
              borderRadius: '6px',
              padding: '0 24px',
              fontSize: '16px'
            }}
          >
            下一步
          </Button>
        </div>
      </div>
    );
  };
  
  // 步骤3: 超胞设置
  const renderSupercellSetup = () => {
    return (
      <div>
        <Card 
          title={
            <div style={{ display: 'flex', alignItems: 'center' }}>
              <BuildOutlined style={{ fontSize: 18, marginRight: 8, color: token.colorPrimary }} />
              <span>超胞设置</span>
            </div>
          }
          className="step-card"
          bordered={false}
          style={{ 
            borderRadius: '8px', 
            boxShadow: '0 4px 12px rgba(0,0,0,0.05)',
            marginBottom: '24px'
          }}
        >
          <Paragraph style={{ fontSize: '16px', marginBottom: '20px' }}>
            设置超胞大小，以确保计算模型包含足够的原子，提高计算精度和结果可靠性。
          </Paragraph>
          
          <Form form={form} layout="vertical">
            <Form.Item 
              name="supercellMethod" 
              label={<Text strong style={{ fontSize: '16px' }}>超胞设置方式</Text>}
              initialValue="auto"
            >
              <Radio.Group 
                onChange={(e) => setSupercellMethod(e.target.value)}
                size="large"
                style={{ marginBottom: '16px' }}
                buttonStyle="solid"
              >
                <Radio.Button value="auto" style={{ borderRadius: '4px 0 0 4px' }}>
                  <span style={{ padding: '0 16px' }}>自动规划</span>
                </Radio.Button>
                <Radio.Button value="manual" style={{ borderRadius: '0 4px 4px 0' }}>
                  <span style={{ padding: '0 16px' }}>手动设置</span>
                </Radio.Button>
              </Radio.Group>
            </Form.Item>
            
            <div style={{ 
              background: token.colorBgContainer, 
              padding: '24px',
              borderRadius: '8px',
              border: `1px solid ${token.colorBorderSecondary}`,
              marginBottom: '24px'
            }}>
              {supercellMethod === 'auto' ? (
                <Form.Item 
                  label={<Text strong style={{ fontSize: '16px' }}>目标原子数</Text>}
                  required
                >
                  <Row gutter={16} align="middle">
                    <Col span={16}>
                      <Slider
                        min={50}
                        max={300}
                        step={10}
                        value={targetAtomCount}
                        onChange={(value) => setTargetAtomCount(value)}
                        marks={{
                          50: '50',
                          100: '100',
                          200: '200',
                          300: '300'
                        }}
                        styles={{ track: { background: token.colorPrimary } }}
                      />
                    </Col>
                    <Col span={8}>
                      <InputNumber
                        min={50}
                        max={300}
                        step={10}
                        value={targetAtomCount}
                        onChange={(value) => setTargetAtomCount(value || 100)}
                        style={{ width: '100px' }}
                        addonAfter="原子"
                      />
                    </Col>
                  </Row>
                  <div style={{ 
                    marginTop: '16px', 
                    color: token.colorTextSecondary, 
                    background: token.colorFillQuaternary,
                    padding: '12px',
                    borderRadius: '6px',
                    fontSize: '14px'
                  }}>
                    系统将自动计算合适的超胞倍数，使最终模型包含约 <Text strong>{targetAtomCount}</Text> 个原子。
                    较多的原子数可提高计算精度，但会增加计算时间。
                  </div>
                </Form.Item>
              ) : (
                <Form.Item 
                  label={<Text strong style={{ fontSize: '16px' }}>超胞尺寸</Text>}
                  required
                >
                  <Row gutter={16} align="middle">
                    <Col span={24} style={{ marginBottom: '16px' }}>
                      <Text style={{ marginRight: '16px' }}>沿晶格矢量的重复次数:</Text>
                    </Col>
                    <Col span={8}>
                      <Form.Item label="a 方向" style={{ marginBottom: '8px' }}>
                        <InputNumber
                          min={1}
                          max={10}
                          value={supercellDimensions.a}
                          onChange={(value) => setSupercellDimensions({...supercellDimensions, a: value || 1})}
                          style={{ width: '100%' }}
                        />
                      </Form.Item>
                    </Col>
                    <Col span={8}>
                      <Form.Item label="b 方向" style={{ marginBottom: '8px' }}>
                        <InputNumber
                          min={1}
                          max={10}
                          value={supercellDimensions.b}
                          onChange={(value) => setSupercellDimensions({...supercellDimensions, b: value || 1})}
                          style={{ width: '100%' }}
                        />
                      </Form.Item>
                    </Col>
                    <Col span={8}>
                      <Form.Item label="c 方向" style={{ marginBottom: '8px' }}>
                        <InputNumber
                          min={1}
                          max={10}
                          value={supercellDimensions.c}
                          onChange={(value) => setSupercellDimensions({...supercellDimensions, c: value || 1})}
                          style={{ width: '100%' }}
                        />
                      </Form.Item>
                    </Col>
                  </Row>
                  <div style={{ 
                    marginTop: '16px', 
                    color: token.colorTextSecondary, 
                    background: token.colorFillQuaternary,
                    padding: '12px',
                    borderRadius: '6px',
                    fontSize: '14px'
                  }}>
                    当前设置将产生 <Text strong>{supercellDimensions.a * supercellDimensions.b * supercellDimensions.c}</Text> 倍大小的超胞，
                    包含约 <Text strong>{supercellDimensions.a * supercellDimensions.b * supercellDimensions.c * 8}</Text> 个原子。
                  </div>
                </Form.Item>
              )}
            </div>
          </Form>
          
          <Alert
            message="超胞设置建议"
            description={
              <div style={{ fontSize: '14px', lineHeight: '1.6' }}>
                <ul style={{ paddingLeft: '20px', margin: '8px 0' }}>
                  <li>对于初步筛选计算，建议使用较小的超胞（50-100个原子）</li>
                  <li>对于精确的电极电位计算，建议使用中等大小的超胞（100-200个原子）</li>
                  <li>对于深入的缺陷和界面研究，建议使用较大的超胞（200-300个原子）</li>
                </ul>
              </div>
            }
            type="info"
            showIcon
            style={{ 
              marginBottom: 16,
              borderRadius: '8px',
              borderLeft: `4px solid ${token.colorInfo}`
            }}
          />
        </Card>
        
        <div style={{ 
          display: 'flex', 
          justifyContent: 'space-between', 
          marginTop: '16px' 
        }}>
          <Button 
            onClick={prevStep} 
            style={{ 
              height: '40px', 
              borderRadius: '6px',
              padding: '0 24px',
              fontSize: '16px'
            }}
          >
            上一步
          </Button>
          <Button 
            type="primary" 
            onClick={nextStep}
            style={{ 
              height: '40px', 
              borderRadius: '6px',
              padding: '0 24px',
              fontSize: '16px'
            }}
          >
            下一步
          </Button>
        </div>
      </div>
    );
  };
  
  // 步骤4: 模型预览
  const renderModelPreview = () => {
    // 简单模拟多个原子位置数据
    const demoModelAtoms = [
      // ... 这里可以放更多的原子数据来模拟超胞
      // 为简洁起见，这里沿用了之前的原子数据
    ];
    
    return (
      <div>
        <Card 
          title={
            <div style={{ display: 'flex', alignItems: 'center' }}>
              <ThunderboltOutlined style={{ fontSize: 18, marginRight: 8, color: token.colorPrimary }} />
              <span>模型预览</span>
            </div>
          }
          className="step-card"
          bordered={false}
          style={{ 
            borderRadius: '8px', 
            boxShadow: '0 4px 12px rgba(0,0,0,0.05)',
            marginBottom: '24px'
          }}
        >
          {loading ? (
            <div style={{ 
              textAlign: 'center', 
              padding: '60px 0',
              background: token.colorBgContainerDisabled,
              borderRadius: '8px'
            }}>
              <Spin size="large" tip="正在生成计算模型..." />
              <p style={{ marginTop: '16px', color: token.colorTextSecondary }}>
                这可能需要一些时间，请耐心等待...
              </p>
            </div>
          ) : (
            <>
              <Paragraph style={{ fontSize: '16px', marginBottom: '20px' }}>
                计算模型已生成完毕。系统将根据您的设置创建不同离子浓度的超胞模型。
              </Paragraph>
              
              <Alert
                message="已生成的计算模型"
                description={
                  <div style={{ fontSize: '14px', lineHeight: '1.6' }}>
                    <p>根据您的设置，系统已生成以下计算模型：</p>
                    <ul style={{ paddingLeft: '20px', margin: '8px 0' }}>
                      <li>总计生成 <Text strong>{concentrationSteps}</Text> 个不同离子浓度的模型</li>
                      <li>浓度范围: <Text strong>{(concentrationRange[0] * 100).toFixed(0)}% - {(concentrationRange[1] * 100).toFixed(0)}%</Text></li>
                      <li>每个模型包含约 <Text strong>{supercellMethod === 'auto' ? targetAtomCount : supercellDimensions.a * supercellDimensions.b * supercellDimensions.c * 8}</Text> 个原子</li>
                    </ul>
                  </div>
                }
                type="success"
                showIcon
                style={{ 
                  marginBottom: 24,
                  borderRadius: '8px',
                  borderLeft: `4px solid ${token.colorSuccess}`
                }}
              />
              
              <div style={{ marginBottom: '24px' }}>
                <Row gutter={[16, 16]}>
                  <Col xs={24} lg={12}>
                    <Card 
                      title="模型结构预览" 
                      bordered
                      style={{ borderRadius: '8px' }}
                    >
                      <div style={{ 
                        height: '300px', 
                        display: 'flex',
                        justifyContent: 'center',
                        alignItems: 'center',
                        background: token.colorBgContainer
                      }}>
                        <StructureViewer 
                          atoms={fileList[0]?.name.startsWith('Si') ? [
                            { element: 'Si', x: 0, y: 0, z: 0 },
                            { element: 'Si', x: 2.7154, y: 0, z: 2.7154 },
                            { element: 'Si', x: 0, y: 2.7154, z: 2.7154 },
                            { element: 'Si', x: 2.7154, y: 2.7154, z: 0 },
                            { element: 'Si', x: 1.3577, y: 1.3577, z: 1.3577 },
                            { element: 'Si', x: 4.0731, y: 1.3577, z: 4.0731 },
                            { element: 'Si', x: 1.3577, y: 4.0731, z: 4.0731 },
                            { element: 'Si', x: 4.0731, y: 4.0731, z: 1.3577 }
                          ] : []} 
                          width={400} 
                          height={250} 
                          showBonds={true}
                          showUnitCell={true}
                          showLabels={true}
                          quality="high"
                          backgroundColor="#f9f9f9"
                          rotationSpeed={0.01}
                        />
                      </div>
                    </Card>
                  </Col>
                  <Col xs={24} lg={12}>
                    <Card 
                      title="预计计算资源" 
                      bordered
                      style={{ borderRadius: '8px' }}
                    >
                      <div style={{ padding: '16px' }}>
                        <p style={{ margin: '12px 0', fontSize: '15px' }}>
                          <Text strong>计算类型:</Text> 
                          <Text style={{ marginLeft: '8px' }}>电极电位计算</Text>
                        </p>
                        <p style={{ margin: '12px 0', fontSize: '15px' }}>
                          <Text strong>预计计算时长:</Text> 
                          <Text style={{ marginLeft: '8px' }}>约 2-3 小时</Text>
                        </p>
                        <p style={{ margin: '12px 0', fontSize: '15px' }}>
                          <Text strong>CPU核心数:</Text> 
                          <Text style={{ marginLeft: '8px' }}>8</Text>
                        </p>
                        <p style={{ margin: '12px 0', fontSize: '15px' }}>
                          <Text strong>内存需求:</Text> 
                          <Text style={{ marginLeft: '8px' }}>约 4GB</Text>
                        </p>
                        <p style={{ margin: '12px 0', fontSize: '15px' }}>
                          <Text strong>磁盘空间:</Text> 
                          <Text style={{ marginLeft: '8px' }}>约 500MB</Text>
                        </p>
                        <p style={{ margin: '12px 0', fontSize: '15px' }}>
                          <Text strong>计算优先级:</Text> 
                          <Text style={{ marginLeft: '8px' }}>正常</Text>
                        </p>
                      </div>
                    </Card>
                  </Col>
                </Row>
              </div>
            </>
          )}
        </Card>
        
        <div style={{ 
          display: 'flex', 
          justifyContent: 'space-between', 
          marginTop: '16px' 
        }}>
          <Button 
            onClick={prevStep} 
            style={{ 
              height: '40px', 
              borderRadius: '6px',
              padding: '0 24px',
              fontSize: '16px'
            }}
            disabled={loading}
          >
            上一步
          </Button>
          <Button 
            type="primary" 
            onClick={nextStep}
            style={{ 
              height: '40px', 
              borderRadius: '6px',
              padding: '0 24px',
              fontSize: '16px'
            }}
            disabled={loading || !modelsGenerated}
          >
            下一步
          </Button>
        </div>
      </div>
    );
  };
  
  // 步骤5: 计算参数设置
  const renderCalculationSetup = () => {
    return (
      <div>
        <Title level={4}>计算参数设置</Title>
        <Paragraph>
          设置计算参数和资源分配。
        </Paragraph>
        
        <Form form={form} layout="vertical">
          <Form.Item 
            name="calculationName" 
            label="计算名称" 
            initialValue="电极电位计算-LiCoO₂"
            rules={[{ required: true, message: '请输入计算名称' }]}
          >
            <Input placeholder="请输入计算任务名称" />
          </Form.Item>
          
          <Form.Item 
            name="calculationPrecision" 
            label="计算精度" 
            initialValue="standard"
          >
            <Radio.Group>
              <Radio value="fast">快速测试</Radio>
              <Radio value="standard">标准精度</Radio>
              <Radio value="high">高精度</Radio>
              <Radio value="custom">自定义</Radio>
            </Radio.Group>
          </Form.Item>
          
          <Row gutter={24}>
            <Col span={12}>
              <Form.Item 
                name="encut" 
                label="平面波截断能 (ENCUT)" 
                initialValue={400}
              >
                <InputNumber min={200} max={800} addonAfter="eV" style={{ width: '100%' }} />
              </Form.Item>
            </Col>
            <Col span={12}>
              <Form.Item 
                name="ediff" 
                label="电子收敛标准 (EDIFF)" 
                initialValue="1E-5"
              >
                <Select>
                  <Option value="1E-4">1E-4 (快速)</Option>
                  <Option value="1E-5">1E-5 (标准)</Option>
                  <Option value="1E-6">1E-6 (精确)</Option>
                </Select>
              </Form.Item>
            </Col>
          </Row>
          
          <Row gutter={24}>
            <Col span={12}>
              <Form.Item 
                name="ismear" 
                label="展宽方法 (ISMEAR)" 
                initialValue={1}
              >
                <Select>
                  <Option value={-5}>-5 (四面体方法)</Option>
                  <Option value={0}>0 (高斯展宽)</Option>
                  <Option value={1}>1 (Methfessel-Paxton)</Option>
                </Select>
              </Form.Item>
            </Col>
            <Col span={12}>
              <Form.Item 
                name="kpoints" 
                label="K点网格" 
                initialValue="3x3x3"
              >
                <Input placeholder="例如: 3x3x3" />
              </Form.Item>
            </Col>
          </Row>
          
          <Divider orientation="left">计算资源</Divider>
          
          <Row gutter={24}>
            <Col span={12}>
              <Form.Item 
                name="cpuCores" 
                label="CPU核心数" 
                initialValue={8}
              >
                <Select>
                  <Option value={4}>4</Option>
                  <Option value={8}>8</Option>
                  <Option value={16}>16</Option>
                  <Option value={32}>32</Option>
                </Select>
              </Form.Item>
            </Col>
            <Col span={12}>
              <Form.Item 
                name="memory" 
                label="内存" 
                initialValue={16}
              >
                <Select>
                  <Option value={8}>8 GB</Option>
                  <Option value={16}>16 GB</Option>
                  <Option value={32}>32 GB</Option>
                  <Option value={64}>64 GB</Option>
                </Select>
              </Form.Item>
            </Col>
          </Row>
          
          <Form.Item 
            name="estimatedTime" 
            label="预计运行时间" 
            initialValue={2}
          >
            <InputNumber min={1} max={48} addonAfter="小时" style={{ width: '100%' }} />
          </Form.Item>
          
          <Form.Item 
            name="notifications" 
            label="计算完成通知" 
            initialValue={['email']}
          >
            <Checkbox.Group>
              <Checkbox value="email">邮件通知</Checkbox>
              <Checkbox value="platform">平台消息</Checkbox>
            </Checkbox.Group>
          </Form.Item>
        </Form>
        
        <Divider />
        
        <Space>
          <Button onClick={prevStep}>上一步</Button>
          <Button type="primary" onClick={submitCalculation}>提交计算</Button>
        </Space>
      </div>
    );
  };
  
  // 步骤6: 提交确认
  const renderSubmission = () => {
    return (
      <div style={{ textAlign: 'center', padding: '20px 0' }}>
        <CheckCircleOutlined style={{ fontSize: 72, color: '#52c41a', marginBottom: 24 }} />
        <Title level={3}>计算任务已成功提交！</Title>
        <Paragraph>
          您的电极电位计算任务已加入计算队列，系统将自动处理您的计算请求。您可以在用户中心的"我的计算任务"中查看任务进度。
        </Paragraph>
        <div style={{ marginTop: 24 }}>
          <Space size="large">
            <Button type="primary" href="/user/tasks">
              查看我的任务
            </Button>
            <Button href="/battery/electrode-material">
              返回电极材料页面
            </Button>
          </Space>
        </div>
      </div>
    );
  };
  
  return (
    <div style={{ padding: '24px', background: token.colorBgLayout, minHeight: '100vh' }}>
      <Card
        title={
          <Title level={3} style={{ margin: 0, color: token.colorTextHeading }}>
            电极电位计算
          </Title>
        }
        bordered={false}
        style={{ 
          borderRadius: '12px', 
          boxShadow: '0 2px 8px rgba(0,0,0,0.08)',
          marginBottom: '24px'
        }}
      >
        <Steps
          current={currentStep}
          labelPlacement="vertical"
          style={{ 
            padding: '16px 0 32px',
            borderBottom: `1px solid ${token.colorBorderSecondary}`,
            marginBottom: '24px'
          }}
          items={[
            {
              title: '结构上传',
              icon: getStepIcon(0)
            },
            {
              title: '离子设置',
              icon: getStepIcon(1)
            },
            {
              title: '超胞设置',
              icon: getStepIcon(2)
            },
            {
              title: '模型预览',
              icon: getStepIcon(3)
            },
            {
              title: '计算设置',
              icon: getStepIcon(4)
            },
            {
              title: '任务提交',
              icon: getStepIcon(5)
            }
          ]}
        />
        <div style={{ padding: '16px 0' }}>
          {renderStepContent()}
        </div>
      </Card>
    </div>
  );
};

export default ElectrodeVoltageCalculation; 