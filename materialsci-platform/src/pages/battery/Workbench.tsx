import React, { useState, useEffect } from 'react';
import { 
  Typography, 
  Card, 
  Row, 
  Col, 
  Divider, 
  Tabs, 
  Form, 
  Input, 
  Select, 
  Button, 
  Upload, 
  Checkbox, 
  Space,
  Steps,
  Radio,
  Alert,
  Spin,
  InputNumber,
  message
} from 'antd';
import { 
  UploadOutlined, 
  ExperimentOutlined, 
  SettingOutlined, 
  FileTextOutlined,
  CheckCircleOutlined,
  InfoCircleOutlined,
  InboxOutlined,
  LineChartOutlined
} from '@ant-design/icons';
import { useParams, useNavigate, useLocation } from 'react-router-dom';
import ElectrodeVoltageCalculation from '../../components/battery/ElectrodeVoltageCalculation';
import FormationEnergyCalculation from '../../components/battery/FormationEnergyCalculation';
import IonMigrationBarrierCalculation from '../../components/battery/IonMigrationBarrierCalculation';
import ElectronicStructureCalculation from '../../components/battery/ElectronicStructureCalculation';

const { Title, Paragraph, Text } = Typography;
const { TabPane } = Tabs;
const { Option } = Select;
const { Step } = Steps;

interface WorkbenchProps {
  calculationType?: string;
}

const Workbench: React.FC<WorkbenchProps> = ({ calculationType }) => {
  const navigate = useNavigate();
  const location = useLocation();
  const [currentStep, setCurrentStep] = useState(0);
  const [selectedCalculations, setSelectedCalculations] = useState<string[]>([]);
  const [form] = Form.useForm();
  const [structurePreview, setStructurePreview] = useState(false);
  const [loading, setLoading] = useState(false);
  const [showCalculation, setShowCalculation] = useState(false);
  
  // 使用传入的calculationType或从URL参数中获取
  const params = useParams();
  const moduleType = calculationType || params.type || 'electrode';
  
  // 从location state中获取计算类型
  useEffect(() => {
    const state = location.state as { calculationType?: string };
    if (state && state.calculationType) {
      setSelectedCalculations([state.calculationType]);
      setShowCalculation(true);
    }
  }, [location]);
  
  // 获取标题
  const getTitle = () => {
    switch(moduleType) {
      case 'electrode':
        return '电极材料计算工作台';
      case 'electrolyte':
        return '电解液计算工作台';
      case 'fullBattery':
        return '全电池系统计算工作台';
      default:
        return '计算工作台';
    }
  };
  
  // 获取可选计算类型
  const getCalculationOptions = () => {
    switch(moduleType) {
      case 'electrode':
        return [
          { value: 'voltage', label: '电极电位计算' },
          { value: 'barrier', label: '离子迁移能垒计算' },
          { value: 'formation', label: '形成能计算' },
          { value: 'electronic', label: '电子结构计算' },
        ];
      case 'electrolyte':
        return [
          { value: 'molecule', label: '单分子性质计算' },
          { value: 'solvation', label: '离子溶剂化结构计算' },
          { value: 'transport', label: '电解液传输性质计算' },
          { value: 'interface', label: '界面反应与SEI膜形成计算' },
        ];
      case 'fullBattery':
        return [
          { value: 'performance', label: '电池性能预测' },
          { value: 'interface', label: '电极/电解质界面模拟' },
          { value: 'lifetime', label: '循环寿命评估' },
          { value: 'safety', label: '安全性评估' },
        ];
      default:
        return [];
    }
  };
  
  // 处理计算类型选择变化
  const handleCalculationChange = (checkedValues: string[]) => {
    setSelectedCalculations(checkedValues);
  };
  
  // 下一步
  const nextStep = () => {
    setCurrentStep(currentStep + 1);
  };
  
  // 上一步
  const prevStep = () => {
    setCurrentStep(currentStep - 1);
  };
  
  // 提交计算
  const submitCalculation = () => {
    // 这里应该处理计算提交逻辑
    console.log('提交计算:', {
      moduleType,
      selectedCalculations,
      // 其他表单数据...
    });
    
    // 显示成功消息并跳转到最后一步
    setCurrentStep(3);
  };
  
  // 文件上传配置
  const uploadProps = {
    name: 'file',
    multiple: false,
    action: 'https://run.mocky.io/v3/435e224c-44fb-4773-9faf-380c5e6a2188',
    beforeUpload: (file: File) => {
      const isSupportedFormat = file.name.endsWith('.cif') || 
                               file.name.endsWith('.poscar') || 
                               file.name.endsWith('.pdb') || 
                               file.name.endsWith('.xyz') || 
                               file.name.endsWith('.cml');
      if (!isSupportedFormat) {
        message.error('只能上传CIF、POSCAR、PDB、XYZ或CML格式的文件!');
      }
      return isSupportedFormat || Upload.LIST_IGNORE;
    },
    onChange(info: any) {
      if (info.file.status === 'done') {
        setLoading(true);
        setTimeout(() => {
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
      return true;
    }
  };
  
  // 开始计算
  const startCalculation = () => {
    setShowCalculation(true);
  };
  
  // 渲染特定计算类型的组件
  const renderCalculationComponent = () => {
    // 如果还没有确认计算类型，显示选择界面
    if (!showCalculation) {
      return (
        <div>
          <Title level={4}>选择需要进行的计算</Title>
          <Paragraph>
            请选择一个或多个计算类型，系统将自动生成计算任务模板。
          </Paragraph>
          
          <Form layout="vertical">
            <Form.Item label="计算类型选择" required>
              <Checkbox.Group 
                options={getCalculationOptions()} 
                value={selectedCalculations}
                onChange={(values) => handleCalculationChange(values as string[])}
              />
            </Form.Item>
            
            <Form.Item>
              <Button 
                type="primary" 
                onClick={startCalculation}
                disabled={selectedCalculations.length === 0}
              >
                开始计算
              </Button>
            </Form.Item>
          </Form>
        </div>
      );
    }

    // 如果只选择了一种计算类型，直接渲染对应组件
    if (selectedCalculations.length === 1) {
      if (moduleType === 'electrode') {
        switch (selectedCalculations[0]) {
          case 'voltage':
            return <ElectrodeVoltageCalculation />;
          case 'formation':
            return <FormationEnergyCalculation />;
          case 'barrier':
            return <IonMigrationBarrierCalculation />;
          case 'electronic':
            return <ElectronicStructureCalculation />;
        }
      }
    }
    
    // 如果选择了多个计算类型，整合步骤
    if (selectedCalculations.length > 1) {
      return (
        <div>
          <Steps current={currentStep} style={{ marginBottom: 40 }}>
            <Step title="结构上传" icon={<UploadOutlined />} />
            <Step title="计算参数设置" icon={<SettingOutlined />} />
            <Step title="计算结果" icon={<LineChartOutlined />} />
            <Step title="完成" icon={<CheckCircleOutlined />} />
          </Steps>

          {currentStep === 0 && (
            <div>
              <Title level={4}>上传结构文件</Title>
              <Paragraph>
                请上传需要计算的材料晶体结构文件。此结构文件将用于所有选定的计算类型。
              </Paragraph>
              
              <Form form={form} layout="vertical">
                <Form.Item 
                  name="structureFile" 
                  rules={[{ required: true, message: '请上传结构文件' }]}
                >
                  <Upload.Dragger
                    {...uploadProps}
                    style={{ marginBottom: 24 }}
                  >
                    <p className="ant-upload-drag-icon">
                      <InboxOutlined />
                    </p>
                    <p className="ant-upload-text">点击或拖拽文件到此区域上传</p>
                    <p className="ant-upload-hint">
                      支持格式: .cif, .poscar, .pdb, .xyz, .cml
                    </p>
                  </Upload.Dragger>
                </Form.Item>
              </Form>

              {structurePreview && (
                <Card title="结构信息" style={{ marginTop: 16 }}>
                  {/* 结构预览内容保持不变 */}
                </Card>
              )}

              <Divider />
              
              <Button 
                type="primary" 
                onClick={nextStep} 
                disabled={!structurePreview}
              >
                下一步
              </Button>
            </div>
          )}

          {currentStep === 1 && (
            <div>
              <Title level={4}>计算参数设置</Title>
              <Paragraph>
                请为选定的计算类型设置相应的参数。
              </Paragraph>

              <Tabs defaultActiveKey="basic">
                <TabPane tab="基本参数" key="basic">
                  <Form layout="vertical">
                    <Row gutter={24}>
                      <Col span={12}>
                        <Form.Item label="计算名称" required>
                          <Input placeholder="请输入计算任务名称" />
                        </Form.Item>
                      </Col>
                      <Col span={12}>
                        <Form.Item label="计算描述">
                          <Input.TextArea placeholder="请简要描述计算目的（选填）" />
                        </Form.Item>
                      </Col>
                    </Row>

                    {selectedCalculations.includes('voltage') && (
                      <Card title="电极电位计算参数" style={{ marginBottom: 16 }}>
                        <Form.Item label="工作离子类型" required>
                          <Select defaultValue="Li" style={{ width: 200 }}>
                            <Option value="Li">Li+</Option>
                            <Option value="Na">Na+</Option>
                            <Option value="K">K+</Option>
                            <Option value="Mg">Mg2+</Option>
                          </Select>
                        </Form.Item>
                      </Card>
                    )}

                    {selectedCalculations.includes('barrier') && (
                      <Card title="离子迁移能垒计算参数" style={{ marginBottom: 16 }}>
                        <Form.Item label="NEB图像数量" required>
                          <InputNumber min={3} max={11} defaultValue={7} step={2} />
                        </Form.Item>
                      </Card>
                    )}

                    {selectedCalculations.includes('electronic') && (
                      <Card title="电子结构计算参数" style={{ marginBottom: 16 }}>
                        <Form.Item label="交换关联泛函" required>
                          <Select defaultValue="pbe">
                            <Option value="pbe">PBE</Option>
                            <Option value="hse">HSE06</Option>
                            <Option value="b3lyp">B3LYP</Option>
                          </Select>
                        </Form.Item>
                      </Card>
                    )}

                    <Form.Item label="计算精度" required>
                      <Radio.Group defaultValue="standard">
                        <Radio value="fast">快速（精度低）</Radio>
                        <Radio value="standard">标准</Radio>
                        <Radio value="accurate">精确（耗时长）</Radio>
                      </Radio.Group>
                    </Form.Item>
                  </Form>
                </TabPane>

                <TabPane tab="高级参数" key="advanced">
                  <Alert
                    message="高级参数说明"
                    description="这些参数适用于有经验的用户，如果您不确定，请保留默认值。"
                    type="info"
                    showIcon
                    style={{ marginBottom: 16 }}
                  />

                  <Form layout="vertical">
                    <Row gutter={24}>
                      <Col span={12}>
                        <Form.Item label="计算引擎">
                          <Select defaultValue="vasp">
                            <Option value="vasp">VASP</Option>
                            <Option value="qe">Quantum ESPRESSO</Option>
                            <Option value="gaussian">Gaussian</Option>
                          </Select>
                        </Form.Item>
                      </Col>
                      <Col span={12}>
                        <Form.Item label="并行核心数">
                          <Select defaultValue="16">
                            <Option value="8">8</Option>
                            <Option value="16">16</Option>
                            <Option value="32">32</Option>
                            <Option value="64">64</Option>
                          </Select>
                        </Form.Item>
                      </Col>
                    </Row>
                  </Form>
                </TabPane>
              </Tabs>

              <Divider />
              
              <Space>
                <Button onClick={prevStep}>上一步</Button>
                <Button type="primary" onClick={nextStep}>开始计算</Button>
              </Space>
            </div>
          )}

          {currentStep === 2 && (
            <div>
              <Title level={4}>计算结果</Title>
              
              {loading ? (
                <div style={{ textAlign: 'center', margin: '40px 0' }}>
                  <Spin size="large" tip="正在进行计算..." />
                </div>
              ) : (
                <div>
                  <Alert
                    message="计算完成"
                    description="所有选定的计算任务已完成，结果如下所示。"
                    type="success"
                    showIcon
                    style={{ marginBottom: 24 }}
                  />

                  <Tabs>
                    {selectedCalculations.includes('voltage') && (
                      <TabPane tab="电极电位计算结果" key="voltage">
                        <Card title="电极电位计算结果">
                          {/* 电极电位计算结果内容 */}
                        </Card>
                      </TabPane>
                    )}

                    {selectedCalculations.includes('barrier') && (
                      <TabPane tab="离子迁移能垒计算结果" key="barrier">
                        <Card title="离子迁移能垒计算结果">
                          {/* 离子迁移能垒计算结果内容 */}
                        </Card>
                      </TabPane>
                    )}

                    {selectedCalculations.includes('electronic') && (
                      <TabPane tab="电子结构计算结果" key="electronic">
                        <Card title="电子结构计算结果">
                          {/* 电子结构计算结果内容 */}
                        </Card>
                      </TabPane>
                    )}
                  </Tabs>

                  <Divider />
                  
                  <Space>
                    <Button onClick={prevStep}>上一步</Button>
                    <Button type="primary" onClick={submitCalculation}>保存结果</Button>
                  </Space>
                </div>
              )}
            </div>
          )}

          {currentStep === 3 && (
            <div style={{ textAlign: 'center', padding: '20px 0' }}>
              <CheckCircleOutlined style={{ fontSize: 72, color: '#52c41a', marginBottom: 24 }} />
              <Title level={3}>计算结果已保存！</Title>
              <Paragraph>
                您选择的所有计算任务已完成并保存。您可以在用户中心的"我的计算任务"中查看详细结果。
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
          )}
        </div>
      );
    }
  };
  
  return (
    <div>
      <div style={{ marginBottom: 24 }}>
        <Title level={2}>{getTitle()}</Title>
        <Paragraph style={{ fontSize: 16 }}>
          使用此工作台配置和执行您的计算任务。选择所需的计算类型，设置参数，提交计算，并查看结果。
        </Paragraph>
      </div>
      
      <Card style={{ marginBottom: 24 }}>
        {renderCalculationComponent()}
      </Card>
    </div>
  );
};

export default Workbench; 