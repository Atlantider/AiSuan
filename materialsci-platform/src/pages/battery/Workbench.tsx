import React, { useState } from 'react';
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
  Alert
} from 'antd';
import { 
  UploadOutlined, 
  ExperimentOutlined, 
  SettingOutlined, 
  FileTextOutlined,
  CheckCircleOutlined,
  InfoCircleOutlined
} from '@ant-design/icons';
import { useParams, useNavigate } from 'react-router-dom';

const { Title, Paragraph, Text } = Typography;
const { TabPane } = Tabs;
const { Option } = Select;
const { Step } = Steps;

interface WorkbenchProps {
  calculationType?: string;
}

const Workbench: React.FC<WorkbenchProps> = ({ calculationType }) => {
  const navigate = useNavigate();
  const [currentStep, setCurrentStep] = useState(0);
  const [selectedCalculations, setSelectedCalculations] = useState<string[]>([]);
  
  // 使用传入的calculationType或从URL参数中获取
  const params = useParams();
  const moduleType = calculationType || params.type || 'electrode';
  
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
  
  return (
    <div>
      <div style={{ marginBottom: 24 }}>
        <Title level={2}>{getTitle()}</Title>
        <Paragraph style={{ fontSize: 16 }}>
          使用此工作台配置和执行您的计算任务。选择所需的计算类型，设置参数，提交计算，并查看结果。
        </Paragraph>
      </div>
      
      <Card style={{ marginBottom: 24 }}>
        <Steps current={currentStep} style={{ marginBottom: 40 }}>
          <Step title="选择计算内容" icon={<ExperimentOutlined />} />
          <Step title="配置参数" icon={<SettingOutlined />} />
          <Step title="提交计算" icon={<FileTextOutlined />} />
          <Step title="完成" icon={<CheckCircleOutlined />} />
        </Steps>
        
        {currentStep === 0 && (
          <div>
            <Title level={4}>选择需要进行的计算</Title>
            <Paragraph>
              您可以选择一个或多个计算类型，系统将自动生成计算任务模板。
            </Paragraph>
            
            <Form layout="vertical">
              <Form.Item label="计算类型选择" required>
                <Checkbox.Group 
                  options={getCalculationOptions()} 
                  onChange={(values) => handleCalculationChange(values as string[])}
                />
              </Form.Item>
              
              <Divider />
              
              <Form.Item>
                <Button 
                  type="primary" 
                  onClick={nextStep} 
                  disabled={selectedCalculations.length === 0}
                >
                  下一步
                </Button>
              </Form.Item>
            </Form>
          </div>
        )}
        
        {currentStep === 1 && (
          <div>
            <Title level={4}>配置计算参数</Title>
            <Paragraph>
              请上传材料结构文件并设置必要的计算参数。
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
                  
                  <Form.Item label="结构文件上传" required>
                    <Upload>
                      <Button icon={<UploadOutlined />}>上传结构文件 (CIF, POSCAR等)</Button>
                    </Upload>
                  </Form.Item>
                  
                  {moduleType === 'electrode' && (
                    <>
                      <Form.Item label="工作离子类型" required>
                        <Select defaultValue="Li" style={{ width: 200 }}>
                          <Option value="Li">Li+</Option>
                          <Option value="Na">Na+</Option>
                          <Option value="K">K+</Option>
                          <Option value="Mg">Mg2+</Option>
                        </Select>
                      </Form.Item>
                      
                      <Form.Item label="计算精度" required>
                        <Radio.Group defaultValue="standard">
                          <Radio value="fast">快速（精度低）</Radio>
                          <Radio value="standard">标准</Radio>
                          <Radio value="accurate">精确（耗时长）</Radio>
                        </Radio.Group>
                      </Form.Item>
                    </>
                  )}
                  
                  {moduleType === 'electrolyte' && (
                    <Form.Item label="溶剂类型" required>
                      <Select mode="multiple" placeholder="选择溶剂类型" style={{ width: '100%' }}>
                        <Option value="EC">碳酸乙烯酯 (EC)</Option>
                        <Option value="DMC">碳酸二甲酯 (DMC)</Option>
                        <Option value="EMC">碳酸甲乙酯 (EMC)</Option>
                        <Option value="DEC">碳酸二乙酯 (DEC)</Option>
                      </Select>
                    </Form.Item>
                  )}
                </Form>
              </TabPane>
              
              <TabPane tab="高级参数" key="advanced">
                <Alert
                  message="高级参数说明"
                  description="这些参数适用于有经验的用户，如果您不确定，请保留默认值。"
                  type="info"
                  showIcon
                  icon={<InfoCircleOutlined />}
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
                        <Select defaultValue="8">
                          <Option value="4">4</Option>
                          <Option value="8">8</Option>
                          <Option value="16">16</Option>
                          <Option value="32">32</Option>
                        </Select>
                      </Form.Item>
                    </Col>
                  </Row>
                  
                  <Form.Item label="其他参数">
                    <Input.TextArea placeholder="请输入其他高级参数（可选）" rows={4} />
                  </Form.Item>
                </Form>
              </TabPane>
            </Tabs>
            
            <Divider />
            
            <Space>
              <Button onClick={prevStep}>上一步</Button>
              <Button type="primary" onClick={nextStep}>下一步</Button>
            </Space>
          </div>
        )}
        
        {currentStep === 2 && (
          <div>
            <Title level={4}>提交计算</Title>
            <Paragraph>
              请仔细检查以下信息，确认无误后提交计算。
            </Paragraph>
            
            <Card title="计算任务概要" style={{ marginBottom: 16 }}>
              <p><Text strong>计算模块：</Text> {moduleType === 'electrode' ? '电极材料计算' : moduleType === 'electrolyte' ? '电解液计算' : '全电池系统计算'}</p>
              <p><Text strong>选择的计算类型：</Text></p>
              <ul>
                {selectedCalculations.map(calc => (
                  <li key={calc}>{getCalculationOptions().find(option => option.value === calc)?.label}</li>
                ))}
              </ul>
              <p><Text strong>计算引擎：</Text> VASP</p>
              <p><Text strong>并行核心数：</Text> 8</p>
              <p><Text strong>预计计算时间：</Text> 2-4小时</p>
            </Card>
            
            <Alert
              message="计算资源提示"
              description="此计算将消耗您账户中的计算资源。当前账户剩余计算时数：40小时。"
              type="info"
              showIcon
              style={{ marginBottom: 16 }}
            />
            
            <Divider />
            
            <Space>
              <Button onClick={prevStep}>上一步</Button>
              <Button type="primary" onClick={submitCalculation}>提交计算</Button>
            </Space>
          </div>
        )}
        
        {currentStep === 3 && (
          <div style={{ textAlign: 'center', padding: '20px 0' }}>
            <CheckCircleOutlined style={{ fontSize: 72, color: '#52c41a', marginBottom: 24 }} />
            <Title level={3}>计算任务已成功提交！</Title>
            <Paragraph>
              您的计算任务已加入计算队列，系统将自动处理您的计算请求。您可以在用户中心的"我的计算任务"中查看任务进度。
            </Paragraph>
            <div style={{ marginTop: 24 }}>
              <Space size="large">
                <Button type="primary" onClick={() => navigate('/user/tasks')}>
                  查看我的任务
                </Button>
                <Button onClick={() => navigate('/')}>
                  返回首页
                </Button>
              </Space>
            </div>
          </div>
        )}
      </Card>
    </div>
  );
};

export default Workbench; 