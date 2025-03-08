import React, { useState } from 'react';
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
  Tabs
} from 'antd';
import { 
  UploadOutlined, 
  InboxOutlined, 
  ThunderboltOutlined, 
  SettingOutlined, 
  LineChartOutlined, 
  CheckCircleOutlined,
  AppstoreOutlined,
  BarChartOutlined
} from '@ant-design/icons';
import type { UploadProps } from 'antd';

const { Title, Paragraph, Text } = Typography;
const { Step } = Steps;
const { Option } = Select;
const { Dragger } = Upload;
const { TabPane } = Tabs;

const ElectronicStructureCalculation: React.FC = () => {
  const [currentStep, setCurrentStep] = useState(0);
  const [form] = Form.useForm();
  const [fileList, setFileList] = useState<any[]>([]);
  const [structurePreview, setStructurePreview] = useState<boolean>(false);
  const [loading, setLoading] = useState<boolean>(false);
  const [calculationType, setCalculationType] = useState<string[]>(['bandstructure', 'dos']);
  const [bandgap, setBandgap] = useState<number | null>(null);
  
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
      return isCIF || Upload.LIST_IGNORE;
    },
    onChange(info) {
      setFileList(info.fileList);
      
      // 当文件上传成功时，模拟结构解析
      if (info.file.status === 'done') {
        setLoading(true);
        // 模拟结构解析过程
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
        if (currentStep === 1) {
          // 模拟计算过程
          setLoading(true);
          setTimeout(() => {
            // 模拟计算结果
            setBandgap(2.1);
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
      setCurrentStep(3);
    }, 2000);
  };
  
  // 渲染步骤内容
  const renderStepContent = () => {
    switch (currentStep) {
      case 0:
        return renderStructureUpload();
      case 1:
        return renderCalculationSetup();
      case 2:
        return renderResults();
      case 3:
        return renderSubmission();
      default:
        return null;
    }
  };
  
  // 步骤1: 结构上传
  const renderStructureUpload = () => {
    return (
      <div>
        <Title level={4}>上传结构文件</Title>
        <Paragraph>
          请上传需要计算电子结构的材料晶体结构文件。支持的格式包括CIF、POSCAR、PDB、XYZ和CML。
        </Paragraph>
        
        <Form form={form} layout="vertical">
          <Form.Item 
            name="structureFile" 
            rules={[{ required: true, message: '请上传结构文件' }]}
          >
            <Dragger {...uploadProps}>
              <p className="ant-upload-drag-icon">
                <InboxOutlined />
              </p>
              <p className="ant-upload-text">点击或拖拽文件到此区域上传</p>
              <p className="ant-upload-hint">
                支持格式: .cif, .poscar, .pdb, .xyz, .cml
              </p>
            </Dragger>
          </Form.Item>
        </Form>
        
        {loading && (
          <div style={{ textAlign: 'center', margin: '20px 0' }}>
            <Spin tip="正在解析结构文件..." />
          </div>
        )}
        
        {structurePreview && (
          <Card title="结构信息" style={{ marginTop: 16 }}>
            <Row gutter={24}>
              <Col span={12}>
                <p><Text strong>化学式:</Text> LiCoO₂</p>
                <p><Text strong>空间群:</Text> R-3m (166)</p>
                <p><Text strong>晶胞参数:</Text></p>
                <ul>
                  <li>a = 2.8156 Å</li>
                  <li>b = 2.8156 Å</li>
                  <li>c = 14.0542 Å</li>
                  <li>α = 90°</li>
                  <li>β = 90°</li>
                  <li>γ = 120°</li>
                </ul>
              </Col>
              <Col span={12}>
                <p><Text strong>体积:</Text> 96.48 Å³</p>
                <p><Text strong>原子数:</Text> 8</p>
                <p><Text strong>密度:</Text> 5.06 g/cm³</p>
                <p><Text strong>元素组成:</Text> Li (25%), Co (25%), O (50%)</p>
                <div style={{ height: 150, background: '#f0f2f5', display: 'flex', justifyContent: 'center', alignItems: 'center' }}>
                  <Text type="secondary">3D结构预览区域</Text>
                </div>
              </Col>
            </Row>
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
    );
  };
  
  // 步骤2: 计算设置
  const renderCalculationSetup = () => {
    return (
      <div>
        <Title level={4}>电子结构计算设置</Title>
        <Paragraph>
          选择需要计算的电子结构属性和计算参数。
        </Paragraph>
        
        <Form form={form} layout="vertical">
          <Form.Item 
            name="calculationType" 
            label="计算类型" 
            initialValue={['bandstructure', 'dos']}
            rules={[{ required: true, message: '请选择至少一种计算类型' }]}
          >
            <Checkbox.Group 
              onChange={(values) => setCalculationType(values as string[])}
              style={{ display: 'flex', flexDirection: 'column', gap: '8px' }}
            >
              <Checkbox value="bandstructure">能带结构</Checkbox>
              <Checkbox value="dos">态密度</Checkbox>
              <Checkbox value="charge">电荷密度</Checkbox>
              <Checkbox value="optical">光学性质</Checkbox>
              <Checkbox value="effective_mass">有效质量</Checkbox>
            </Checkbox.Group>
          </Form.Item>
          
          <Form.Item 
            name="functional" 
            label="交换关联泛函" 
            initialValue="pbe"
          >
            <Select style={{ width: 200 }}>
              <Option value="pbe">PBE</Option>
              <Option value="pbesol">PBEsol</Option>
              <Option value="b3lyp">B3LYP</Option>
              <Option value="hse06">HSE06</Option>
              <Option value="gw">GW近似</Option>
            </Select>
          </Form.Item>
          
          <Form.Item 
            name="spinPolarized" 
            label="自旋极化计算" 
            initialValue={true}
          >
            <Radio.Group>
              <Radio value={true}>是</Radio>
              <Radio value={false}>否</Radio>
            </Radio.Group>
          </Form.Item>
          
          <Form.Item 
            name="kpointsDensity" 
            label="K点网格密度" 
            initialValue={40}
          >
            <Slider
              min={10}
              max={100}
              marks={{
                10: '粗略',
                40: '标准',
                70: '精细',
                100: '超精细'
              }}
            />
          </Form.Item>
          
          <Form.Item 
            name="encut" 
            label="平面波截断能 (ENCUT)" 
            initialValue={520}
          >
            <InputNumber min={300} max={800} addonAfter="eV" style={{ width: 200 }} />
          </Form.Item>
          
          <Form.Item 
            name="hubbardU" 
            label="Hubbard U校正"
            initialValue={['Co']}
          >
            <Checkbox.Group>
              <Checkbox value="Co">Co (U = 3.32 eV)</Checkbox>
              <Checkbox value="Fe">Fe (U = 4.0 eV)</Checkbox>
              <Checkbox value="Mn">Mn (U = 3.9 eV)</Checkbox>
            </Checkbox.Group>
          </Form.Item>
          
          <Divider orientation="left">计算资源</Divider>
          
          <Row gutter={24}>
            <Col span={12}>
              <Form.Item 
                name="cpuCores" 
                label="CPU核心数" 
                initialValue={16}
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
                initialValue={32}
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
        </Form>
        
        <Divider />
        
        <Space>
          <Button onClick={prevStep}>上一步</Button>
          <Button type="primary" onClick={nextStep}>开始计算</Button>
        </Space>
      </div>
    );
  };
  
  // 步骤3: 计算结果
  const renderResults = () => {
    return (
      <div>
        <Title level={4}>电子结构计算结果</Title>
        
        {loading ? (
          <div style={{ textAlign: 'center', margin: '40px 0' }}>
            <Spin size="large" tip="正在计算电子结构..." />
          </div>
        ) : (
          <>
            <Alert
              message="计算完成"
              description="电子结构计算已成功完成，结果如下所示。"
              type="success"
              showIcon
              style={{ marginBottom: 24 }}
            />
            
            <Card title="计算结果概要" style={{ marginBottom: 24 }}>
              <Row gutter={[24, 24]}>
                <Col span={12}>
                  <Card type="inner" title="带隙">
                    <div style={{ textAlign: 'center', padding: '20px 0' }}>
                      <Title level={2} style={{ color: '#1C64F2' }}>{bandgap} eV</Title>
                      <Text type="secondary">半导体材料</Text>
                    </div>
                  </Card>
                </Col>
                <Col span={12}>
                  <Card type="inner" title="电子性质评估">
                    <div style={{ padding: '10px 0' }}>
                      <p><Text strong>导电性:</Text> 半导体</p>
                      <p><Text strong>载流子类型:</Text> p型</p>
                      <p><Text strong>有效质量:</Text> 电子: 0.8 m₀, 空穴: 1.2 m₀</p>
                      <p><Text strong>光学吸收边:</Text> 590 nm</p>
                    </div>
                  </Card>
                </Col>
              </Row>
            </Card>
            
            <Tabs defaultActiveKey="bandstructure">
              {calculationType.includes('bandstructure') && (
                <TabPane tab="能带结构" key="bandstructure">
                  <div style={{ height: 400, background: '#f0f2f5', display: 'flex', justifyContent: 'center', alignItems: 'center', marginBottom: 16 }}>
                    <Text type="secondary">能带结构图表区域</Text>
                  </div>
                  <p><Text type="secondary">能带结构显示了电子在晶体中允许的能量状态。费米能级被设置为0 eV。</Text></p>
                </TabPane>
              )}
              
              {calculationType.includes('dos') && (
                <TabPane tab="态密度" key="dos">
                  <div style={{ height: 400, background: '#f0f2f5', display: 'flex', justifyContent: 'center', alignItems: 'center', marginBottom: 16 }}>
                    <Text type="secondary">态密度图表区域</Text>
                  </div>
                  <p><Text type="secondary">态密度显示了每个能量水平上可用的电子态数量。</Text></p>
                </TabPane>
              )}
              
              {calculationType.includes('charge') && (
                <TabPane tab="电荷密度" key="charge">
                  <div style={{ height: 400, background: '#f0f2f5', display: 'flex', justifyContent: 'center', alignItems: 'center', marginBottom: 16 }}>
                    <Text type="secondary">电荷密度图表区域</Text>
                  </div>
                  <p><Text type="secondary">电荷密度显示了电子在晶体中的空间分布。</Text></p>
                </TabPane>
              )}
              
              {calculationType.includes('optical') && (
                <TabPane tab="光学性质" key="optical">
                  <div style={{ height: 400, background: '#f0f2f5', display: 'flex', justifyContent: 'center', alignItems: 'center', marginBottom: 16 }}>
                    <Text type="secondary">光学性质图表区域</Text>
                  </div>
                  <p><Text type="secondary">光学性质包括吸收系数、介电函数和反射率等。</Text></p>
                </TabPane>
              )}
            </Tabs>
            
            <Divider />
            
            <Card title="计算详情" style={{ marginBottom: 24 }}>
              <p><Text strong>交换关联泛函:</Text> PBE</p>
              <p><Text strong>K点网格:</Text> 8×8×8</p>
              <p><Text strong>平面波截断能:</Text> 520 eV</p>
              <p><Text strong>自旋极化:</Text> 是</p>
              <p><Text strong>Hubbard U校正:</Text> Co (U = 3.32 eV)</p>
              <p><Text strong>收敛状态:</Text> 已收敛 (力收敛标准: 0.01 eV/Å)</p>
            </Card>
            
            <Divider />
            
            <Space>
              <Button onClick={prevStep}>上一步</Button>
              <Button type="primary" onClick={submitCalculation}>保存结果</Button>
            </Space>
          </>
        )}
      </div>
    );
  };
  
  // 步骤4: 提交确认
  const renderSubmission = () => {
    return (
      <div style={{ textAlign: 'center', padding: '20px 0' }}>
        <CheckCircleOutlined style={{ fontSize: 72, color: '#52c41a', marginBottom: 24 }} />
        <Title level={3}>计算结果已保存！</Title>
        <Paragraph>
          您的电子结构计算结果已成功保存。您可以在用户中心的"我的计算任务"中查看详细结果。
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
    <div>
      <Card style={{ marginBottom: 24 }}>
        <Steps current={currentStep} style={{ marginBottom: 40 }}>
          <Step title="结构上传" icon={<UploadOutlined />} />
          <Step title="计算设置" icon={<SettingOutlined />} />
          <Step title="计算结果" icon={<BarChartOutlined />} />
          <Step title="完成" icon={<CheckCircleOutlined />} />
        </Steps>
        
        {renderStepContent()}
      </Card>
    </div>
  );
};

export default ElectronicStructureCalculation; 