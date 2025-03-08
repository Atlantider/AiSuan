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
  Checkbox
} from 'antd';
import { 
  UploadOutlined, 
  InboxOutlined, 
  ThunderboltOutlined, 
  SettingOutlined, 
  LineChartOutlined, 
  CheckCircleOutlined 
} from '@ant-design/icons';
import type { UploadProps } from 'antd';

const { Title, Paragraph, Text } = Typography;
const { Step } = Steps;
const { Option } = Select;
const { Dragger } = Upload;

const ElectrodeVoltageCalculation: React.FC = () => {
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
    return (
      <div>
        <Title level={4}>上传结构文件</Title>
        <Paragraph>
          请上传电极材料的晶体结构文件。支持的格式包括CIF、POSCAR、PDB、XYZ和CML。
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
  
  // 步骤2: 离子配置
  const renderIonConfiguration = () => {
    return (
      <div>
        <Title level={4}>离子与浓度设置</Title>
        <Paragraph>
          选择工作离子类型并设置浓度范围。
        </Paragraph>
        
        <Form form={form} layout="vertical">
          <Form.Item 
            name="ionType" 
            label="工作离子类型" 
            initialValue="Li" 
            rules={[{ required: true, message: '请选择工作离子类型' }]}
          >
            <Radio.Group>
              <Radio value="Li">Li⁺</Radio>
              <Radio value="Na">Na⁺</Radio>
              <Radio value="K">K⁺</Radio>
              <Radio value="Mg">Mg²⁺</Radio>
              <Radio value="Ca">Ca²⁺</Radio>
              <Radio value="Al">Al³⁺</Radio>
            </Radio.Group>
          </Form.Item>
          
          <Form.Item 
            label="离子插入比例范围" 
            required
          >
            <Row gutter={16}>
              <Col span={16}>
                <Slider
                  range
                  min={0}
                  max={1}
                  step={0.01}
                  value={concentrationRange}
                  onChange={(value) => setConcentrationRange(value as [number, number])}
                />
              </Col>
              <Col span={8}>
                <InputNumber
                  min={0}
                  max={1}
                  step={0.01}
                  value={concentrationRange[0]}
                  onChange={(value) => setConcentrationRange([value || 0, concentrationRange[1]])}
                  style={{ marginRight: 8 }}
                />
                <span style={{ margin: '0 8px' }}>-</span>
                <InputNumber
                  min={0}
                  max={1}
                  step={0.01}
                  value={concentrationRange[1]}
                  onChange={(value) => setConcentrationRange([concentrationRange[0], value || 1])}
                />
              </Col>
            </Row>
          </Form.Item>
          
          <Form.Item 
            label="浓度步数" 
            required
          >
            <Row gutter={16}>
              <Col span={16}>
                <Slider
                  min={2}
                  max={10}
                  value={concentrationSteps}
                  onChange={(value) => setConcentrationSteps(value)}
                />
              </Col>
              <Col span={8}>
                <InputNumber
                  min={2}
                  max={10}
                  value={concentrationSteps}
                  onChange={(value) => setConcentrationSteps(value || 5)}
                />
              </Col>
            </Row>
          </Form.Item>
        </Form>
        
        <Alert
          message="浓度设置说明"
          description={`系统将在${concentrationRange[0]}到${concentrationRange[1]}的范围内生成${concentrationSteps}个不同浓度的模型进行计算。`}
          type="info"
          showIcon
          style={{ marginBottom: 16 }}
        />
        
        <Divider />
        
        <Space>
          <Button onClick={prevStep}>上一步</Button>
          <Button type="primary" onClick={nextStep}>下一步</Button>
        </Space>
      </div>
    );
  };
  
  // 步骤3: 超胞设置
  const renderSupercellSetup = () => {
    return (
      <div>
        <Title level={4}>超胞设置</Title>
        <Paragraph>
          设置超胞大小，以确保计算模型包含足够的原子。
        </Paragraph>
        
        <Form form={form} layout="vertical">
          <Form.Item 
            name="supercellMethod" 
            label="超胞设置方式" 
            initialValue="auto"
          >
            <Radio.Group onChange={(e) => setSupercellMethod(e.target.value)}>
              <Radio value="auto">自动规划</Radio>
              <Radio value="manual">手动设置</Radio>
            </Radio.Group>
          </Form.Item>
          
          {supercellMethod === 'auto' ? (
            <Form.Item 
              label="目标原子数" 
              required
            >
              <Row gutter={16}>
                <Col span={16}>
                  <Slider
                    min={50}
                    max={300}
                    step={10}
                    value={targetAtomCount}
                    onChange={(value) => setTargetAtomCount(value)}
                  />
                </Col>
                <Col span={8}>
                  <InputNumber
                    min={50}
                    max={300}
                    step={10}
                    value={targetAtomCount}
                    onChange={(value) => setTargetAtomCount(value || 100)}
                  />
                </Col>
              </Row>
              <Text type="secondary">推荐值: 100-200个原子</Text>
            </Form.Item>
          ) : (
            <Form.Item 
              label="超胞尺寸" 
              required
            >
              <Row gutter={16}>
                <Col span={8}>
                  <Form.Item label="a方向" noStyle>
                    <InputNumber
                      min={1}
                      max={10}
                      value={supercellDimensions.a}
                      onChange={(value) => setSupercellDimensions({...supercellDimensions, a: value || 1})}
                      addonAfter="倍"
                    />
                  </Form.Item>
                </Col>
                <Col span={8}>
                  <Form.Item label="b方向" noStyle>
                    <InputNumber
                      min={1}
                      max={10}
                      value={supercellDimensions.b}
                      onChange={(value) => setSupercellDimensions({...supercellDimensions, b: value || 1})}
                      addonAfter="倍"
                    />
                  </Form.Item>
                </Col>
                <Col span={8}>
                  <Form.Item label="c方向" noStyle>
                    <InputNumber
                      min={1}
                      max={10}
                      value={supercellDimensions.c}
                      onChange={(value) => setSupercellDimensions({...supercellDimensions, c: value || 1})}
                      addonAfter="倍"
                    />
                  </Form.Item>
                </Col>
              </Row>
            </Form.Item>
          )}
        </Form>
        
        <Card title="超胞预览" style={{ marginTop: 16 }}>
          <Row gutter={24}>
            <Col span={12}>
              <p><Text strong>超胞尺寸:</Text> {supercellMethod === 'auto' ? '3 × 3 × 1' : `${supercellDimensions.a} × ${supercellDimensions.b} × ${supercellDimensions.c}`}</p>
              <p><Text strong>原子数:</Text> {supercellMethod === 'auto' ? 72 : 8 * supercellDimensions.a * supercellDimensions.b * supercellDimensions.c}</p>
              <p><Text strong>体积:</Text> {supercellMethod === 'auto' ? '868.32 Å³' : `${(96.48 * supercellDimensions.a * supercellDimensions.b * supercellDimensions.c).toFixed(2)} Å³`}</p>
            </Col>
            <Col span={12}>
              <div style={{ height: 150, background: '#f0f2f5', display: 'flex', justifyContent: 'center', alignItems: 'center' }}>
                <Text type="secondary">超胞结构预览区域</Text>
              </div>
            </Col>
          </Row>
        </Card>
        
        <Divider />
        
        <Space>
          <Button onClick={prevStep}>上一步</Button>
          <Button type="primary" onClick={nextStep}>下一步</Button>
        </Space>
      </div>
    );
  };
  
  // 步骤4: 模型预览
  const renderModelPreview = () => {
    // 生成模拟的浓度模型数据
    const concentrationStep = (concentrationRange[1] - concentrationRange[0]) / (concentrationSteps - 1);
    const models = Array.from({ length: concentrationSteps }, (_, i) => {
      const concentration = concentrationRange[0] + i * concentrationStep;
      return {
        concentration,
        formula: `Li${concentration.toFixed(2)}CoO₂`,
        atomCount: supercellMethod === 'auto' ? 72 : 8 * supercellDimensions.a * supercellDimensions.b * supercellDimensions.c,
      };
    });
    
    return (
      <div>
        <Title level={4}>模型生成与预览</Title>
        <Paragraph>
          系统已根据您的设置生成以下浓度模型。请检查并确认这些模型。
        </Paragraph>
        
        {loading && (
          <div style={{ textAlign: 'center', margin: '20px 0' }}>
            <Spin tip="正在生成模型..." />
          </div>
        )}
        
        {modelsGenerated && (
          <>
            <div style={{ marginBottom: 16 }}>
              <Radio.Group defaultValue="0">
                {models.map((model, index) => (
                  <Radio.Button value={index.toString()} key={index}>
                    x = {model.concentration.toFixed(2)}
                  </Radio.Button>
                ))}
              </Radio.Group>
            </div>
            
            <Card title="当前模型预览" style={{ marginBottom: 16 }}>
              <Row gutter={24}>
                <Col span={12}>
                  <p><Text strong>化学式:</Text> {models[0].formula}</p>
                  <p><Text strong>浓度/比例:</Text> x = {models[0].concentration.toFixed(2)}</p>
                  <p><Text strong>原子数:</Text> {models[0].atomCount}</p>
                  <p><Text strong>体积:</Text> 868.32 Å³</p>
                  <p><Text strong>理论密度:</Text> 5.06 g/cm³</p>
                </Col>
                <Col span={12}>
                  <div style={{ height: 200, background: '#f0f2f5', display: 'flex', justifyContent: 'center', alignItems: 'center' }}>
                    <Text type="secondary">模型结构预览区域</Text>
                  </div>
                </Col>
              </Row>
            </Card>
            
            <Alert
              message="模型生成完成"
              description={`已成功生成${models.length}个不同浓度的模型，点击上方按钮可切换查看不同浓度的模型。`}
              type="success"
              showIcon
              style={{ marginBottom: 16 }}
            />
          </>
        )}
        
        <Divider />
        
        <Space>
          <Button onClick={prevStep}>上一步</Button>
          <Button type="primary" onClick={nextStep}>下一步</Button>
        </Space>
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
    <div>
      <Card style={{ marginBottom: 24 }}>
        <Steps current={currentStep} style={{ marginBottom: 40 }}>
          <Step title="结构上传" icon={<UploadOutlined />} />
          <Step title="离子设置" icon={<ThunderboltOutlined />} />
          <Step title="超胞配置" icon={<SettingOutlined />} />
          <Step title="模型预览" icon={<LineChartOutlined />} />
          <Step title="计算参数" icon={<SettingOutlined />} />
          <Step title="完成" icon={<CheckCircleOutlined />} />
        </Steps>
        
        {renderStepContent()}
      </Card>
    </div>
  );
};

export default ElectrodeVoltageCalculation; 