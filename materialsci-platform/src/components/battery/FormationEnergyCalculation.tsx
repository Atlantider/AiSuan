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
  Table
} from 'antd';
import { 
  UploadOutlined, 
  InboxOutlined, 
  ThunderboltOutlined, 
  SettingOutlined, 
  LineChartOutlined, 
  CheckCircleOutlined,
  CalculatorOutlined,
  DatabaseOutlined
} from '@ant-design/icons';
import type { UploadProps } from 'antd';
import type { ColumnsType } from 'antd/es/table';

const { Title, Paragraph, Text } = Typography;
const { Step } = Steps;
const { Option } = Select;
const { Dragger } = Upload;

interface ElementReference {
  element: string;
  structure: string;
  energy: number;
}

const FormationEnergyCalculation: React.FC = () => {
  const [currentStep, setCurrentStep] = useState(0);
  const [form] = Form.useForm();
  const [fileList, setFileList] = useState<any[]>([]);
  const [structurePreview, setStructurePreview] = useState<boolean>(false);
  const [loading, setLoading] = useState<boolean>(false);
  const [referenceElements, setReferenceElements] = useState<ElementReference[]>([
    { element: 'Li', structure: '体心立方', energy: -1.90 },
    { element: 'Co', structure: '六方密堆', energy: -7.10 },
    { element: 'O', structure: '气态分子', energy: -4.94 },
  ]);
  const [calculationMethod, setCalculationMethod] = useState<string>('dft');
  const [formationEnergy, setFormationEnergy] = useState<number | null>(null);
  
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
          // 模拟计算过程
          setLoading(true);
          setTimeout(() => {
            // 模拟计算结果
            setFormationEnergy(-2.34);
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
        return renderReferenceElements();
      case 2:
        return renderCalculationMethod();
      case 3:
        return renderResults();
      case 4:
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
          请上传需要计算形成能的材料晶体结构文件。支持的格式包括CIF、POSCAR、PDB、XYZ和CML。
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
  
  // 步骤2: 参考元素设置
  const renderReferenceElements = () => {
    const columns: ColumnsType<ElementReference> = [
      {
        title: '元素',
        dataIndex: 'element',
        key: 'element',
      },
      {
        title: '参考结构',
        dataIndex: 'structure',
        key: 'structure',
        render: (text, record) => (
          <Select 
            defaultValue={record.structure} 
            style={{ width: 150 }}
            onChange={(value) => {
              const newReferences = [...referenceElements];
              const index = newReferences.findIndex(item => item.element === record.element);
              if (index !== -1) {
                newReferences[index].structure = value;
                setReferenceElements(newReferences);
              }
            }}
          >
            <Option value="体心立方">体心立方</Option>
            <Option value="面心立方">面心立方</Option>
            <Option value="六方密堆">六方密堆</Option>
            <Option value="气态分子">气态分子</Option>
            <Option value="金刚石结构">金刚石结构</Option>
          </Select>
        ),
      },
      {
        title: '参考能量 (eV/atom)',
        dataIndex: 'energy',
        key: 'energy',
        render: (text, record) => (
          <InputNumber
            defaultValue={record.energy}
            step={0.01}
            style={{ width: 120 }}
            onChange={(value) => {
              const newReferences = [...referenceElements];
              const index = newReferences.findIndex(item => item.element === record.element);
              if (index !== -1 && value !== null) {
                newReferences[index].energy = value;
                setReferenceElements(newReferences);
              }
            }}
          />
        ),
      },
    ];
    
    return (
      <div>
        <Title level={4}>参考元素设置</Title>
        <Paragraph>
          形成能计算需要设置组成元素的参考态和能量。系统已根据您上传的结构自动识别组成元素，
          您可以根据需要调整参考态和能量值。
        </Paragraph>
        
        <Alert
          message="参考元素说明"
          description="参考元素的选择会影响形成能的计算结果。通常选择元素的标准态作为参考态，例如金属元素选择其稳定晶体结构，气体元素选择其分子形式。"
          type="info"
          showIcon
          style={{ marginBottom: 16 }}
        />
        
        <Card title="参考元素设置" style={{ marginBottom: 16 }}>
          <Table 
            dataSource={referenceElements} 
            columns={columns} 
            pagination={false}
            rowKey="element"
          />
        </Card>
        
        <Form layout="vertical">
          <Form.Item label="从材料数据库导入参考能量">
            <Button icon={<DatabaseOutlined />}>从Materials Project导入</Button>
          </Form.Item>
        </Form>
        
        <Divider />
        
        <Space>
          <Button onClick={prevStep}>上一步</Button>
          <Button type="primary" onClick={nextStep}>下一步</Button>
        </Space>
      </div>
    );
  };
  
  // 步骤3: 计算方法设置
  const renderCalculationMethod = () => {
    return (
      <div>
        <Title level={4}>计算方法设置</Title>
        <Paragraph>
          选择形成能计算方法和相关参数。
        </Paragraph>
        
        <Form form={form} layout="vertical">
          <Form.Item 
            name="calculationMethod" 
            label="计算方法" 
            initialValue="dft"
          >
            <Radio.Group onChange={(e) => setCalculationMethod(e.target.value)}>
              <Radio value="dft">密度泛函理论 (DFT)</Radio>
              <Radio value="gga">广义梯度近似 (GGA)</Radio>
              <Radio value="lda">局域密度近似 (LDA)</Radio>
              <Radio value="hybrid">杂化泛函</Radio>
            </Radio.Group>
          </Form.Item>
          
          {calculationMethod === 'dft' && (
            <>
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
                </Select>
              </Form.Item>
              
              <Form.Item 
                name="hubbardU" 
                label="Hubbard U校正"
                initialValue={false}
              >
                <Checkbox.Group>
                  <Checkbox value="Co">Co (U = 3.32 eV)</Checkbox>
                  <Checkbox value="Fe">Fe (U = 4.0 eV)</Checkbox>
                  <Checkbox value="Mn">Mn (U = 3.9 eV)</Checkbox>
                </Checkbox.Group>
              </Form.Item>
            </>
          )}
          
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
            name="spinPolarized" 
            label="自旋极化计算" 
            initialValue={true}
          >
            <Radio.Group>
              <Radio value={true}>是</Radio>
              <Radio value={false}>否</Radio>
            </Radio.Group>
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
  
  // 步骤4: 计算结果
  const renderResults = () => {
    return (
      <div>
        <Title level={4}>形成能计算结果</Title>
        
        {loading ? (
          <div style={{ textAlign: 'center', margin: '40px 0' }}>
            <Spin size="large" tip="正在计算形成能..." />
          </div>
        ) : (
          <>
            <Alert
              message="计算完成"
              description="形成能计算已成功完成，结果如下所示。"
              type="success"
              showIcon
              style={{ marginBottom: 24 }}
            />
            
            <Card title="计算结果概要" style={{ marginBottom: 24 }}>
              <Row gutter={[24, 24]}>
                <Col span={12}>
                  <Card type="inner" title="形成能">
                    <div style={{ textAlign: 'center', padding: '20px 0' }}>
                      <Title level={2} style={{ color: '#1C64F2' }}>{formationEnergy} eV/atom</Title>
                      <Text type="secondary">负值表示材料热力学稳定</Text>
                    </div>
                  </Card>
                </Col>
                <Col span={12}>
                  <Card type="inner" title="稳定性评估">
                    <div style={{ padding: '10px 0' }}>
                      <p><Text strong>稳定性:</Text> 热力学稳定</p>
                      <p><Text strong>相对稳定性:</Text> 高 (在类似材料中排名前20%)</p>
                      <p><Text strong>分解能:</Text> 0.12 eV/atom</p>
                      <p><Text strong>可能的分解产物:</Text> Li₂O, CoO</p>
                    </div>
                  </Card>
                </Col>
              </Row>
              
              <Divider />
              
              <Row gutter={[24, 24]}>
                <Col span={24}>
                  <Card type="inner" title="计算详情">
                    <p><Text strong>总能量:</Text> -45.67 eV</p>
                    <p><Text strong>每原子能量:</Text> -5.71 eV/atom</p>
                    <p><Text strong>参考元素总能量:</Text> -27.12 eV</p>
                    <p><Text strong>计算方法:</Text> DFT (PBE)</p>
                    <p><Text strong>收敛状态:</Text> 已收敛 (力收敛标准: 0.01 eV/Å)</p>
                  </Card>
                </Col>
              </Row>
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
  
  // 步骤5: 提交确认
  const renderSubmission = () => {
    return (
      <div style={{ textAlign: 'center', padding: '20px 0' }}>
        <CheckCircleOutlined style={{ fontSize: 72, color: '#52c41a', marginBottom: 24 }} />
        <Title level={3}>计算结果已保存！</Title>
        <Paragraph>
          您的形成能计算结果已成功保存。您可以在用户中心的"我的计算任务"中查看详细结果。
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
          <Step title="参考元素" icon={<DatabaseOutlined />} />
          <Step title="计算方法" icon={<SettingOutlined />} />
          <Step title="计算结果" icon={<CalculatorOutlined />} />
          <Step title="完成" icon={<CheckCircleOutlined />} />
        </Steps>
        
        {renderStepContent()}
      </Card>
    </div>
  );
};

export default FormationEnergyCalculation; 