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
  Tooltip
} from 'antd';
import { 
  UploadOutlined, 
  InboxOutlined, 
  ThunderboltOutlined, 
  SettingOutlined, 
  LineChartOutlined, 
  CheckCircleOutlined,
  EnvironmentOutlined,
  ArrowsAltOutlined
} from '@ant-design/icons';
import type { UploadProps } from 'antd';

const { Title, Paragraph, Text } = Typography;
const { Step } = Steps;
const { Option } = Select;
const { Dragger } = Upload;

interface MigrationPath {
  id: number;
  startSite: string;
  endSite: string;
  distance: number;
  selected: boolean;
}

const IonMigrationBarrierCalculation: React.FC = () => {
  const [currentStep, setCurrentStep] = useState(0);
  const [form] = Form.useForm();
  const [fileList, setFileList] = useState<any[]>([]);
  const [structurePreview, setStructurePreview] = useState<boolean>(false);
  const [loading, setLoading] = useState<boolean>(false);
  const [ionType, setIonType] = useState<string>('Li');
  const [migrationPaths, setMigrationPaths] = useState<MigrationPath[]>([
    { id: 1, startSite: 'Li1', endSite: 'Li2', distance: 2.85, selected: true },
    { id: 2, startSite: 'Li2', endSite: 'Li3', distance: 3.12, selected: false },
    { id: 3, startSite: 'Li3', endSite: 'Li4', distance: 2.97, selected: false },
  ]);
  const [nebImages, setNebImages] = useState<number>(7);
  const [barrierHeight, setBarrierHeight] = useState<number | null>(null);
  
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
            setBarrierHeight(0.38);
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
        return renderIonSelection();
      case 2:
        return renderPathSelection();
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
          请上传需要计算离子迁移能垒的材料晶体结构文件。支持的格式包括CIF、POSCAR、PDB、XYZ和CML。
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
  
  // 步骤2: 离子选择
  const renderIonSelection = () => {
    return (
      <div>
        <Title level={4}>迁移离子设置</Title>
        <Paragraph>
          选择需要计算迁移能垒的离子类型和相关参数。
        </Paragraph>
        
        <Form form={form} layout="vertical">
          <Form.Item 
            name="ionType" 
            label="迁移离子类型" 
            initialValue="Li"
            rules={[{ required: true, message: '请选择迁移离子类型' }]}
          >
            <Radio.Group onChange={(e) => setIonType(e.target.value)}>
              <Radio value="Li">Li⁺</Radio>
              <Radio value="Na">Na⁺</Radio>
              <Radio value="K">K⁺</Radio>
              <Radio value="Mg">Mg²⁺</Radio>
              <Radio value="Ca">Ca²⁺</Radio>
              <Radio value="Al">Al³⁺</Radio>
            </Radio.Group>
          </Form.Item>
          
          <Form.Item 
            name="pathIdentificationMethod" 
            label="迁移路径识别方法" 
            initialValue="auto"
          >
            <Radio.Group>
              <Radio value="auto">自动识别</Radio>
              <Radio value="manual">手动设置</Radio>
            </Radio.Group>
          </Form.Item>
          
          <Form.Item 
            name="supercellSize" 
            label="超胞大小" 
            initialValue="2x2x1"
          >
            <Select style={{ width: 200 }}>
              <Option value="1x1x1">1 × 1 × 1</Option>
              <Option value="2x2x1">2 × 2 × 1</Option>
              <Option value="2x2x2">2 × 2 × 2</Option>
              <Option value="3x3x1">3 × 3 × 1</Option>
              <Option value="3x3x3">3 × 3 × 3</Option>
            </Select>
          </Form.Item>
          
          <Form.Item 
            name="nebImages" 
            label={
              <Tooltip title="NEB计算中的图像数量，数量越多计算越精确但耗时更长">
                NEB图像数量
              </Tooltip>
            } 
            initialValue={7}
          >
            <Row gutter={16}>
              <Col span={16}>
                <Slider
                  min={3}
                  max={11}
                  step={2}
                  value={nebImages}
                  onChange={(value) => setNebImages(value)}
                />
              </Col>
              <Col span={8}>
                <InputNumber
                  min={3}
                  max={11}
                  step={2}
                  value={nebImages}
                  onChange={(value) => setNebImages(value || 7)}
                />
              </Col>
            </Row>
          </Form.Item>
          
          <Form.Item 
            name="climbingImage" 
            label="使用爬坡图像(CI-NEB)" 
            initialValue={true}
          >
            <Radio.Group>
              <Radio value={true}>是</Radio>
              <Radio value={false}>否</Radio>
            </Radio.Group>
          </Form.Item>
        </Form>
        
        <Alert
          message="NEB计算说明"
          description="弹性带方法(Nudged Elastic Band, NEB)是计算离子迁移能垒的常用方法。通过在初始态和末态之间插入一系列图像，并优化这些图像的能量，可以找到最小能量路径和过渡态能量。"
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
  
  // 步骤3: 路径选择
  const renderPathSelection = () => {
    return (
      <div>
        <Title level={4}>迁移路径选择</Title>
        <Paragraph>
          系统已自动识别可能的离子迁移路径。请选择需要计算能垒的路径。
        </Paragraph>
        
        <Card title="可能的迁移路径" style={{ marginBottom: 16 }}>
          <div style={{ height: 200, background: '#f0f2f5', display: 'flex', justifyContent: 'center', alignItems: 'center', marginBottom: 16 }}>
            <Text type="secondary">结构预览区域（显示迁移路径）</Text>
          </div>
          
          <div>
            {migrationPaths.map(path => (
              <Card 
                key={path.id} 
                size="small" 
                style={{ 
                  marginBottom: 8, 
                  borderColor: path.selected ? '#1C64F2' : undefined,
                  backgroundColor: path.selected ? '#f0f5ff' : undefined
                }}
                onClick={() => {
                  const newPaths = migrationPaths.map(p => ({
                    ...p,
                    selected: p.id === path.id
                  }));
                  setMigrationPaths(newPaths);
                }}
              >
                <Row align="middle">
                  <Col span={1}>
                    <Radio checked={path.selected} />
                  </Col>
                  <Col span={18}>
                    <Space>
                      <Text strong>路径 {path.id}:</Text>
                      <Text>{path.startSite}</Text>
                      <ArrowsAltOutlined />
                      <Text>{path.endSite}</Text>
                      <Text type="secondary">({path.distance.toFixed(2)} Å)</Text>
                    </Space>
                  </Col>
                  <Col span={5} style={{ textAlign: 'right' }}>
                    <Button 
                      type="link" 
                      icon={<EnvironmentOutlined />}
                      onClick={(e) => {
                        e.stopPropagation();
                        message.info(`在结构中显示路径 ${path.id}`);
                      }}
                    >
                      在结构中显示
                    </Button>
                  </Col>
                </Row>
              </Card>
            ))}
          </div>
        </Card>
        
        <Form layout="vertical">
          <Form.Item 
            label="计算精度" 
            initialValue="standard"
          >
            <Radio.Group>
              <Radio value="fast">快速（精度低）</Radio>
              <Radio value="standard">标准</Radio>
              <Radio value="accurate">精确（耗时长）</Radio>
            </Radio.Group>
          </Form.Item>
          
          <Row gutter={24}>
            <Col span={12}>
              <Form.Item 
                label="计算引擎" 
                initialValue="vasp"
              >
                <Select>
                  <Option value="vasp">VASP</Option>
                  <Option value="qe">Quantum ESPRESSO</Option>
                  <Option value="cp2k">CP2K</Option>
                </Select>
              </Form.Item>
            </Col>
            <Col span={12}>
              <Form.Item 
                label="并行核心数" 
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
        <Title level={4}>离子迁移能垒计算结果</Title>
        
        {loading ? (
          <div style={{ textAlign: 'center', margin: '40px 0' }}>
            <Spin size="large" tip="正在计算离子迁移能垒..." />
          </div>
        ) : (
          <>
            <Alert
              message="计算完成"
              description="离子迁移能垒计算已成功完成，结果如下所示。"
              type="success"
              showIcon
              style={{ marginBottom: 24 }}
            />
            
            <Card title="计算结果概要" style={{ marginBottom: 24 }}>
              <Row gutter={[24, 24]}>
                <Col span={12}>
                  <Card type="inner" title="能垒高度">
                    <div style={{ textAlign: 'center', padding: '20px 0' }}>
                      <Title level={2} style={{ color: '#1C64F2' }}>{barrierHeight} eV</Title>
                      <Text type="secondary">能垒越低，离子迁移越容易</Text>
                    </div>
                  </Card>
                </Col>
                <Col span={12}>
                  <Card type="inner" title="扩散性能评估">
                    <div style={{ padding: '10px 0' }}>
                      <p><Text strong>扩散系数(300K):</Text> 1.2 × 10⁻⁹ cm²/s</p>
                      <p><Text strong>离子电导率:</Text> 2.5 mS/cm</p>
                      <p><Text strong>相对性能:</Text> 优 (在类似材料中排名前30%)</p>
                      <p><Text strong>适用温度范围:</Text> 250K - 400K</p>
                    </div>
                  </Card>
                </Col>
              </Row>
              
              <Divider />
              
              <div style={{ marginBottom: 16 }}>
                <Title level={5}>能量-反应坐标曲线</Title>
                <div style={{ height: 300, background: '#f0f2f5', display: 'flex', justifyContent: 'center', alignItems: 'center' }}>
                  <Text type="secondary">能量-反应坐标曲线图表区域</Text>
                </div>
              </div>
              
              <Row gutter={[24, 24]}>
                <Col span={24}>
                  <Card type="inner" title="计算详情">
                    <p><Text strong>迁移离子:</Text> {ionType}⁺</p>
                    <p><Text strong>迁移路径:</Text> {migrationPaths.find(p => p.selected)?.startSite} → {migrationPaths.find(p => p.selected)?.endSite} ({migrationPaths.find(p => p.selected)?.distance.toFixed(2)} Å)</p>
                    <p><Text strong>NEB图像数:</Text> {nebImages}</p>
                    <p><Text strong>计算方法:</Text> CI-NEB (VASP)</p>
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
          您的离子迁移能垒计算结果已成功保存。您可以在用户中心的"我的计算任务"中查看详细结果。
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
          <Step title="路径选择" icon={<EnvironmentOutlined />} />
          <Step title="计算结果" icon={<LineChartOutlined />} />
          <Step title="完成" icon={<CheckCircleOutlined />} />
        </Steps>
        
        {renderStepContent()}
      </Card>
    </div>
  );
};

export default IonMigrationBarrierCalculation; 