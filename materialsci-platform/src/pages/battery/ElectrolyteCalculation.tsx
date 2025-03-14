import React, { useState } from 'react';
import { Typography, Card, Row, Col, Button, List, Divider, Space, Tag, Tabs, Steps, Image, Form, Input, InputNumber, Select, Upload, message, Tooltip, Checkbox, theme } from 'antd';
import { useNavigate } from 'react-router-dom';
import { 
  ExperimentOutlined, 
  NodeIndexOutlined,
  BranchesOutlined,
  InteractionOutlined,
  SolutionOutlined,
  AppstoreOutlined,
  BarChartOutlined,
  ClusterOutlined,
  SettingOutlined,
  UploadOutlined,
  DownloadOutlined,
  PlusOutlined,
  MinusCircleOutlined,
  InfoCircleOutlined
} from '@ant-design/icons';

const { Title, Paragraph, Text } = Typography;
const { TabPane } = Tabs;
const { Option } = Select;
const { Step } = Steps;
const { useToken } = theme;

const ElectrolyteCalculation: React.FC = () => {
  const navigate = useNavigate();
  const [activeTab, setActiveTab] = useState('1');
  const [form] = Form.useForm();
  const [currentStep, setCurrentStep] = useState(0);
  const { token } = useToken();

  // 阳离子选项
  const cationOptions = [
    { label: 'Li⁺', value: 'Li', charge: 1 },
    { label: 'Na⁺', value: 'Na', charge: 1 },
    { label: 'K⁺', value: 'K', charge: 1 },
    { label: 'Mg²⁺', value: 'Mg', charge: 2 },
    { label: 'Ca²⁺', value: 'Ca', charge: 2 },
  ];

  // 阴离子选项
  const anionOptions = [
    { label: 'PF₆⁻', value: 'PF6', charge: -1 },
    { label: 'BF₄⁻', value: 'BF4', charge: -1 },
    { label: 'TFSI⁻', value: 'TFSI', charge: -1 },
    { label: 'FSI⁻', value: 'FSI', charge: -1 },
  ];

  // 溶剂选项
  const solventOptions = [
    { 
      label: 'EC (碳酸乙烯酯)', 
      value: 'EC',
      smile: 'C1COC(=O)O1',
      description: '环状碳酸酯'
    },
    { 
      label: 'DMC (碳酸二甲酯)', 
      value: 'DMC',
      smile: 'COC(=O)OC',
      description: '链状碳酸酯'
    },
    { 
      label: 'EMC (碳酸甲乙酯)', 
      value: 'EMC',
      smile: 'CCOC(=O)OC',
      description: '链状碳酸酯'
    },
    { 
      label: 'PC (碳酸丙烯酯)', 
      value: 'PC',
      smile: 'CC1COC(=O)O1',
      description: '环状碳酸酯'
    },
    { 
      label: 'DEC (碳酸二乙酯)', 
      value: 'DEC',
      smile: 'CCOC(=O)OCC',
      description: '链状碳酸酯'
    },
    { 
      label: 'FEC (氟代碳酸乙烯酯)', 
      value: 'FEC',
      smile: 'C1OC(=O)OC1F',
      description: '氟代环状碳酸酯'
    },
    {
      label: '自定义溶剂',
      value: 'custom',
      smile: '',
      description: '自定义溶剂类型'
    }
  ];

  // 计算类型选项
  const calculationOptions = [
    {
      label: '单分子性质计算',
      value: 'single_molecule',
      description: '计算电解液组分分子的几何结构、能量、前线轨道分布等基本性质',
      icon: <NodeIndexOutlined style={{ fontSize: 24, color: token.colorPrimary }} />,
      requires: ['geometry_optimization', 'energy_calculation', 'orbital_analysis']
    },
    {
      label: '离子溶剂化结构计算',
      value: 'solvation',
      description: '计算锂离子与溶剂分子和阴离子的溶剂化结构和结合能',
      icon: <BranchesOutlined style={{ fontSize: 24, color: token.colorPrimary }} />,
      requires: ['rdf_analysis', 'coordination_analysis', 'binding_energy']
    },
    {
      label: '电解液传输性质计算',
      value: 'transport',
      description: '计算电解液体系的离子电导率、扩散系数等传输性质',
      icon: <InteractionOutlined style={{ fontSize: 24, color: token.colorPrimary }} />,
      requires: ['msd_analysis', 'conductivity_calculation', 'viscosity_calculation']
    },
    {
      label: '界面反应与SEI膜形成计算',
      value: 'interface',
      description: '模拟电解液在电极表面的分解反应和SEI膜形成过程',
      icon: <ExperimentOutlined style={{ fontSize: 24, color: token.colorPrimary }} />,
      requires: ['interface_reaction', 'decomposition_analysis', 'sei_formation']
    }
  ];

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
  const handleSubmit = (values: any) => {
    console.log('Form values:', values);
    // 根据选择的计算类型生成不同的LAMMPS输入文件
    generateLammpsInputs(values);
    // 显示成功消息
    message.success('已成功生成LAMMPS输入文件！');
  };

  // 生成LAMMPS输入文件
  const generateLammpsInputs = (values: any) => {
    // 根据选择的计算类型和参数生成相应的LAMMPS输入文件
    const { calculations, temperatures, cations, anions, solvents, concentration } = values;
    console.log({ calculations, temperatures, cations, anions, solvents, concentration });
    
    // 这里添加生成LAMMPS输入文件的逻辑
    // ...
  };

  // 渲染步骤内容
  const renderStepContent = () => {
    switch (currentStep) {
      case 0:
        return renderFormulaDesign();
      case 1:
        return renderConditionsSettings();
      case 2:
        return renderCalculationSelection();
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
          style={{ marginBottom: 24, boxShadow: '0 1px 2px rgba(0,0,0,0.05)' }}
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
          <Row gutter={[24, 24]}>
            <Col span={12}>
              <Card 
                type="inner" 
                title={<span style={{ fontWeight: 'bold' }}>阳离子选择</span>} 
                style={{ height: '100%' }}
              >
                <Form.List name="cations">
                  {(fields, { add, remove }) => (
                    <>
                      {fields.map(({ key, name, ...restField }) => (
                        <Row key={key} gutter={16} align="middle" style={{ marginBottom: 16 }}>
                          <Col span={12}>
                            <Form.Item
                              {...restField}
                              name={[name, 'type']}
                              rules={[{ required: true, message: '请选择阳离子' }]}
                            >
                              <Select 
                                placeholder="选择阳离子"
                                showSearch
                                optionFilterProp="children"
                              >
                                {cationOptions.map(option => (
                                  <Option key={option.value} value={option.value}>
                                    <Tooltip title={`电荷: ${option.charge}+`}>
                                      {option.label}
                                    </Tooltip>
                                  </Option>
                                ))}
                              </Select>
                            </Form.Item>
                          </Col>
                          <Col span={10}>
                            <Form.Item
                              {...restField}
                              name={[name, 'ratio']}
                              rules={[{ required: true, message: '请输入比例' }]}
                            >
                              <InputNumber 
                                min={0.01} 
                                max={10} 
                                step={0.01} 
                                placeholder="比例" 
                                style={{ width: '100%' }}
                              />
                            </Form.Item>
                          </Col>
                          <Col span={2}>
                            <Button 
                              type="text" 
                              danger 
                              icon={<MinusCircleOutlined />} 
                              onClick={() => remove(name)}
                            />
                          </Col>
                        </Row>
                      ))}
                      <Form.Item>
                        <Button
                          type="dashed"
                          onClick={() => add()}
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
            </Col>
            
            <Col span={12}>
              <Card 
                type="inner" 
                title={<span style={{ fontWeight: 'bold' }}>阴离子选择</span>} 
                style={{ height: '100%' }}
              >
                <Form.List name="anions">
                  {(fields, { add, remove }) => (
                    <>
                      {fields.map(({ key, name, ...restField }) => (
                        <Row key={key} gutter={16} align="middle" style={{ marginBottom: 16 }}>
                          <Col span={12}>
                            <Form.Item
                              {...restField}
                              name={[name, 'type']}
                              rules={[{ required: true, message: '请选择阴离子' }]}
                            >
                              <Select 
                                placeholder="选择阴离子"
                                showSearch
                                optionFilterProp="children"
                              >
                                {anionOptions.map(option => (
                                  <Option key={option.value} value={option.value}>
                                    <Tooltip title={`电荷: ${option.charge}`}>
                                      {option.label}
                                    </Tooltip>
                                  </Option>
                                ))}
                              </Select>
                            </Form.Item>
                          </Col>
                          <Col span={10}>
                            <Form.Item
                              {...restField}
                              name={[name, 'ratio']}
                              rules={[{ required: true, message: '请输入比例' }]}
                            >
                              <InputNumber 
                                min={0.01} 
                                max={10} 
                                step={0.01} 
                                placeholder="比例" 
                                style={{ width: '100%' }}
                              />
                            </Form.Item>
                          </Col>
                          <Col span={2}>
                            <Button 
                              type="text" 
                              danger 
                              icon={<MinusCircleOutlined />} 
                              onClick={() => remove(name)}
                            />
                          </Col>
                        </Row>
                      ))}
                      <Form.Item>
                        <Button
                          type="dashed"
                          onClick={() => add()}
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
            </Col>
          </Row>

          <Divider />

          <Card 
            type="inner" 
            title={<span style={{ fontWeight: 'bold' }}>溶剂组成</span>}
            style={{ marginTop: 16 }}
          >
            <Form.List name="solvents">
              {(fields, { add, remove }) => (
                <>
                  {fields.map(({ key, name, ...restField }) => (
                    <Row key={key} gutter={16} align="middle" style={{ marginBottom: 16 }}>
                      <Col span={10}>
                        <Form.Item
                          {...restField}
                          name={[name, 'type']}
                          rules={[{ required: true, message: '请选择溶剂' }]}
                        >
                          <Select 
                            placeholder="选择溶剂"
                            showSearch
                            optionFilterProp="children"
                            onChange={(value) => {
                              if (value === 'custom') {
                                // 如果选择自定义，清空SMILE字段让用户输入
                                form.setFieldValue(['solvents', name, 'smile'], '');
                                // 设置自定义标记为true
                                form.setFieldValue(['solvents', name, 'isCustom'], true);
                              } else {
                                const solvent = solventOptions.find(s => s.value === value);
                                if (solvent) {
                                  form.setFieldValue(['solvents', name, 'smile'], solvent.smile);
                                  // 设置自定义标记为false
                                  form.setFieldValue(['solvents', name, 'isCustom'], false);
                                }
                              }
                            }}
                          >
                            {solventOptions.map(option => (
                              <Option key={option.value} value={option.value}>
                                <Tooltip title={option.description}>
                                  {option.label}
                                </Tooltip>
                              </Option>
                            ))}
                          </Select>
                        </Form.Item>
                        {form.getFieldValue(['solvents', name, 'isCustom']) && (
                          <Form.Item
                            {...restField}
                            name={[name, 'customName']}
                            rules={[{ required: true, message: '请输入溶剂名称' }]}
                            style={{ marginTop: 8 }}
                          >
                            <Input placeholder="自定义溶剂名称" />
                          </Form.Item>
                        )}
                      </Col>
                      <Col span={8}>
                        <Form.Item
                          {...restField}
                          name={[name, 'smile']}
                          rules={[
                            { 
                              required: true, 
                              message: '请输入SMILE结构' 
                            }
                          ]}
                        >
                          <Input 
                            placeholder="SMILE结构" 
                            readOnly={!form.getFieldValue(['solvents', name, 'isCustom'])}
                            style={{ 
                              backgroundColor: form.getFieldValue(['solvents', name, 'isCustom']) 
                                ? 'white' 
                                : '#f5f5f5' 
                            }}
                          />
                        </Form.Item>
                      </Col>
                      <Col span={4}>
                        <Form.Item
                          {...restField}
                          name={[name, 'ratio']}
                          rules={[{ required: true, message: '请输入比例' }]}
                        >
                          <InputNumber 
                            min={0.01} 
                            max={10} 
                            step={0.01} 
                            placeholder="比例" 
                            style={{ width: '100%' }}
                          />
                        </Form.Item>
                      </Col>
                      <Col span={2}>
                        <Button 
                          type="text" 
                          danger 
                          icon={<MinusCircleOutlined />} 
                          onClick={() => remove(name)}
                        />
                      </Col>
                    </Row>
                  ))}
                  <Form.Item>
                    <Button
                      type="dashed"
                      onClick={() => add()}
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
        </Card>
      </>
    );
  };

  // 渲染条件设置（温度和摩尔浓度）
  const renderConditionsSettings = () => {
    return (
      <Card 
        title={<Title level={4}>模拟条件设置</Title>} 
        style={{ marginBottom: 24, boxShadow: '0 1px 2px rgba(0,0,0,0.05)' }}
      >
        <Row gutter={[24, 24]}>
          <Col span={12}>
            <Card 
              type="inner" 
              title={<span style={{ fontWeight: 'bold' }}>温度设置</span>}
            >
              <Form.List name="temperatures">
                {(fields, { add, remove }) => (
                  <>
                    {fields.map(({ key, name, ...restField }) => (
                      <Row key={key} gutter={16} align="middle" style={{ marginBottom: 16 }}>
                        <Col span={20}>
                          <Form.Item
                            {...restField}
                            name={[name, 'value']}
                            rules={[{ required: true, message: '请输入温度' }]}
                          >
                            <InputNumber 
                              min={-40} 
                              max={80} 
                              placeholder="温度" 
                              style={{ width: '100%' }}
                              addonAfter="°C"
                            />
                          </Form.Item>
                        </Col>
                        <Col span={4}>
                          <Button 
                            type="text" 
                            danger 
                            icon={<MinusCircleOutlined />} 
                            onClick={() => remove(name)}
                          />
                        </Col>
                      </Row>
                    ))}
                    <Form.Item>
                      <Button
                        type="dashed"
                        onClick={() => add()}
                        block
                        icon={<PlusOutlined />}
                      >
                        添加温度点
                      </Button>
                    </Form.Item>
                  </>
                )}
              </Form.List>
            </Card>
          </Col>
          
          <Col span={12}>
            <Card 
              type="inner" 
              title={
                <Space>
                  <span style={{ fontWeight: 'bold' }}>摩尔浓度设置</span>
                  <Tooltip title="设置电解液中盐的摩尔浓度，单位为mol/L">
                    <InfoCircleOutlined />
                  </Tooltip>
                </Space>
              }
            >
              <Form.List name="concentrations">
                {(fields, { add, remove }) => (
                  <>
                    {fields.map(({ key, name, ...restField }) => (
                      <Row key={key} gutter={16} align="middle" style={{ marginBottom: 16 }}>
                        <Col span={20}>
                          <Form.Item
                            {...restField}
                            name={[name, 'value']}
                            rules={[{ required: true, message: '请输入摩尔浓度' }]}
                          >
                            <InputNumber 
                              min={0.1} 
                              max={5.0} 
                              step={0.1} 
                              placeholder="盐浓度" 
                              style={{ width: '100%' }}
                              addonAfter="mol/L"
                            />
                          </Form.Item>
                        </Col>
                        <Col span={4}>
                          <Button 
                            type="text" 
                            danger 
                            icon={<MinusCircleOutlined />} 
                            onClick={() => remove(name)}
                          />
                        </Col>
                      </Row>
                    ))}
                    <Form.Item>
                      <Button
                        type="dashed"
                        onClick={() => add()}
                        block
                        icon={<PlusOutlined />}
                      >
                        添加浓度点
                      </Button>
                    </Form.Item>
                  </>
                )}
              </Form.List>
            </Card>
          </Col>
        </Row>

        <Card 
          type="inner" 
          title={<span style={{ fontWeight: 'bold' }}>其它模拟参数</span>}
          style={{ marginTop: 24 }}
        >
          <Row gutter={[24, 16]}>
            <Col span={8}>
              <Form.Item
                name="simulationTime"
                label="模拟时间"
                rules={[{ required: true, message: '请输入模拟时间' }]}
              >
                <InputNumber 
                  min={1} 
                  max={1000} 
                  placeholder="模拟时间" 
                  style={{ width: '100%' }}
                  addonAfter="ns"
                />
              </Form.Item>
            </Col>
            <Col span={8}>
              <Form.Item
                name="timeStep"
                label="时间步长"
                rules={[{ required: true, message: '请输入时间步长' }]}
                initialValue={1}
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
                name="equilibrationTime"
                label="平衡时间"
                rules={[{ required: true, message: '请输入平衡时间' }]}
              >
                <InputNumber 
                  min={0} 
                  max={100} 
                  placeholder="平衡时间" 
                  style={{ width: '100%' }}
                  addonAfter="ns"
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
        style={{ marginBottom: 24, boxShadow: '0 1px 2px rgba(0,0,0,0.05)' }}
      >
        <Paragraph style={{ marginBottom: 16 }}>
          请选择需要进行的计算类型，系统将根据选择生成相应的LAMMPS输入文件。
        </Paragraph>
        
        <Form.Item 
          name="calculations" 
          rules={[{ required: true, message: '请选择至少一种计算类型' }]}
        >
          <Checkbox.Group style={{ width: '100%' }}>
            <Row gutter={[24, 24]}>
              {calculationOptions.map(option => (
                <Col span={12} key={option.value}>
                  <Card 
                    hoverable 
                    style={{ height: '100%', transition: 'all 0.3s' }}
                  >
                    <div style={{ display: 'flex', alignItems: 'flex-start' }}>
                      <div style={{ marginRight: 16 }}>
                        {option.icon}
                      </div>
                      <div>
                        <Checkbox value={option.value}>
                          <Title level={5}>{option.label}</Title>
                        </Checkbox>
                        <Paragraph style={{ marginTop: 8 }}>
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

        <Card 
          type="inner" 
          title={<span style={{ fontWeight: 'bold' }}>输出设置</span>}
          style={{ marginTop: 16 }}
        >
          <Row gutter={[24, 16]}>
            <Col span={12}>
              <Form.Item
                name="outputFrequency"
                label="输出频率"
                initialValue={1000}
              >
                <InputNumber 
                  min={100} 
                  max={10000} 
                  placeholder="输出频率" 
                  style={{ width: '100%' }}
                  addonAfter="步"
                />
              </Form.Item>
            </Col>
            <Col span={12}>
              <Form.Item
                name="outputFolder"
                label="输出文件夹名称"
              >
                <Input placeholder="默认使用日期作为文件夹名" />
              </Form.Item>
            </Col>
          </Row>
        </Card>
      </Card>
    );
  };

  return (
    <div style={{ maxWidth: 1200, margin: '0 auto', padding: '24px' }}>
      <div style={{ marginBottom: 24 }}>
        <Title level={2} style={{ marginBottom: 8 }}>电解液计算工作台</Title>
        <Paragraph style={{ fontSize: 16, color: token.colorTextSecondary }}>
          设计电解液配方，选择计算类型，生成LAMMPS输入文件进行分子动力学模拟。
        </Paragraph>
      </div>

      <Card
        style={{ 
          marginBottom: 24, 
          boxShadow: '0 2px 8px rgba(0,0,0,0.08)',
          borderRadius: token.borderRadiusLG 
        }}
      >
        <Steps current={currentStep} style={{ marginBottom: 24, padding: '0 16px' }}>
          <Step title="配方设计" description="阴阳离子和溶剂" />
          <Step title="条件设置" description="温度和摩尔浓度" />
          <Step title="计算选择" description="选择计算类型" />
        </Steps>
      </Card>

      <Form
        form={form}
        layout="vertical"
        onFinish={handleSubmit}
        initialValues={{
          cations: [{}],
          anions: [{}],
          solvents: [{}],
          temperatures: [{}],
          concentrations: [{ value: 1.0 }],
          simulationTime: 10,
          timeStep: 1,
          equilibrationTime: 2,
          outputFrequency: 1000
        }}
      >
        {renderStepContent()}

        <div style={{ marginTop: 24, textAlign: 'right' }}>
          {currentStep > 0 && (
            <Button 
              style={{ marginRight: 8 }} 
              onClick={() => setCurrentStep(currentStep - 1)}
            >
              上一步
            </Button>
          )}
          {currentStep < 2 && (
            <Button 
              type="primary" 
              onClick={() => setCurrentStep(currentStep + 1)}
            >
              下一步
            </Button>
          )}
          {currentStep === 2 && (
            <Button 
              type="primary" 
              size="large"
              icon={<SettingOutlined />}
              onClick={() => form.submit()}
            >
              生成输入文件
            </Button>
          )}
        </div>
      </Form>
    </div>
  );
};

export default ElectrolyteCalculation; 