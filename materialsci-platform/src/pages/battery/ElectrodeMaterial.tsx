import React from 'react';
import { Typography, Card, Row, Col, Button, List, Divider, Space, Tag } from 'antd';
import { useNavigate } from 'react-router-dom';
import { 
  ExperimentOutlined, 
  ThunderboltOutlined,
  LineChartOutlined,
  BuildOutlined
} from '@ant-design/icons';

const { Title, Paragraph } = Typography;

const ElectrodeMaterial: React.FC = () => {
  const navigate = useNavigate();

  // 可选计算类型
  const calculationTypes = [
    {
      title: '电极电位计算',
      description: '计算材料在充放电过程中的电极电位，预测电池的电压平台和能量密度',
      icon: <ThunderboltOutlined style={{ fontSize: 24, color: '#1C64F2' }} />,
      type: 'voltage'
    },
    {
      title: '离子迁移能垒计算',
      description: '计算离子在材料中迁移的能量障碍，预测电池的功率性能',
      icon: <LineChartOutlined style={{ fontSize: 24, color: '#1C64F2' }} />,
      type: 'barrier'
    },
    {
      title: '形成能计算',
      description: '计算材料的形成能，评估材料的热力学稳定性和合成可行性',
      icon: <BuildOutlined style={{ fontSize: 24, color: '#1C64F2' }} />,
      type: 'formation'
    },
    {
      title: '电子结构计算',
      description: '计算材料的能带结构、态密度等电子结构特性，理解材料的电子传输性质',
      icon: <ExperimentOutlined style={{ fontSize: 24, color: '#1C64F2' }} />,
      type: 'electronic'
    },
  ];

  // 处理计算类型选择
  const handleCalculationSelect = (calculationType: string) => {
    navigate(`/battery/electrode-material/workbench`, { state: { calculationType } });
  };

  return (
    <div>
      <div style={{ marginBottom: 24 }}>
        <Title level={2}>电极材料计算</Title>
        <Paragraph style={{ fontSize: 16 }}>
          电极材料计算模块提供电池电极材料性能的各种计算服务，包括电极电位、离子迁移能垒、形成能等多种计算，
          帮助研究人员深入了解电极材料的结构-性能关系，设计开发高性能电池材料。
        </Paragraph>
        <Space>
          <Button 
            type="primary" 
            size="large" 
            icon={<ExperimentOutlined />}
            onClick={() => navigate('/battery/electrode-material/workbench')}
          >
            进入综合计算工作台
          </Button>
          <Button 
            size="large"
            onClick={() => navigate('/documentation')}
          >
            查看文档
          </Button>
        </Space>
      </div>
      
      <Divider orientation="left">可选计算类型</Divider>
      
      <Row gutter={[24, 24]}>
        {calculationTypes.map((type, index) => (
          <Col xs={24} md={12} key={index}>
            <Card 
              hoverable 
              onClick={() => handleCalculationSelect(type.type)}
              style={{ cursor: 'pointer' }}
            >
              <div style={{ display: 'flex', alignItems: 'flex-start' }}>
                <div style={{ marginRight: 16 }}>
                  {type.icon}
                </div>
                <div>
                  <Title level={4}>{type.title}</Title>
                  <Paragraph>{type.description}</Paragraph>
                  <Button type="link" style={{ padding: 0 }}>开始计算 &gt;</Button>
                </div>
              </div>
            </Card>
          </Col>
        ))}
      </Row>
      
      <Divider orientation="left">典型应用案例</Divider>
      
      <List
        itemLayout="vertical"
        dataSource={[
          {
            title: 'LiFePO4 电极性能计算',
            description: '通过电极电位和Li离子迁移能垒计算，优化LiFePO4正极材料的电子结构和离子传输性能，提高其倍率性能。',
            tags: ['电极电位', '离子迁移', '电子结构'],
          },
          {
            title: 'Na离子电池层状氧化物正极筛选',
            description: '计算一系列层状氧化物的Na嵌入电位和结构稳定性，筛选高容量、高电压的钠离子电池正极材料。',
            tags: ['钠离子电池', '层状氧化物', '材料筛选'],
          },
        ]}
        renderItem={item => (
          <List.Item>
            <Card style={{ width: '100%' }}>
              <Title level={4}>{item.title}</Title>
              <div style={{ marginBottom: 12 }}>
                {item.tags.map(tag => (
                  <Tag color="blue" key={tag}>{tag}</Tag>
                ))}
              </div>
              <Paragraph>{item.description}</Paragraph>
            </Card>
          </List.Item>
        )}
      />
      
      <div style={{ marginTop: 40, textAlign: 'center' }}>
        <Button 
          type="primary" 
          size="large" 
          onClick={() => navigate('/battery/electrode-material/workbench')}
        >
          进入综合计算工作台
        </Button>
      </div>
    </div>
  );
};

export default ElectrodeMaterial; 