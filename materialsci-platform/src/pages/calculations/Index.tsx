import React from 'react';
import { Typography, Row, Col, Card, Button, Space } from 'antd';
import { useNavigate } from 'react-router-dom';
import { 
  ThunderboltOutlined, 
  ExperimentOutlined, 
  DotChartOutlined,
  LineChartOutlined, 
  ApartmentOutlined
} from '@ant-design/icons';

const { Title, Paragraph } = Typography;

const CalculationsIndex: React.FC = () => {
  const navigate = useNavigate();

  // 计算模块列表
  const calculationModules = [
    {
      title: '电池材料计算',
      description: '提供电极材料、电解液和全电池系统的计算功能，帮助设计和优化高性能电池材料',
      icon: <ThunderboltOutlined style={{ fontSize: 48, color: '#1C64F2' }} />,
      path: '/battery'
    },
    {
      title: '催化材料计算',
      description: '支持表面催化、电催化和光催化反应机理与性能评估，加速催化剂筛选和设计',
      icon: <ExperimentOutlined style={{ fontSize: 48, color: '#1C64F2' }} />,
      path: '/catalysis'
    },
    {
      title: '半导体材料计算',
      description: '计算半导体材料的电子结构、光学性质和载流子动力学，辅助设计新型电子和光电材料',
      icon: <DotChartOutlined style={{ fontSize: 48, color: '#1C64F2' }} />,
      path: '/semiconductor'
    },
    {
      title: '能源材料计算',
      description: '专注于太阳能、热电和储能材料的性能预测和优化，促进可再生能源技术发展',
      icon: <LineChartOutlined style={{ fontSize: 48, color: '#1C64F2' }} />,
      path: '/energy'
    },
    {
      title: '结构材料计算',
      description: '分析结构材料的力学性能、热稳定性和失效机制，支持先进结构材料的开发',
      icon: <ApartmentOutlined style={{ fontSize: 48, color: '#1C64F2' }} />,
      path: '/structural'
    }
  ];

  return (
    <div>
      <div style={{ marginBottom: 40 }}>
        <Title level={2}>选择计算类型</Title>
        <Paragraph>
          选择您要进行的计算类型，平台将为您提供专业的计算工具和分析能力。
          我们的计算模块涵盖材料科学的多个领域，帮助您高效开展研究工作。
        </Paragraph>
      </div>

      <Row gutter={[24, 24]}>
        {calculationModules.map((module, index) => (
          <Col xs={24} sm={12} md={8} key={index}>
            <Card 
              hoverable 
              style={{ height: '100%', display: 'flex', flexDirection: 'column' }}
              bodyStyle={{ flex: 1, display: 'flex', flexDirection: 'column' }}
            >
              <div style={{ textAlign: 'center', marginBottom: 20 }}>
                {module.icon}
              </div>
              <Title level={3} style={{ textAlign: 'center', fontSize: 20 }}>
                {module.title}
              </Title>
              <Paragraph style={{ flexGrow: 1 }}>
                {module.description}
              </Paragraph>
              <Button 
                type="primary" 
                block 
                size="large"
                onClick={() => navigate(module.path)}
              >
                进入模块
              </Button>
            </Card>
          </Col>
        ))}
      </Row>

      <div style={{ marginTop: 40, textAlign: 'center' }}>
        <Space>
          <Button onClick={() => navigate('/user/dashboard')}>返回控制台</Button>
          <Button type="primary" onClick={() => navigate('/documentation')}>查看使用文档</Button>
        </Space>
      </div>
    </div>
  );
};

export default CalculationsIndex; 