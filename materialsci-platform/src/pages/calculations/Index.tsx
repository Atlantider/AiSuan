import React, { useState } from 'react';
import { Typography, Row, Col, Card, Button, Space } from 'antd';
import { useNavigate, useLocation } from 'react-router-dom';
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
  const location = useLocation();
  
  // 计算模块卡片
  const calculationModules = [
    {
      title: '电池材料计算',
      description: '电池电极材料、电解质和全电池系统的计算与分析',
      icon: <ThunderboltOutlined style={{ fontSize: '36px', color: '#1890ff' }} />,
      path: '/battery'
    },
    {
      title: '催化材料计算',
      description: '催化材料性能预测和反应机理分析',
      icon: <ExperimentOutlined style={{ fontSize: '36px', color: '#52c41a' }} />,
      path: '/catalysis'
    },
    {
      title: '半导体材料计算',
      description: '半导体材料能带结构和电子性质分析',
      icon: <DotChartOutlined style={{ fontSize: '36px', color: '#722ed1' }} />,
      path: '/semiconductor'
    },
    {
      title: '能源材料计算',
      description: '太阳能、热电等能源材料性能计算',
      icon: <LineChartOutlined style={{ fontSize: '36px', color: '#fa8c16' }} />,
      path: '/energy'
    },
    {
      title: '结构材料计算',
      description: '金属、陶瓷、高分子等结构材料的力学性能分析',
      icon: <ApartmentOutlined style={{ fontSize: '36px', color: '#eb2f96' }} />,
      path: '/structural'
    }
  ];

  return (
    <div style={{ padding: '24px' }}>
      <Title level={2}>计算模块</Title>
      <Paragraph>
        选择以下计算模块之一开始您的材料计算任务。任务完成后，您可以在"用户中心 &gt; 任务管理"中查看任务列表和结果。
      </Paragraph>

      <Row gutter={[24, 24]} style={{ marginTop: '24px' }}>
        {calculationModules.map((module, index) => (
          <Col xs={24} sm={12} lg={8} key={index}>
            <Card 
              hoverable 
              style={{ height: '100%' }}
              onClick={() => navigate(module.path)}
            >
              <Space direction="vertical" size="large" style={{ width: '100%', textAlign: 'center' }}>
                {module.icon}
                <div>
                  <Title level={4}>{module.title}</Title>
                  <Paragraph>{module.description}</Paragraph>
                </div>
                <Button type="primary">开始计算</Button>
              </Space>
            </Card>
          </Col>
        ))}
      </Row>
    </div>
  );
};

export default CalculationsIndex; 