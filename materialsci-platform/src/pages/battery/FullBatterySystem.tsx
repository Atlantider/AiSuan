import React from 'react';
import { Typography, Card, Row, Col, Button, List, Divider, Space, Tag } from 'antd';
import { useNavigate } from 'react-router-dom';
import { 
  ExperimentOutlined, 
  RocketOutlined,
  SafetyCertificateOutlined,
  BarChartOutlined
} from '@ant-design/icons';

const { Title, Paragraph } = Typography;

const FullBatterySystem: React.FC = () => {
  const navigate = useNavigate();

  // 可选计算类型
  const calculationTypes = [
    {
      title: '电池性能预测',
      description: '基于电极材料和电解液特性，预测电池的容量、电压、功率等核心性能',
      icon: <BarChartOutlined style={{ fontSize: 24, color: '#1C64F2' }} />,
    },
    {
      title: '电极/电解质界面模拟',
      description: '模拟电极与电解质界面的相互作用，研究SEI膜形成和界面稳定性',
      icon: <ExperimentOutlined style={{ fontSize: 24, color: '#1C64F2' }} />,
    },
    {
      title: '循环寿命评估',
      description: '评估电池在不同循环条件下的容量衰减和寿命特性',
      icon: <RocketOutlined style={{ fontSize: 24, color: '#1C64F2' }} />,
    },
    {
      title: '安全性评估',
      description: '模拟高温、过充等极端条件下的电池反应，评估安全风险',
      icon: <SafetyCertificateOutlined style={{ fontSize: 24, color: '#1C64F2' }} />,
    },
  ];

  return (
    <div>
      <div style={{ marginBottom: 24 }}>
        <Title level={2}>全电池系统计算</Title>
        <Paragraph style={{ fontSize: 16 }}>
          全电池系统计算模块整合电极和电解液参数，进行完整电池系统的性能预测，帮助研究人员理解电池组件间的相互作用，
          优化电池设计，提高整体性能。
        </Paragraph>
        <Space>
          <Button 
            type="primary" 
            size="large" 
            icon={<ExperimentOutlined />}
            onClick={() => navigate('/battery/full-battery/workbench')}
          >
            进入工作台
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
            <Card hoverable>
              <div style={{ display: 'flex', alignItems: 'flex-start' }}>
                <div style={{ marginRight: 16 }}>
                  {type.icon}
                </div>
                <div>
                  <Title level={4}>{type.title}</Title>
                  <Paragraph>{type.description}</Paragraph>
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
            title: '固态锂金属电池优化',
            description: '通过组合模拟电极材料、固态电解质和界面特性，优化固态锂金属电池的设计，提高能量密度和安全性。',
            tags: ['固态电池', '锂金属负极', '界面稳定性'],
          },
          {
            title: '高温工作电池设计',
            description: '模拟评估不同电极-电解液组合在高温环境下的稳定性和性能，设计适用于极端环境的电池系统。',
            tags: ['高温环境', '热稳定性', '安全评估'],
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
          onClick={() => navigate('/battery/full-battery/workbench')}
        >
          开始计算
        </Button>
      </div>
    </div>
  );
};

export default FullBatterySystem; 