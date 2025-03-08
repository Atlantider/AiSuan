import React from 'react';
import { Typography, Card, Row, Col, Button, List, Divider, Space, Tag } from 'antd';
import { useNavigate } from 'react-router-dom';
import { 
  ExperimentOutlined, 
  NodeIndexOutlined,
  BranchesOutlined,
  InteractionOutlined,
} from '@ant-design/icons';

const { Title, Paragraph } = Typography;

const ElectrolyteCalculation: React.FC = () => {
  const navigate = useNavigate();

  // 可选计算类型
  const calculationTypes = [
    {
      title: '单分子性质计算',
      description: '计算电解液组分分子的几何结构、能量、前线轨道分布等基本性质',
      icon: <NodeIndexOutlined style={{ fontSize: 24, color: '#1C64F2' }} />,
    },
    {
      title: '离子溶剂化结构计算',
      description: '计算锂离子与溶剂分子和阴离子的溶剂化结构和结合能',
      icon: <BranchesOutlined style={{ fontSize: 24, color: '#1C64F2' }} />,
    },
    {
      title: '电解液传输性质计算',
      description: '计算电解液体系的离子电导率、扩散系数等传输性质',
      icon: <InteractionOutlined style={{ fontSize: 24, color: '#1C64F2' }} />,
    },
    {
      title: '界面反应与SEI膜形成计算',
      description: '模拟电解液在电极表面的分解反应和SEI膜形成过程',
      icon: <ExperimentOutlined style={{ fontSize: 24, color: '#1C64F2' }} />,
    },
  ];

  return (
    <div>
      <div style={{ marginBottom: 24 }}>
        <Title level={2}>电解液计算</Title>
        <Paragraph style={{ fontSize: 16 }}>
          电解液计算模块提供电池电解液性能的各种计算服务，包括单分子性质、锂离子溶剂化结构、电解液传输性质等计算，
          帮助研究人员设计开发高性能、高安全性的电池电解液体系。
        </Paragraph>
        <Space>
          <Button 
            type="primary" 
            size="large" 
            icon={<ExperimentOutlined />}
            onClick={() => navigate('/battery/electrolyte/workbench')}
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
            title: '高电压电解液设计',
            description: '通过计算不同添加剂的HOMO-LUMO能级和氧化分解电位，筛选能够提高电池电压窗口上限的电解液添加剂。',
            tags: ['高电压电解液', '添加剂筛选', '氧化稳定性'],
          },
          {
            title: '锂盐-溶剂相互作用研究',
            description: '通过量子化学方法计算不同锂盐与碳酸酯溶剂的相互作用强度和配位结构，理解锂盐解离机制和传输特性。',
            tags: ['锂盐-溶剂相互作用', '配位结构', '离子解离'],
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
          onClick={() => navigate('/battery/electrolyte/workbench')}
        >
          开始计算
        </Button>
      </div>
    </div>
  );
};

export default ElectrolyteCalculation; 