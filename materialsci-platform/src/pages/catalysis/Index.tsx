import React from 'react';
import { Typography, Button, Card, Row, Col, Space } from 'antd';
import { useNavigate } from 'react-router-dom';

const { Title, Paragraph } = Typography;

const CatalysisIndex: React.FC = () => {
  const navigate = useNavigate();

  // 催化计算子模块
  const catalysisModules = [
    {
      title: '表面催化计算',
      description: '计算各类催化剂表面的吸附能、反应路径与能垒等参数',
      path: '/catalysis/surface',
    },
    {
      title: '电催化计算',
      description: '计算电催化反应机理和性能，包括HER、ORR、CO₂RR等反应',
      path: '/catalysis/electro',
    },
    {
      title: '光催化计算',
      description: '光催化材料设计和光生载流子动力学计算及反应机理分析',
      path: '/catalysis/photo',
    },
  ];

  return (
    <div>
      <div style={{ marginBottom: 24 }}>
        <Title level={2}>催化材料计算</Title>
        <Paragraph style={{ fontSize: 16 }}>
          催化材料计算模块提供表面催化、电催化和光催化反应机理与性能评估功能，帮助研究人员深入理解催化过程中的电子结构变化和反应机理，设计开发高效催化材料。
        </Paragraph>
      </div>

      <Row gutter={[24, 24]}>
        {catalysisModules.map((module, index) => (
          <Col xs={24} md={8} key={index}>
            <Card 
              hoverable
              title={<span style={{ fontSize: 18 }}>{module.title}</span>}
              style={{ height: '100%' }}
              onClick={() => navigate(module.path)}
            >
              <Paragraph>{module.description}</Paragraph>
              <Button type="primary">开始计算</Button>
            </Card>
          </Col>
        ))}
      </Row>

      <div style={{ marginTop: 48 }}>
        <Title level={3}>常见应用场景</Title>
        <Row gutter={[24, 24]}>
          <Col xs={24} md={8}>
            <Card title="CO₂电催化还原优化">
              <Paragraph>设计CO₂转化为高附加值产品（如甲醇、乙醇）的高效催化剂。</Paragraph>
            </Card>
          </Col>
          <Col xs={24} md={8}>
            <Card title="析氢反应(HER)催化剂开发">
              <Paragraph>寻找高效、低成本的氢能源催化材料，降低氢生产成本。</Paragraph>
            </Card>
          </Col>
          <Col xs={24} md={8}>
            <Card title="光催化水分解材料设计">
              <Paragraph>设计利用太阳能实现水分解制氢的高效光催化材料。</Paragraph>
            </Card>
          </Col>
        </Row>
      </div>

      <div style={{ marginTop: 48, textAlign: 'center' }}>
        <Title level={4}>需要帮助？</Title>
        <Space>
          <Button onClick={() => navigate('/documentation')}>查看文档</Button>
          <Button type="primary">联系我们</Button>
        </Space>
      </div>
    </div>
  );
};

export default CatalysisIndex; 