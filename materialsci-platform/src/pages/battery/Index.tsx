import React from 'react';
import { Typography, Button, Card, Row, Col, Space } from 'antd';
import { useNavigate } from 'react-router-dom';

const { Title, Paragraph } = Typography;

const BatteryIndex: React.FC = () => {
  const navigate = useNavigate();

  // 电池计算子模块
  const batteryModules = [
    {
      title: '电极材料计算',
      description: '计算电极材料的电极电位、离子迁移能垒、形成能等性质',
      path: '/battery/electrode-material',
    },
    {
      title: '电解液计算',
      description: '计算电解液的单分子性质、锂离子溶剂化结构和传输特性',
      path: '/battery/electrolyte',
    },
    {
      title: '全电池系统计算',
      description: '整合电极和电解液参数，进行全电池系统性能预测',
      path: '/battery/full-battery',
    },
  ];

  return (
    <div>
      <div style={{ marginBottom: 24 }}>
        <Title level={2}>电池材料计算</Title>
        <Paragraph style={{ fontSize: 16 }}>
          电池材料计算模块提供电极材料、电解液和全电池系统的性能预测与优化功能，帮助研究人员深入了解电池工作机理，设计开发高性能电池材料。
        </Paragraph>
      </div>

      <Row gutter={[24, 24]}>
        {batteryModules.map((module, index) => (
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
            <Card title="锂离子电池正极材料优化">
              <Paragraph>提高正极材料的容量和循环稳定性，减少Li迁移能垒，优化电极电位。</Paragraph>
            </Card>
          </Col>
          <Col xs={24} md={8}>
            <Card title="固态电解质设计">
              <Paragraph>预测新型固态电解质材料的离子电导率、电化学窗口和界面稳定性。</Paragraph>
            </Card>
          </Col>
          <Col xs={24} md={8}>
            <Card title="新型电池体系探索">
              <Paragraph>Na离子、K离子、Mg离子等新型电池体系的材料设计与性能评估。</Paragraph>
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

export default BatteryIndex; 