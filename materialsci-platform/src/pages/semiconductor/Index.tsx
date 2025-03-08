import React from 'react';
import { Typography, Card, Button, Alert, Space } from 'antd';
import { useNavigate } from 'react-router-dom';
import { DotChartOutlined, ArrowLeftOutlined } from '@ant-design/icons';

const { Title, Paragraph } = Typography;

const SemiconductorIndex: React.FC = () => {
  const navigate = useNavigate();

  return (
    <div>
      <div style={{ marginBottom: 24 }}>
        <Button 
          icon={<ArrowLeftOutlined />} 
          onClick={() => navigate('/calculations')}
        >
          返回计算类型选择
        </Button>
      </div>

      <Card>
        <div style={{ textAlign: 'center', marginBottom: 24 }}>
          <DotChartOutlined style={{ fontSize: 64, color: '#1C64F2' }} />
          <Title level={2} style={{ marginTop: 16 }}>半导体材料计算</Title>
        </div>

        <Paragraph style={{ fontSize: 16 }}>
          半导体材料计算模块提供全面的计算工具，帮助您设计和优化半导体材料性能。
          您可以计算电子结构、带隙、迁移率、掺杂效应等关键参数，为开发新型半导体器件提供理论依据。
        </Paragraph>

        <Alert
          message="模块开发中"
          description="半导体材料计算模块正在积极开发中，即将推出。敬请期待！"
          type="info"
          showIcon
          style={{ marginBottom: 24 }}
        />

        <div style={{ textAlign: 'center' }}>
          <Space>
            <Button onClick={() => navigate('/calculations')}>返回选择</Button>
            <Button onClick={() => navigate('/documentation')}>查看文档</Button>
          </Space>
        </div>
      </Card>
    </div>
  );
};

export default SemiconductorIndex; 