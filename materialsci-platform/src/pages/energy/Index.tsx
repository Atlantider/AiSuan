import React from 'react';
import { Typography, Card, Button, Alert, Space } from 'antd';
import { useNavigate } from 'react-router-dom';
import { LineChartOutlined, ArrowLeftOutlined } from '@ant-design/icons';

const { Title, Paragraph } = Typography;

const EnergyIndex: React.FC = () => {
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
          <LineChartOutlined style={{ fontSize: 64, color: '#1C64F2' }} />
          <Title level={2} style={{ marginTop: 16 }}>能源材料计算</Title>
        </div>

        <Paragraph style={{ fontSize: 16 }}>
          能源材料计算模块致力于新能源材料的设计与优化，包括太阳能电池材料、热电材料、
          储能材料等。通过第一性原理和多尺度模拟方法，预测材料性能并指导实验合成。
        </Paragraph>

        <Alert
          message="模块开发中"
          description="能源材料计算模块正在积极开发中，即将推出。敬请期待！"
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

export default EnergyIndex; 