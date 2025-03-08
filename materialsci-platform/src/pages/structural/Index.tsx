import React from 'react';
import { Typography, Card, Button, Alert, Space } from 'antd';
import { useNavigate } from 'react-router-dom';
import { ApartmentOutlined, ArrowLeftOutlined } from '@ant-design/icons';

const { Title, Paragraph } = Typography;

const StructuralIndex: React.FC = () => {
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
          <ApartmentOutlined style={{ fontSize: 64, color: '#1C64F2' }} />
          <Title level={2} style={{ marginTop: 16 }}>结构材料计算</Title>
        </div>

        <Paragraph style={{ fontSize: 16 }}>
          结构材料计算模块专注于研究材料的力学性能、热学性能和微观结构特性。
          通过分子动力学、有限元分析和多尺度模拟等方法，预测材料的强度、韧性、热稳定性等关键参数，
          为高性能结构材料的设计和优化提供理论依据。
        </Paragraph>

        <Alert
          message="模块开发中"
          description="结构材料计算模块正在积极开发中，即将推出。敬请期待！"
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

export default StructuralIndex; 