import React from 'react';
import { Typography, Button, Space, Card, Row, Col } from 'antd';
import { ExperimentOutlined, DatabaseOutlined, ReadOutlined, TeamOutlined } from '@ant-design/icons';
import { useNavigate } from 'react-router-dom';

const { Title, Paragraph } = Typography;

const Home: React.FC = () => {
  const navigate = useNavigate();

  const features = [
    {
      title: '计算工具',
      icon: <ExperimentOutlined style={{ fontSize: 36, color: '#1C64F2' }} />,
      description: '提供电池、半导体、催化剂等多种材料的计算功能，快速得到可靠结果。',
      action: () => navigate('/calculations')
    },
    {
      title: '数据管理',
      icon: <DatabaseOutlined style={{ fontSize: 36, color: '#1C64F2' }} />,
      description: '集中管理您的计算结果、实验数据和分析模型，轻松查找和分享。',
      action: () => navigate('/user/data')
    },
    {
      title: '学习资源',
      icon: <ReadOutlined style={{ fontSize: 36, color: '#1C64F2' }} />,
      description: '获取最新的材料科学研究方法、教程和计算技巧，提升研究效率。',
      action: () => navigate('/documentation')
    },
    {
      title: '协作空间',
      icon: <TeamOutlined style={{ fontSize: 36, color: '#1C64F2' }} />,
      description: '与团队成员共享项目和数据，实现高效协作和研究成果共享。',
      action: () => navigate('/collaboration')
    }
  ];

  return (
    <div style={{ padding: '40px 24px' }}>
      <div>
        {/* 头部区域 */}
        <div style={{ textAlign: 'center', marginBottom: 60 }}>
          <Title level={1} style={{ fontSize: 42 }}>
            材料科学计算平台
          </Title>
          <Paragraph style={{ fontSize: 18, color: '#666', margin: '24px auto', maxWidth: '800px' }}>
            一站式解决方案，为材料科学研究提供高效、精准的计算工具和数据分析能力
          </Paragraph>
          <Space size="large">
            <Button type="primary" size="large" onClick={() => navigate('/user/dashboard')}>
              开始使用
            </Button>
            <Button size="large" onClick={() => navigate('/documentation')}>
              了解更多
            </Button>
          </Space>
        </div>

        {/* 特性区域 */}
        <div style={{ margin: '60px 0' }}>
          <Title level={2} style={{ textAlign: 'center', marginBottom: 40 }}>
            平台特性
          </Title>
          <Row gutter={[24, 24]}>
            {features.map((feature, index) => (
              <Col xs={24} sm={12} md={6} lg={6} key={index}>
                <Card
                  hoverable
                  style={{ height: 320 }}
                  onClick={feature.action}
                >
                  <div style={{ textAlign: 'center', marginBottom: 24 }}>
                    {feature.icon}
                  </div>
                  <Title level={4} style={{ textAlign: 'center' }}>
                    {feature.title}
                  </Title>
                  <Paragraph style={{ textAlign: 'center', marginTop: 12 }}>
                    {feature.description}
                  </Paragraph>
                </Card>
              </Col>
            ))}
          </Row>
        </div>
      </div>
    </div>
  );
};

export default Home; 