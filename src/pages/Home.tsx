import React from 'react';
import { Typography, Row, Col, Card, Button, Space, Carousel, Divider, List, Tag } from 'antd';
import { useNavigate } from 'react-router-dom';
import { 
  ExperimentOutlined, 
  RocketOutlined, 
  FileTextOutlined, 
  TeamOutlined,
  BulbOutlined,
  ThunderboltOutlined,
  BarChartOutlined
} from '@ant-design/icons';
import { globalStyles } from '@/theme';

const { Title, Paragraph } = Typography;

const Home: React.FC = () => {
  const navigate = useNavigate();
  
  // 功能模块数据
  const moduleItems = [
    {
      title: '电池材料计算',
      description: '电极材料、电解液和全电池系统的性能预测与优化',
      icon: <ThunderboltOutlined style={{ fontSize: 36, color: '#1C64F2' }} />,
      route: '/battery',
    },
    {
      title: '催化材料计算',
      description: '表面催化、电催化和光催化反应机理与性能评估',
      icon: <ExperimentOutlined style={{ fontSize: 36, color: '#1C64F2' }} />,
      route: '/catalysis',
    },
    {
      title: '文档中心',
      description: '详细的使用教程、计算方法说明和科学参考资料',
      icon: <FileTextOutlined style={{ fontSize: 36, color: '#1C64F2' }} />,
      route: '/documentation',
    },
  ];
  
  // 研究案例数据
  const caseStudies = [
    {
      title: '高性能锂电池正极材料设计',
      description: '通过第一性原理计算和机器学习方法，从5000个候选材料中筛选出3个高容量、长循环寿命的正极材料。',
      image: 'https://picsum.photos/600/300',
      tags: ['电池材料', '第一性原理', '机器学习'],
    },
    {
      title: 'CO₂电催化还原反应机理研究',
      description: '揭示了Cu-Sn合金催化剂上CO₂还原为甲醇的反应路径，为设计高效催化剂提供理论指导。',
      image: 'https://picsum.photos/600/300',
      tags: ['电催化', '反应机理', 'CO₂还原'],
    },
    {
      title: '纳米多孔硅的热输运性能预测',
      description: '采用分子动力学模拟方法，系统研究了纳米多孔硅的热导率随孔隙率、孔径和表面化学修饰的变化关系。',
      image: 'https://picsum.photos/600/300',
      tags: ['热输运', '分子动力学', '纳米材料'],
    },
  ];
  
  // 平台功能数据
  const features = [
    {
      title: '多尺度计算',
      description: '从原子尺度到连续介质的多层次材料模拟',
      icon: <BulbOutlined />,
    },
    {
      title: '高通量计算',
      description: '自动化的材料筛选与优化工作流',
      icon: <RocketOutlined />,
    },
    {
      title: '交互式可视化',
      description: '先进的3D结构与属性可视化系统',
      icon: <BarChartOutlined />,
    },
    {
      title: '协作与共享',
      description: '团队协作与计算结果共享功能',
      icon: <TeamOutlined />,
    },
  ];
  
  return (
    <div>
      {/* 头部Banner */}
      <div style={{ 
        padding: '60px 0', 
        background: 'linear-gradient(120deg, #1c64f2 0%, #3c83f6 100%)',
        borderRadius: 8,
        marginBottom: 40,
        color: 'white',
        textAlign: 'center'
      }}>
        <Title level={1} style={{ color: 'white', marginBottom: 24 }}>
          计算材料科学平台
        </Title>
        <Paragraph style={{ fontSize: 18, color: 'white', maxWidth: 800, margin: '0 auto 32px' }}>
          面向材料科学研究的在线计算服务平台，提供电池材料、催化材料等多种计算服务，
          简化复杂的材料计算流程，帮助研究人员高效进行材料设计与性能预测。
        </Paragraph>
        <Space size="large">
          <Button type="primary" size="large" ghost onClick={() => navigate('/battery')}>
            开始计算
          </Button>
          <Button size="large" style={{ background: 'white', color: '#1C64F2' }} onClick={() => navigate('/documentation')}>
            了解更多
          </Button>
        </Space>
      </div>
      
      {/* 功能模块区域 */}
      <div style={{ marginBottom: 60 }}>
        <Title level={2} style={{ textAlign: 'center', marginBottom: 40 }}>
          主要功能模块
        </Title>
        <Row gutter={[24, 24]}>
          {moduleItems.map((item, index) => (
            <Col xs={24} sm={24} md={8} key={index}>
              <Card 
                hoverable
                style={{ height: '100%', textAlign: 'center', borderRadius: 8, boxShadow: globalStyles.cardShadow }}
                onClick={() => navigate(item.route)}
              >
                <div style={{ margin: '16px 0' }}>
                  {item.icon}
                </div>
                <Title level={4}>{item.title}</Title>
                <Paragraph>{item.description}</Paragraph>
                <Button type="primary">进入模块</Button>
              </Card>
            </Col>
          ))}
        </Row>
      </div>
      
      {/* 研究案例 */}
      <div style={{ marginBottom: 60 }}>
        <Title level={2} style={{ textAlign: 'center', marginBottom: 40 }}>
          研究案例展示
        </Title>
        <Carousel autoplay dots={{ className: 'custom-carousel-dots' }}>
          {caseStudies.map((study, index) => (
            <div key={index}>
              <Row gutter={24} align="middle">
                <Col xs={24} sm={24} md={12}>
                  <img 
                    src={study.image} 
                    alt={study.title} 
                    style={{ width: '100%', borderRadius: 8, boxShadow: globalStyles.cardShadow }}
                  />
                </Col>
                <Col xs={24} sm={24} md={12}>
                  <Card style={{ borderRadius: 8, boxShadow: 'none' }}>
                    <Title level={3}>{study.title}</Title>
                    <Space style={{ marginBottom: 16 }}>
                      {study.tags.map((tag, idx) => (
                        <Tag color="blue" key={idx}>{tag}</Tag>
                      ))}
                    </Space>
                    <Paragraph style={{ fontSize: 16 }}>{study.description}</Paragraph>
                    <Button type="link" style={{ paddingLeft: 0 }}>查看详情</Button>
                  </Card>
                </Col>
              </Row>
            </div>
          ))}
        </Carousel>
      </div>
      
      {/* 平台特点 */}
      <div style={{ marginBottom: 60, background: '#f5f8ff', padding: '40px 20px', borderRadius: 8 }}>
        <Title level={2} style={{ textAlign: 'center', marginBottom: 40 }}>
          平台特点
        </Title>
        <Row gutter={[32, 32]}>
          {features.map((feature, index) => (
            <Col xs={24} sm={12} md={6} key={index}>
              <Card 
                style={{ height: '100%', textAlign: 'center', borderRadius: 8 }}
                bordered={false}
              >
                <div style={{ fontSize: 36, color: '#1C64F2', marginBottom: 16 }}>
                  {feature.icon}
                </div>
                <Title level={4}>{feature.title}</Title>
                <Paragraph>{feature.description}</Paragraph>
              </Card>
            </Col>
          ))}
        </Row>
      </div>
      
      {/* 最新动态 */}
      <div style={{ marginBottom: 60 }}>
        <Title level={2} style={{ textAlign: 'center', marginBottom: 20 }}>
          最新动态
        </Title>
        <Divider />
        <List
          itemLayout="horizontal"
          dataSource={[
            {
              title: '平台更新 | 新增量子化学计算模块',
              date: '2023-10-15',
            },
            {
              title: '学术动态 | 基于人工智能的电池材料设计研讨会',
              date: '2023-10-10',
            },
            {
              title: '用户公告 | 计算资源扩容完成，任务处理能力提升50%',
              date: '2023-10-05',
            },
          ]}
          renderItem={(item) => (
            <List.Item>
              <List.Item.Meta
                title={<a href="#" style={{ fontSize: 16 }}>{item.title}</a>}
                description={item.date}
              />
            </List.Item>
          )}
        />
        <div style={{ textAlign: 'center', marginTop: 20 }}>
          <Button>查看更多动态</Button>
        </div>
      </div>
      
      {/* 用户行动区 */}
      <div style={{ 
        background: 'linear-gradient(120deg, #1c64f2 0%, #3c83f6 100%)',
        padding: '40px 20px',
        borderRadius: 8,
        textAlign: 'center',
        color: 'white'
      }}>
        <Title level={2} style={{ color: 'white', marginBottom: 16 }}>
          立即开始您的材料计算之旅
        </Title>
        <Paragraph style={{ fontSize: 16, color: 'white', maxWidth: 600, margin: '0 auto 24px' }}>
          注册账号，享受功能完备的计算材料科学平台服务，加速您的研究进程
        </Paragraph>
        <Space size="large">
          <Button type="primary" size="large" style={{ background: 'white', color: '#1C64F2' }} onClick={() => navigate('/register')}>
            免费注册
          </Button>
          <Button ghost size="large" onClick={() => navigate('/documentation')}>
            查看文档
          </Button>
        </Space>
      </div>
    </div>
  );
};

export default Home; 