import React from 'react';
import { Typography, Row, Col, Card, Avatar, Divider, Space } from 'antd';
import { 
  ExperimentOutlined, 
  SafetyOutlined, 
  TeamOutlined, 
  RocketOutlined, 
  MailOutlined, 
  GithubOutlined, 
  LinkedinOutlined 
} from '@ant-design/icons';

const { Title, Paragraph, Text } = Typography;

const About: React.FC = () => {
  const features = [
    {
      icon: <ExperimentOutlined style={{ fontSize: 32, color: '#1C64F2' }} />,
      title: '高性能计算',
      description: '利用先进的并行计算技术和优化的算法，提供高效、准确的材料科学计算能力。'
    },
    {
      icon: <SafetyOutlined style={{ fontSize: 32, color: '#1C64F2' }} />,
      title: '数据安全',
      description: '采用多重安全机制保护您的研究数据和计算结果，确保数据的私密性和完整性。'
    },
    {
      icon: <TeamOutlined style={{ fontSize: 32, color: '#1C64F2' }} />,
      title: '专家支持',
      description: '由材料科学和计算机专家团队提供技术支持，解决您在研究中遇到的问题。'
    },
    {
      icon: <RocketOutlined style={{ fontSize: 32, color: '#1C64F2' }} />,
      title: '持续创新',
      description: '不断更新和完善计算方法和模型，跟踪材料科学领域的最新研究进展。'
    }
  ];

  const team = [
    {
      name: '张教授',
      title: '首席科学家',
      avatar: 'https://xsgames.co/randomusers/avatar.php?g=male&seed=1',
      bio: '材料科学博士，专注于新能源材料研究20年，发表学术论文100余篇。',
      social: [
        { icon: <MailOutlined />, link: 'mailto:zhang@example.com' },
        { icon: <GithubOutlined />, link: 'https://github.com/' },
        { icon: <LinkedinOutlined />, link: 'https://linkedin.com/' }
      ]
    },
    {
      name: '李博士',
      title: '计算方法专家',
      avatar: 'https://xsgames.co/randomusers/avatar.php?g=female&seed=2',
      bio: '计算物理学博士，开发多种高效计算算法，在计算材料学领域有丰富经验。',
      social: [
        { icon: <MailOutlined />, link: 'mailto:li@example.com' },
        { icon: <GithubOutlined />, link: 'https://github.com/' },
        { icon: <LinkedinOutlined />, link: 'https://linkedin.com/' }
      ]
    },
    {
      name: '王工程师',
      title: '软件开发负责人',
      avatar: 'https://xsgames.co/randomusers/avatar.php?g=male&seed=3',
      bio: '10年科学计算软件开发经验，专注于高性能计算和用户友好界面设计。',
      social: [
        { icon: <MailOutlined />, link: 'mailto:wang@example.com' },
        { icon: <GithubOutlined />, link: 'https://github.com/' },
        { icon: <LinkedinOutlined />, link: 'https://linkedin.com/' }
      ]
    }
  ];

  return (
    <div style={{ padding: '40px 24px' }}>
      <div>
        {/* 平台介绍 */}
        <div style={{ textAlign: 'center', marginBottom: 60 }}>
          <Title level={1}>关于我们</Title>
          <Paragraph style={{ fontSize: 16, maxWidth: 800, margin: '24px auto' }}>
            材料科学计算平台是一个致力于为材料科学研究人员提供高效计算工具的在线服务平台。
            我们的使命是通过先进的计算方法和用户友好的界面，简化复杂的材料科学计算过程，
            帮助研究人员更快地设计、分析和优化新材料。
          </Paragraph>
        </div>

        {/* 特色功能 */}
        <div style={{ marginBottom: 60 }}>
          <Title level={2} style={{ textAlign: 'center', marginBottom: 40 }}>平台特色</Title>
          <Row gutter={[32, 32]}>
            {features.map((feature, index) => (
              <Col xs={24} sm={12} md={6} lg={6} key={index}>
                <Card style={{ height: 240, display: 'flex', flexDirection: 'column', alignItems: 'center' }}>
                  <div style={{ marginBottom: 16 }}>{feature.icon}</div>
                  <Title level={4} style={{ marginTop: 0 }}>{feature.title}</Title>
                  <Paragraph style={{ textAlign: 'center' }}>{feature.description}</Paragraph>
                </Card>
              </Col>
            ))}
          </Row>
        </div>

        <Divider />

        {/* 团队介绍 */}
        <div style={{ marginTop: 60 }}>
          <Title level={2} style={{ textAlign: 'center', marginBottom: 40 }}>核心团队</Title>
          <Row gutter={[32, 32]}>
            {team.map((member, index) => (
              <Col xs={24} sm={12} md={8} lg={8} key={index}>
                <Card style={{ textAlign: 'center' }}>
                  <Avatar size={100} src={member.avatar} />
                  <Title level={4} style={{ marginTop: 16, marginBottom: 4 }}>{member.name}</Title>
                  <Text type="secondary" style={{ display: 'block', marginBottom: 16 }}>{member.title}</Text>
                  <Paragraph>{member.bio}</Paragraph>
                  <Space>
                    {member.social.map((social, idx) => (
                      <a href={social.link} target="_blank" rel="noopener noreferrer" key={idx} style={{ color: '#1C64F2' }}>
                        {social.icon}
                      </a>
                    ))}
                  </Space>
                </Card>
              </Col>
            ))}
          </Row>
        </div>

        {/* 联系信息 */}
        <div style={{ marginTop: 60, textAlign: 'center' }}>
          <Title level={2} style={{ marginBottom: 24 }}>联系我们</Title>
          <Paragraph>
            邮箱: <a href="mailto:contact@materialsci-platform.com">contact@materialsci-platform.com</a>
          </Paragraph>
          <Paragraph>
            地址: 北京市海淀区中关村科技园区
          </Paragraph>
        </div>
      </div>
    </div>
  );
};

export default About; 