import React from 'react';
import { Typography, Row, Col, Card, Statistic, List, Tag, Button, Progress, Avatar } from 'antd';
import { useNavigate } from 'react-router-dom';
import { 
  ExperimentOutlined, 
  RocketOutlined, 
  HistoryOutlined, 
  ThunderboltOutlined, 
  DatabaseOutlined,
  UserOutlined
} from '@ant-design/icons';

const { Title, Paragraph } = Typography;

const Dashboard: React.FC = () => {
  const navigate = useNavigate();

  // 模拟数据
  const recentTasks = [
    {
      id: '1',
      name: '锂离子电池电极材料优化',
      status: 'completed',
      date: '2023-03-01 10:45:23',
      module: '电池-电极材料',
    },
    {
      id: '2',
      name: 'NMC正极材料离子迁移研究',
      status: 'running',
      date: '2023-03-05 14:22:10',
      module: '电池-电极材料',
      progress: 68,
    },
    {
      id: '3',
      name: 'CO2还原反应机理分析',
      status: 'queued',
      date: '2023-03-06 11:10:30',
      module: '催化-电催化',
    },
  ];
  
  const getStatusTag = (status: string) => {
    switch(status) {
      case 'completed':
        return <Tag color="success">已完成</Tag>;
      case 'running':
        return <Tag color="processing">运行中</Tag>;
      case 'queued':
        return <Tag color="default">排队中</Tag>;
      case 'failed':
        return <Tag color="error">失败</Tag>;
      default:
        return <Tag color="default">{status}</Tag>;
    }
  };

  return (
    <div>
      <div style={{ marginBottom: 24 }}>
        <Row gutter={24} align="middle">
          <Col>
            <Avatar size={64} icon={<UserOutlined />} />
          </Col>
          <Col>
            <Title level={4} style={{ margin: 0 }}>欢迎回来，用户</Title>
            <Paragraph style={{ marginBottom: 0 }}>
              今天是 {new Date().toLocaleDateString('zh-CN', { weekday: 'long', year: 'numeric', month: 'long', day: 'numeric' })}
            </Paragraph>
          </Col>
        </Row>
      </div>
      
      {/* 统计卡片 */}
      <Row gutter={[16, 16]} style={{ marginBottom: 24 }}>
        <Col xs={24} sm={12} md={6}>
          <Card>
            <Statistic 
              title="运行中任务" 
              value={1} 
              prefix={<RocketOutlined />} 
              valueStyle={{ color: '#1C64F2' }}
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} md={6}>
          <Card>
            <Statistic 
              title="已完成任务" 
              value={8} 
              prefix={<ExperimentOutlined />} 
              valueStyle={{ color: '#52C41A' }}
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} md={6}>
          <Card>
            <Statistic 
              title="计算时长 (小时)" 
              value={27.5} 
              prefix={<HistoryOutlined />} 
              precision={1}
              valueStyle={{ color: '#1C64F2' }}
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} md={6}>
          <Card>
            <Statistic 
              title="存储空间使用" 
              value={34} 
              suffix="%" 
              prefix={<DatabaseOutlined />} 
              valueStyle={{ color: '#1C64F2' }}
            />
            <Progress percent={34} showInfo={false} />
          </Card>
        </Col>
      </Row>
      
      {/* 最近任务 */}
      <Card title="最近计算任务" extra={<Button type="link" onClick={() => navigate('/user/tasks')}>查看全部</Button>}>
        <List
          itemLayout="horizontal"
          dataSource={recentTasks}
          renderItem={(item) => (
            <List.Item
              actions={[
                <Button type="link" onClick={() => navigate('/user/tasks')}>详情</Button>
              ]}
            >
              <List.Item.Meta
                avatar={<ThunderboltOutlined style={{ color: '#1C64F2', fontSize: 24 }} />}
                title={<a onClick={() => navigate('/user/tasks')}>{item.name}</a>}
                description={
                  <span>
                    {item.module} | {item.date} | {getStatusTag(item.status)}
                    {item.status === 'running' && item.progress && (
                      <Progress percent={item.progress} size="small" style={{ maxWidth: 100, display: 'inline-block', marginLeft: 8 }} />
                    )}
                  </span>
                }
              />
            </List.Item>
          )}
        />
      </Card>
      
      {/* 快速操作 */}
      <Title level={4} style={{ margin: '24px 0 16px' }}>快速操作</Title>
      <Row gutter={[16, 16]}>
        <Col xs={24} sm={12} md={8} lg={6}>
          <Card hoverable onClick={() => navigate('/calculations')}>
            <Statistic 
              title="新建计算任务" 
              value={" "} 
              prefix={<RocketOutlined style={{ fontSize: 36, color: '#1C64F2' }} />} 
              valueStyle={{ fontSize: 0 }}
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} md={8} lg={6}>
          <Card hoverable onClick={() => navigate('/user/tasks')}>
            <Statistic 
              title="查看计算任务" 
              value={" "} 
              prefix={<ExperimentOutlined style={{ fontSize: 36, color: '#1C64F2' }} />} 
              valueStyle={{ fontSize: 0 }}
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} md={8} lg={6}>
          <Card hoverable onClick={() => navigate('/user/data')}>
            <Statistic 
              title="数据管理" 
              value={" "} 
              prefix={<DatabaseOutlined style={{ fontSize: 36, color: '#1C64F2' }} />} 
              valueStyle={{ fontSize: 0 }}
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} md={8} lg={6}>
          <Card hoverable onClick={() => navigate('/documentation')}>
            <Statistic 
              title="查看文档" 
              value={" "} 
              prefix={<HistoryOutlined style={{ fontSize: 36, color: '#1C64F2' }} />} 
              valueStyle={{ fontSize: 0 }}
            />
          </Card>
        </Col>
      </Row>
    </div>
  );
};

export default Dashboard; 