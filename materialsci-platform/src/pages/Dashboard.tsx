import React from 'react';
import { Typography, Row, Col, Card, Statistic, List, Tag, Button, Progress, Avatar } from 'antd';
import {
  RocketOutlined,
  ExperimentOutlined,
  DatabaseOutlined,
  HistoryOutlined,
  CheckCircleOutlined,
  ClockCircleOutlined,
  SyncOutlined,
  StopOutlined,
  UserOutlined
} from '@ant-design/icons';
import { useNavigate } from 'react-router-dom';

const { Title } = Typography;

const Dashboard: React.FC = () => {
  const navigate = useNavigate();
  
  // 模拟用户数据
  const userInfo = {
    name: '张三',
    avatar: 'https://xsgames.co/randomusers/avatar.php?g=male',
    runningTasks: 3,
    completedTasks: 12,
    computationTime: '124小时',
    storageUsage: 75 // 百分比
  };
  
  // 模拟最近任务数据
  const recentTasks = [
    {
      id: 'TASK-2023-001',
      name: 'LiCoO2电池材料计算',
      status: 'running',
      date: '2023-06-05 10:30',
      module: 'battery'
    },
    {
      id: 'TASK-2023-002',
      name: 'Si-C复合负极材料',
      status: 'completed',
      date: '2023-06-03 15:45',
      module: 'battery'
    },
    {
      id: 'TASK-2023-003',
      name: 'LiFePO4充放电分析',
      status: 'queued',
      date: '2023-06-05 08:15',
      module: 'battery'
    },
    {
      id: 'TASK-2023-004',
      name: 'Pt催化剂表面吸附',
      status: 'failed',
      date: '2023-06-02 14:20',
      module: 'catalyst'
    }
  ];
  
  // 获取状态标签
  const getStatusTag = (status: string) => {
    switch (status) {
      case 'completed':
        return <Tag icon={<CheckCircleOutlined />} color="success">已完成</Tag>;
      case 'running':
        return <Tag icon={<SyncOutlined spin />} color="processing">运行中</Tag>;
      case 'queued':
        return <Tag icon={<ClockCircleOutlined />} color="default">排队中</Tag>;
      case 'failed':
        return <Tag icon={<StopOutlined />} color="error">失败</Tag>;
      default:
        return <Tag>未知</Tag>;
    }
  };

  // 跳转到任务页面
  const goToTasks = () => {
    navigate('/user/tasks');
  };

  return (
    <div>
      {/* 欢迎信息 */}
      <Row gutter={[24, 24]} style={{ marginBottom: 24 }}>
        <Col span={24}>
          <Card>
            <div style={{ display: 'flex', alignItems: 'center' }}>
              <Avatar size={64} src={userInfo.avatar} icon={<UserOutlined />} />
              <div style={{ marginLeft: 24 }}>
                <Title level={4} style={{ margin: 0 }}>欢迎回来，{userInfo.name}</Title>
                <p style={{ color: '#8c8c8c', margin: '8px 0 0 0' }}>
                  今天是 {new Date().toLocaleDateString('zh-CN', { year: 'numeric', month: 'long', day: 'numeric' })}
                </p>
              </div>
            </div>
          </Card>
        </Col>
      </Row>

      {/* 统计数据 */}
      <Row gutter={[24, 24]} style={{ marginBottom: 24 }}>
        <Col xs={24} sm={12} md={6}>
          <Card>
            <Statistic 
              title="运行中任务" 
              value={userInfo.runningTasks} 
              valueStyle={{ color: '#1890ff' }} 
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} md={6}>
          <Card>
            <Statistic 
              title="已完成任务" 
              value={userInfo.completedTasks} 
              valueStyle={{ color: '#52c41a' }} 
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} md={6}>
          <Card>
            <Statistic 
              title="计算总时长" 
              value={userInfo.computationTime} 
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} md={6}>
          <Card>
            <Statistic 
              title="存储使用量" 
              value={userInfo.storageUsage} 
              suffix="%" 
            />
            <Progress percent={userInfo.storageUsage} status={userInfo.storageUsage > 80 ? "exception" : "normal"} />
          </Card>
        </Col>
      </Row>

      {/* 最近任务 */}
      <Row gutter={[24, 24]} style={{ marginBottom: 24 }}>
        <Col span={24}>
          <Card 
            title="最近任务" 
            extra={<Button type="link" onClick={goToTasks}>查看全部</Button>}
          >
            <List
              dataSource={recentTasks}
              renderItem={item => (
                <List.Item 
                  key={item.id}
                  actions={[
                    <Button type="link" onClick={() => console.log('查看详情', item.id)}>查看详情</Button>
                  ]}
                >
                  <List.Item.Meta
                    title={item.name}
                    description={`${item.date} | 任务ID: ${item.id}`}
                  />
                  {getStatusTag(item.status)}
                </List.Item>
              )}
            />
          </Card>
        </Col>
      </Row>

      {/* 快速操作 */}
      <Row gutter={[24, 24]}>
        <Col xs={24} sm={12} md={6}>
          <Card hoverable onClick={() => navigate('/battery')}>
            <Statistic 
              title="新建计算任务" 
              value={" "} 
              prefix={<RocketOutlined style={{ fontSize: 36, color: '#1C64F2' }} />} 
              valueStyle={{ fontSize: 0 }}
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} md={6}>
          <Card hoverable onClick={goToTasks}>
            <Statistic 
              title="查看计算任务" 
              value={" "} 
              prefix={<ExperimentOutlined style={{ fontSize: 36, color: '#1C64F2' }} />} 
              valueStyle={{ fontSize: 0 }}
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} md={6}>
          <Card hoverable onClick={() => navigate('/user/data')}>
            <Statistic 
              title="数据管理" 
              value={" "} 
              prefix={<DatabaseOutlined style={{ fontSize: 36, color: '#1C64F2' }} />} 
              valueStyle={{ fontSize: 0 }}
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} md={6}>
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