import React, { useState, useEffect } from 'react';
import { Typography, Row, Col, Card, Statistic, List, Tag, Button, Progress, Avatar, Spin, message } from 'antd';
import { useNavigate } from 'react-router-dom';
import { 
  ExperimentOutlined, 
  RocketOutlined, 
  HistoryOutlined, 
  ThunderboltOutlined, 
  DatabaseOutlined,
  UserOutlined
} from '@ant-design/icons';
import * as electrolyteService from '../../services/electrolyteService';

const { Title, Paragraph } = Typography;

const Dashboard: React.FC = () => {
  const navigate = useNavigate();
  const [loading, setLoading] = useState(false);
  const [calculations, setCalculations] = useState<any[]>([]);
  
  // 获取计算任务列表
  const fetchCalculations = async () => {
    try {
      setLoading(true);
      const response = await electrolyteService.getCalculations();
      console.log('获取到的计算任务列表:', response.data);
      setCalculations(response.data || []);
    } catch (error) {
      console.error('获取计算任务列表失败:', error);
      message.error('获取计算任务列表失败，请稍后重试');
    } finally {
      setLoading(false);
    }
  };

  // 组件加载时获取计算任务
  useEffect(() => {
    fetchCalculations();
  }, []);
  
  // 获取最近的计算任务
  const getRecentTasks = () => {
    // 按创建时间排序，取最近的3个任务
    return [...calculations]
      .sort((a, b) => new Date(b.created_at).getTime() - new Date(a.created_at).getTime())
      .slice(0, 3)
      .map(calc => ({
        id: calc.id,
        name: calc.name,
        status: calc.status,
        date: new Date(calc.created_at).toLocaleString('zh-CN'),
        module: '电解液计算', // 目前只有电解液计算
        progress: calc.status === 'running' ? 50 : undefined // 模拟进度
      }));
  };

  // 获取统计数据
  const getStatistics = () => {
    return {
      running: calculations.filter(calc => calc.status === 'running').length,
      completed: calculations.filter(calc => calc.status === 'completed').length,
      total: calculations.length,
      // 计算任务总时长（如果有开始和结束时间）
      totalComputeHours: calculations.reduce((total, calc) => {
        if (calc.started_at && calc.finished_at) {
          const startTime = new Date(calc.started_at).getTime();
          const endTime = new Date(calc.finished_at).getTime();
          return total + (endTime - startTime) / (1000 * 60 * 60); // 转换为小时
        }
        return total;
      }, 0)
    };
  };

  const getStatusTag = (status: string) => {
    switch(status) {
      case 'completed':
        return <Tag color="success">已完成</Tag>;
      case 'running':
        return <Tag color="processing">运行中</Tag>;
      case 'pending':
        return <Tag color="default">排队中</Tag>;
      case 'failed':
        return <Tag color="error">失败</Tag>;
      default:
        return <Tag color="default">{status}</Tag>;
    }
  };

  // 获取统计数据
  const stats = getStatistics();

  return (
    <div>
      <div style={{ marginBottom: 24 }}>
        <Row gutter={24} align="middle">
          <Col>
            <Avatar size={64} icon={<UserOutlined />} />
          </Col>
          <Col>
            <Title level={4} style={{ margin: 0 }}>欢迎回来，{localStorage.getItem('username') || '用户'}</Title>
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
              value={stats.running} 
              prefix={<RocketOutlined />} 
              valueStyle={{ color: '#1C64F2' }}
              loading={loading}
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} md={6}>
          <Card>
            <Statistic 
              title="已完成任务" 
              value={stats.completed} 
              prefix={<ExperimentOutlined />} 
              valueStyle={{ color: '#52C41A' }}
              loading={loading}
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} md={6}>
          <Card>
            <Statistic 
              title="计算时长 (小时)" 
              value={stats.totalComputeHours} 
              prefix={<HistoryOutlined />} 
              precision={1}
              valueStyle={{ color: '#1C64F2' }}
              loading={loading}
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} md={6}>
          <Card>
            <Statistic 
              title="存储空间使用" 
              value={stats.total > 0 ? Math.min(stats.total * 5, 95) : 0} // 模拟数据：每个任务占用5%，最高95%
              suffix="%" 
              prefix={<DatabaseOutlined />} 
              valueStyle={{ color: '#1C64F2' }}
              loading={loading}
            />
            <Progress percent={stats.total > 0 ? Math.min(stats.total * 5, 95) : 0} showInfo={false} />
          </Card>
        </Col>
      </Row>
      
      {/* 最近任务 */}
      <Card 
        title="最近计算任务" 
        extra={<Button type="link" onClick={() => navigate('/user/tasks')}>查看全部</Button>}
      >
        {loading ? (
          <div style={{ textAlign: 'center', padding: '20px 0' }}>
            <Spin />
          </div>
        ) : calculations.length === 0 ? (
          <div style={{ textAlign: 'center', padding: '20px 0' }}>
            暂无计算任务，请先创建新的计算任务
          </div>
        ) : (
          <List
            itemLayout="horizontal"
            dataSource={getRecentTasks()}
            renderItem={(item) => (
              <List.Item
                actions={[
                  <Button type="link" onClick={() => navigate(`/calculations/${item.id}`)}>详情</Button>
                ]}
              >
                <List.Item.Meta
                  avatar={<ThunderboltOutlined style={{ color: '#1C64F2', fontSize: 24 }} />}
                  title={<a onClick={() => navigate(`/calculations/${item.id}`)}>{item.name}</a>}
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
        )}
      </Card>
      
      {/* 快速操作 */}
      <Title level={4} style={{ margin: '24px 0 16px' }}>快速操作</Title>
      <Row gutter={[16, 16]}>
        <Col xs={24} sm={12} md={8} lg={6}>
          <Card hoverable onClick={() => navigate('/calculations', { state: { activeTab: 'modules' } })}>
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