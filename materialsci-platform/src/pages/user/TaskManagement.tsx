import React, { useState, useEffect } from 'react';
import { useNavigate } from 'react-router-dom';
import { Card, Typography, Button, Spin, Table, Tag, Space } from 'antd';
import * as electrolyteService from '../../services/electrolyteService';

const { Title, Paragraph } = Typography;

const TaskManagement: React.FC = () => {
  const navigate = useNavigate();
  const [loading, setLoading] = useState(false);
  const [calculations, setCalculations] = useState<any[]>([]);

  useEffect(() => {
    fetchCalculations();
  }, []);

  const fetchCalculations = async () => {
    setLoading(true);
    try {
      const response = await electrolyteService.getCalculations();
      setCalculations(response.data || []);
    } catch (error) {
      console.error('获取计算任务列表失败:', error);
    } finally {
      setLoading(false);
    }
  };

  // 处理任务点击事件
  const handleTaskClick = (taskId: number) => {
    // 导航到任务详情页面
    navigate(`/user/tasks/${taskId}`);
  };

  const columns = [
    {
      title: '任务名称',
      dataIndex: 'name',
      key: 'name',
      render: (text: string, record: any) => (
        <a onClick={() => handleTaskClick(record.id)}>{text}</a>
      ),
    },
    {
      title: '状态',
      dataIndex: 'status',
      key: 'status',
      render: (status: string) => {
        let color = 'default';
        let text = status;
        switch (status) {
          case 'completed':
            color = 'success';
            text = '已完成';
            break;
          case 'running':
            color = 'processing';
            text = '运行中';
            break;
          case 'failed':
            color = 'error';
            text = '失败';
            break;
          case 'pending':
            color = 'warning';
            text = '等待中';
            break;
          default:
            color = 'default';
            text = '未知';
        }
        return <Tag color={color}>{text}</Tag>;
      },
    },
    {
      title: '创建时间',
      dataIndex: 'created_at',
      key: 'created_at',
      render: (date: string) => new Date(date).toLocaleString(),
    },
    {
      title: '完成时间',
      dataIndex: 'finished_at',
      key: 'finished_at',
      render: (date: string) => date ? new Date(date).toLocaleString() : '-',
    },
    {
      title: '操作',
      key: 'action',
      render: (text: string, record: any) => (
        <Space size="middle">
          <Button type="link" onClick={() => handleTaskClick(record.id)}>
            查看详情
          </Button>
          {record.status === 'failed' && (
            <Button type="link" onClick={() => electrolyteService.restartCalculation(record.id)}>
              重新运行
            </Button>
          )}
        </Space>
      ),
    },
  ];

  return (
    <div className="task-management-page" style={{ padding: '24px', backgroundColor: '#f5f7fa', minHeight: 'calc(100vh - 64px)' }}>
      <Card style={{ marginBottom: '20px', borderRadius: '8px' }}>
        <Title level={2} style={{ color: '#1890ff' }}>任务管理</Title>
        <Paragraph>
          在这里您可以查看和管理所有计算任务。点击任务可以查看详细信息和结果。
        </Paragraph>
      </Card>

      <Card 
        style={{ borderRadius: '8px' }}
        title="计算任务列表" 
        extra={
          <Button type="primary" onClick={() => navigate('/calculations')}>
            创建新任务
          </Button>
        }
      >
        <Table
          loading={loading}
          dataSource={calculations}
          columns={columns}
          rowKey="id"
          pagination={{
            defaultPageSize: 10,
            showSizeChanger: true,
            showQuickJumper: true,
          }}
        />
      </Card>
    </div>
  );
};

export default TaskManagement; 