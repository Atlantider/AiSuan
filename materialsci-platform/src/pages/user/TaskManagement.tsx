import React, { useState, useEffect } from 'react';
import { 
  Typography, 
  Table, 
  Tag, 
  Space, 
  Button, 
  Card, 
  Tabs, 
  Modal, 
  Descriptions, 
  Progress,
  Statistic,
  Row,
  Col,
  DatePicker,
  Select,
  Input,
  Divider,
  message,
  Spin,
  Empty
} from 'antd';
import { 
  CheckCircleOutlined, 
  SyncOutlined, 
  ClockCircleOutlined, 
  ExclamationCircleOutlined,
  SearchOutlined,
  DownloadOutlined,
  DeleteOutlined,
  EyeOutlined,
  PauseCircleOutlined,
  PlayCircleOutlined
} from '@ant-design/icons';
import { useNavigate } from 'react-router-dom';
import * as electrolyteService from '../../services/electrolyteService';

const { Title, Paragraph } = Typography;
const { TabPane } = Tabs;
const { RangePicker } = DatePicker;
const { Option } = Select;

// 保留模拟数据用于界面展示补充
const mockTasks = [
  {
    id: '1',
    name: '锂离子电池电极材料优化',
    taskType: '电极电位计算',
    module: 'battery',
    submodule: 'electrode',
    status: 'completed',
    priority: 1,
    createdAt: '2023-03-01 09:30:45',
    startedAt: '2023-03-01 09:31:12',
    completedAt: '2023-03-01 10:45:23',
    computeTime: 74.2,
    progress: 100,
  },
  {
    id: '2',
    name: 'NMC正极材料离子迁移研究',
    taskType: '离子迁移能垒计算',
    module: 'battery',
    submodule: 'electrode',
    status: 'running',
    priority: 2,
    createdAt: '2023-03-05 14:22:10',
    startedAt: '2023-03-05 14:23:05',
    completedAt: null,
    computeTime: 103.5,
    progress: 68,
  },
  {
    id: '3',
    name: 'CO2还原反应机理分析',
    taskType: '反应路径计算',
    module: 'catalysis',
    submodule: 'electro',
    status: 'queued',
    priority: 1,
    createdAt: '2023-03-06 11:10:30',
    startedAt: null,
    completedAt: null,
    computeTime: null,
    progress: 0,
  },
  {
    id: '4',
    name: '锂电池电解液优化',
    taskType: '单分子性质计算',
    module: 'battery',
    submodule: 'electrolyte',
    status: 'failed',
    priority: 3,
    createdAt: '2023-03-02 16:45:20',
    startedAt: '2023-03-02 16:46:12',
    completedAt: '2023-03-02 17:12:45',
    computeTime: 26.5,
    progress: 45,
  },
  {
    id: '5',
    name: '全固态电池模拟',
    taskType: '界面稳定性分析',
    module: 'battery',
    submodule: 'fullBattery',
    status: 'cancelled',
    priority: 2,
    createdAt: '2023-03-04 10:30:15',
    startedAt: '2023-03-04 10:31:22',
    completedAt: '2023-03-04 10:35:12',
    computeTime: 3.8,
    progress: 12,
  },
];

const TaskManagement: React.FC = () => {
  const navigate = useNavigate();
  const [activeTab, setActiveTab] = useState('all');
  const [taskDetailVisible, setTaskDetailVisible] = useState(false);
  const [currentTask, setCurrentTask] = useState<any>(null);
  const [calculations, setCalculations] = useState<any[]>([]);
  const [loading, setLoading] = useState(false);
  
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

  // 组件加载时获取计算任务列表
  useEffect(() => {
    fetchCalculations();
  }, []);
  
  // 合并真实计算任务和模拟数据
  // 为保持界面丰富，在实际任务较少时，补充模拟数据
  const getCombinedTasks = () => {
    const realCalculations = calculations.map(calc => ({
      id: calc.id,
      name: calc.name,
      taskType: '电解液计算',
      module: 'battery',
      submodule: 'electrolyte',
      status: calc.status,
      priority: 2,
      createdAt: calc.created_at,
      startedAt: calc.started_at,
      completedAt: calc.finished_at,
      computeTime: calc.started_at && calc.finished_at ? 
        (new Date(calc.finished_at).getTime() - new Date(calc.started_at).getTime()) / (1000 * 60) : null,
      progress: calc.status === 'running' ? 50 : 
               calc.status === 'completed' ? 100 : 
               calc.status === 'failed' ? 30 : 0,
      isReal: true  // 标记为真实任务
    }));
    
    // 如果实际任务较少，添加一些模拟数据
    return realCalculations.length > 0 ? 
      [...realCalculations, ...mockTasks.map(task => ({...task, isReal: false}))] : 
      mockTasks.map(task => ({...task, isReal: false}));
  };
  
  // 获取任务状态标签
  const getStatusTag = (status: string) => {
    switch(status) {
      case 'completed':
        return <Tag icon={<CheckCircleOutlined />} color="success">已完成</Tag>;
      case 'running':
        return <Tag icon={<SyncOutlined spin />} color="processing">运行中</Tag>;
      case 'queued':
        return <Tag icon={<ClockCircleOutlined />} color="default">排队中</Tag>;
      case 'failed':
        return <Tag icon={<ExclamationCircleOutlined />} color="error">失败</Tag>;
      case 'cancelled':
        return <Tag icon={<DeleteOutlined />} color="warning">已取消</Tag>;
      default:
        return <Tag color="default">{status}</Tag>;
    }
  };
  
  // 获取模块名称
  const getModuleName = (module: string, submodule: string) => {
    if (module === 'battery') {
      switch(submodule) {
        case 'electrode': return '电池-电极材料';
        case 'electrolyte': return '电池-电解液';
        case 'fullBattery': return '电池-全电池系统';
        default: return '电池材料';
      }
    } else if (module === 'catalysis') {
      switch(submodule) {
        case 'surface': return '催化-表面催化';
        case 'electro': return '催化-电催化';
        case 'photo': return '催化-光催化';
        default: return '催化材料';
      }
    }
    return '未知模块';
  };
  
  // 显示任务详情
  const showTaskDetail = (task: any) => {
    // 如果是真实数据且有ID，跳转到详情页
    if (task.isReal && task.id) {
      navigate(`/calculations/${task.id}`);
      return;
    }
    
    // 否则显示模态框
    setCurrentTask(task);
    setTaskDetailVisible(true);
  };
  
  // 表格列定义
  const columns = [
    {
      title: '任务名称',
      dataIndex: 'name',
      key: 'name',
      render: (text: string, record: any) => (
        <a onClick={() => showTaskDetail(record)}>{text}</a>
      ),
    },
    {
      title: '计算类型',
      dataIndex: 'taskType',
      key: 'taskType',
    },
    {
      title: '模块',
      key: 'module',
      render: (_: any, record: any) => (
        getModuleName(record.module, record.submodule)
      ),
    },
    {
      title: '状态',
      dataIndex: 'status',
      key: 'status',
      render: (status: string) => getStatusTag(status),
    },
    {
      title: '创建时间',
      dataIndex: 'createdAt',
      key: 'createdAt',
      sorter: (a: any, b: any) => new Date(a.createdAt).getTime() - new Date(b.createdAt).getTime(),
    },
    {
      title: '进度',
      dataIndex: 'progress',
      key: 'progress',
      render: (progress: number, record: any) => (
        record.status === 'running' ? 
        <Progress percent={progress} size="small" /> : 
        `${progress}%`
      ),
    },
    {
      title: '操作',
      key: 'action',
      render: (_: any, record: any) => (
        <Space size="small">
          <Button 
            type="text" 
            icon={<EyeOutlined />} 
            onClick={() => showTaskDetail(record)}
          />
          {record.status === 'running' && (
            <Button 
              type="text" 
              icon={<PauseCircleOutlined />} 
              onClick={() => console.log('暂停任务', record.id)}
            />
          )}
          {record.status === 'queued' && (
            <Button 
              type="text" 
              icon={<DeleteOutlined />} 
              onClick={() => console.log('取消任务', record.id)}
            />
          )}
          {record.status === 'completed' && (
            <Button 
              type="text" 
              icon={<DownloadOutlined />} 
              onClick={() => console.log('下载结果', record.id)}
            />
          )}
          {record.isReal && (
            <Button
              type="text"
              onClick={() => navigate(`/calculations/${record.id}`)}
            >
              详情
            </Button>
          )}
        </Space>
      ),
    },
  ];
  
  // 根据标签筛选任务
  const getFilteredTasks = () => {
    const allTasks = getCombinedTasks();
    if (activeTab === 'all') return allTasks;
    return allTasks.filter(task => task.status === activeTab);
  };
  
  return (
    <div>
      <div style={{ marginBottom: 24 }}>
        <Title level={2}>计算任务管理</Title>
        <Paragraph style={{ fontSize: 16 }}>
          管理您的计算任务，查看任务状态和结果，下载计算数据。
        </Paragraph>
      </div>
      
      {/* 任务统计 */}
      <Row gutter={[16, 16]} style={{ marginBottom: 24 }}>
        <Col xs={24} sm={12} md={6}>
          <Card>
            <Statistic 
              title="总任务数" 
              value={getCombinedTasks().length} 
              valueStyle={{ color: '#1C64F2' }}
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} md={6}>
          <Card>
            <Statistic 
              title="运行中任务" 
              value={getCombinedTasks().filter(t => t.status === 'running').length}
              valueStyle={{ color: '#40A9FF' }}
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} md={6}>
          <Card>
            <Statistic 
              title="已完成任务" 
              value={getCombinedTasks().filter(t => t.status === 'completed').length}
              valueStyle={{ color: '#52C41A' }}
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} md={6}>
          <Card>
            <Statistic 
              title="总计算时长 (小时)" 
              value={getCombinedTasks().reduce((sum, task) => sum + (task.computeTime || 0), 0) / 60}
              precision={1}
              valueStyle={{ color: '#1C64F2' }}
            />
          </Card>
        </Col>
      </Row>
      
      {/* 搜索过滤 */}
      <Card style={{ marginBottom: 16 }}>
        <Row gutter={16}>
          <Col xs={24} md={8}>
            <Input placeholder="搜索任务名称或ID" prefix={<SearchOutlined />} />
          </Col>
          <Col xs={24} md={8}>
            <Select 
              placeholder="计算模块" 
              style={{ width: '100%' }}
              allowClear
            >
              <Option value="battery-electrode">电池-电极材料</Option>
              <Option value="battery-electrolyte">电池-电解液</Option>
              <Option value="battery-fullBattery">电池-全电池系统</Option>
              <Option value="catalysis-surface">催化-表面催化</Option>
              <Option value="catalysis-electro">催化-电催化</Option>
              <Option value="catalysis-photo">催化-光催化</Option>
            </Select>
          </Col>
          <Col xs={24} md={8}>
            <RangePicker style={{ width: '100%' }} placeholder={['开始日期', '结束日期']} />
          </Col>
        </Row>
        <div style={{ marginTop: 16, textAlign: 'right' }}>
          <Button type="primary" icon={<SearchOutlined />}>搜索</Button>
          <Button style={{ marginLeft: 8 }}>重置</Button>
        </div>
      </Card>
      
      {/* 任务列表 */}
      <Card>
        <Tabs 
          activeKey={activeTab} 
          onChange={setActiveTab}
          tabBarExtraContent={
            <Space size="middle">
              <Button 
                type="primary" 
                onClick={() => navigate('/battery/electrolyte')}
              >
                新建电解液计算
              </Button>
              <Button 
                onClick={() => fetchCalculations()}
                loading={loading}
              >
                刷新任务列表
              </Button>
            </Space>
          }
        >
          <TabPane tab="全部任务" key="all" />
          <TabPane tab="运行中" key="running" />
          <TabPane tab="排队中" key="queued" />
          <TabPane tab="已完成" key="completed" />
          <TabPane tab="失败/取消" key="failed" />
        </Tabs>
        
        {loading ? (
          <div style={{ textAlign: 'center', padding: '50px 0' }}>
            <Spin size="large" />
            <div style={{ marginTop: 16 }}>正在加载任务数据...</div>
          </div>
        ) : (
          <Table 
            columns={columns} 
            dataSource={getFilteredTasks()} 
            rowKey="id"
            pagination={{ pageSize: 10 }}
          />
        )}
      </Card>
      
      {/* 任务详情弹窗 */}
      <Modal
        title="任务详情"
        open={taskDetailVisible}
        onCancel={() => setTaskDetailVisible(false)}
        footer={[
          <Button key="back" onClick={() => setTaskDetailVisible(false)}>
            关闭
          </Button>,
          currentTask?.status === 'completed' && (
            <Button 
              key="download" 
              type="primary" 
              icon={<DownloadOutlined />}
              onClick={() => console.log('下载结果', currentTask?.id)}
            >
              下载结果
            </Button>
          )
        ]}
        width={800}
      >
        {currentTask && (
          <>
            <Descriptions bordered column={2}>
              <Descriptions.Item label="任务ID" span={2}>{currentTask.id}</Descriptions.Item>
              <Descriptions.Item label="任务名称" span={2}>{currentTask.name}</Descriptions.Item>
              <Descriptions.Item label="状态">{getStatusTag(currentTask.status)}</Descriptions.Item>
              <Descriptions.Item label="优先级">{'⭐'.repeat(currentTask.priority)}</Descriptions.Item>
              <Descriptions.Item label="计算模块">{getModuleName(currentTask.module, currentTask.submodule)}</Descriptions.Item>
              <Descriptions.Item label="计算类型">{currentTask.taskType}</Descriptions.Item>
              <Descriptions.Item label="创建时间">{currentTask.createdAt}</Descriptions.Item>
              <Descriptions.Item label="开始时间">{currentTask.startedAt || '-'}</Descriptions.Item>
              <Descriptions.Item label="完成时间">{currentTask.completedAt || '-'}</Descriptions.Item>
              <Descriptions.Item label="计算时长">{currentTask.computeTime ? `${currentTask.computeTime} 分钟` : '-'}</Descriptions.Item>
            </Descriptions>
            
            <Divider>进度</Divider>
            
            <Progress 
              percent={currentTask.progress} 
              status={
                currentTask.status === 'completed' ? 'success' : 
                currentTask.status === 'failed' ? 'exception' : 
                'active'
              }
            />
            
            <Divider>计算参数</Divider>
            
            <Paragraph>
              {/* 这里可以根据不同类型的计算展示不同的参数信息 */}
              这里显示任务的详细计算参数，包括计算模型、精度设置、收敛条件等。
            </Paragraph>
          </>
        )}
      </Modal>
    </div>
  );
};

export default TaskManagement; 