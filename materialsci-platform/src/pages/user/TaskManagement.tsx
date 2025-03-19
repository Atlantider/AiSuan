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
import type { Dayjs } from 'dayjs';

const { Title, Paragraph } = Typography;
const { TabPane } = Tabs;
const { RangePicker } = DatePicker;
const { Option } = Select;

const TaskManagement: React.FC = () => {
  const navigate = useNavigate();
  const [activeTab, setActiveTab] = useState('all');
  const [taskDetailVisible, setTaskDetailVisible] = useState(false);
  const [currentTask, setCurrentTask] = useState<any>(null);
  const [calculations, setCalculations] = useState<any[]>([]);
  const [loading, setLoading] = useState(false);
  const [searchText, setSearchText] = useState('');
  const [moduleFilter, setModuleFilter] = useState('');
  const [dateRange, setDateRange] = useState<[Dayjs | null, Dayjs | null]>([null, null]);
  
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
  
  // 格式化任务数据以符合表格需要
  const getFormattedTasks = () => {
    return calculations.map(calc => ({
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
    }));
  };
  
  // 获取任务状态标签
  const getStatusTag = (status: string) => {
    switch(status) {
      case 'completed':
        return <Tag icon={<CheckCircleOutlined />} color="success">已完成</Tag>;
      case 'running':
        return <Tag icon={<SyncOutlined spin />} color="processing">运行中</Tag>;
      case 'pending':
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
    // 跳转到详情页
    navigate(`/calculations/${task.id}`);
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
            onClick={() => navigate(`/calculations/${record.id}`)}
          />
          {record.status === 'completed' && (
            <Button 
              type="text" 
              icon={<DownloadOutlined />} 
              onClick={() => console.log('下载结果', record.id)}
            />
          )}
          <Button
            type="link"
            onClick={() => navigate(`/calculations/${record.id}`)}
          >
            详情
          </Button>
        </Space>
      ),
    },
  ];
  
  // 根据标签筛选任务
  const getFilteredTasks = () => {
    let filteredTasks = getFormattedTasks();
    
    // 根据状态标签过滤
    if (activeTab !== 'all') {
      filteredTasks = filteredTasks.filter(task => task.status === activeTab);
    }
    
    // 根据搜索文本过滤
    if (searchText) {
      const lowercaseSearchText = searchText.toLowerCase();
      filteredTasks = filteredTasks.filter(task => 
        task.name.toLowerCase().includes(lowercaseSearchText) || 
        task.id.toString().includes(lowercaseSearchText)
      );
    }
    
    // 根据模块过滤
    if (moduleFilter) {
      const [moduleType, submoduleType] = moduleFilter.split('-');
      filteredTasks = filteredTasks.filter(task => 
        task.module === moduleType && task.submodule === submoduleType
      );
    }
    
    // 根据日期范围过滤
    if (dateRange[0] && dateRange[1]) {
      const startDate = dateRange[0].startOf('day').valueOf();
      const endDate = dateRange[1].endOf('day').valueOf();
      filteredTasks = filteredTasks.filter(task => {
        const taskDate = new Date(task.createdAt).getTime();
        return taskDate >= startDate && taskDate <= endDate;
      });
    }
    
    return filteredTasks;
  };
  
  // 重置过滤条件
  const resetFilters = () => {
    setSearchText('');
    setModuleFilter('');
    setDateRange([null, null]);
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
              value={calculations.length} 
              valueStyle={{ color: '#1C64F2' }}
              loading={loading}
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} md={6}>
          <Card>
            <Statistic 
              title="运行中任务" 
              value={calculations.filter(t => t.status === 'running').length}
              valueStyle={{ color: '#40A9FF' }}
              loading={loading}
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} md={6}>
          <Card>
            <Statistic 
              title="已完成任务" 
              value={calculations.filter(t => t.status === 'completed').length}
              valueStyle={{ color: '#52C41A' }}
              loading={loading}
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} md={6}>
          <Card>
            <Statistic 
              title="总计算时长 (小时)" 
              value={calculations.reduce((sum, task) => {
                if (task.started_at && task.finished_at) {
                  const computeTime = (new Date(task.finished_at).getTime() - new Date(task.started_at).getTime()) / (1000 * 60 * 60);
                  return sum + computeTime;
                }
                return sum;
              }, 0)}
              precision={1}
              valueStyle={{ color: '#1C64F2' }}
              loading={loading}
            />
          </Card>
        </Col>
      </Row>
      
      {/* 搜索过滤 */}
      <Card style={{ marginBottom: 16 }}>
        <Row gutter={16}>
          <Col xs={24} md={8}>
            <Input 
              placeholder="搜索任务名称或ID" 
              prefix={<SearchOutlined />}
              value={searchText}
              onChange={e => setSearchText(e.target.value)}
            />
          </Col>
          <Col xs={24} md={8}>
            <Select 
              placeholder="计算模块" 
              style={{ width: '100%' }}
              allowClear
              value={moduleFilter || undefined}
              onChange={value => setModuleFilter(value)}
            >
              <Option value="battery-electrolyte">电池-电解液</Option>
            </Select>
          </Col>
          <Col xs={24} md={8}>
            <RangePicker 
              style={{ width: '100%' }} 
              placeholder={['开始日期', '结束日期']}
              value={dateRange as any}
              onChange={dates => setDateRange(dates as [Dayjs | null, Dayjs | null])}
            />
          </Col>
        </Row>
        <div style={{ marginTop: 16, textAlign: 'right' }}>
          <Button type="primary" icon={<SearchOutlined />} onClick={() => getFilteredTasks()}>搜索</Button>
          <Button style={{ marginLeft: 8 }} onClick={resetFilters}>重置</Button>
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
                onClick={() => navigate('/calculations', { state: { activeTab: 'modules' } })}
              >
                新建计算任务
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
          <TabPane tab="排队中" key="pending" />
          <TabPane tab="已完成" key="completed" />
          <TabPane tab="失败/取消" key="failed" />
        </Tabs>
        
        {loading ? (
          <div style={{ textAlign: 'center', padding: '50px 0' }}>
            <Spin size="large" />
            <div style={{ marginTop: 16 }}>正在加载任务数据...</div>
          </div>
        ) : calculations.length === 0 ? (
          <Empty 
            description="暂无计算任务，请先创建新的计算任务" 
            style={{ margin: '50px 0' }}
          />
        ) : (
          <Table 
            columns={columns} 
            dataSource={getFilteredTasks()} 
            rowKey="id"
            pagination={{ pageSize: 10 }}
          />
        )}
      </Card>
    </div>
  );
};

export default TaskManagement; 