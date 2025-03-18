import React, { useState, useEffect } from 'react';
import { Typography, Row, Col, Card, Button, Space, Table, Tag, Spin, message, Tabs, Empty } from 'antd';
import { useNavigate, useLocation } from 'react-router-dom';
import { 
  ThunderboltOutlined, 
  ExperimentOutlined, 
  DotChartOutlined,
  LineChartOutlined, 
  ApartmentOutlined,
  ReloadOutlined,
  EyeOutlined,
  CheckCircleOutlined,
  ClockCircleOutlined,
  SyncOutlined,
  CloseCircleOutlined
} from '@ant-design/icons';
import * as electrolyteService from '../../services/electrolyteService';

const { Title, Paragraph } = Typography;
const { TabPane } = Tabs;

const CalculationsIndex: React.FC = () => {
  const navigate = useNavigate();
  const location = useLocation();
  const [activeTab, setActiveTab] = useState('tasks');
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
    
    // 如果是从其他页面带着刷新标志跳转过来的，自动激活任务列表并刷新
    if (location.state && (location.state as any).refresh) {
      setActiveTab('tasks');
      // 清除state，防止重复刷新
      navigate(location.pathname, { replace: true });
    }
  }, [location]);

  // 表格列配置
  const columns = [
    {
      title: '任务名称',
      dataIndex: 'name',
      key: 'name',
      render: (text: string, record: any) => (
        <a onClick={() => navigate(`/calculations/${record.id}`)}>{text}</a>
      )
    },
    {
      title: '状态',
      dataIndex: 'status',
      key: 'status',
      render: (status: string) => {
        let color = 'default';
        let icon = null;
        
        switch (status) {
          case 'pending':
            color = 'warning';
            icon = <ClockCircleOutlined />;
            break;
          case 'running':
            color = 'processing';
            icon = <SyncOutlined spin />;
            break;
          case 'completed':
            color = 'success';
            icon = <CheckCircleOutlined />;
            break;
          case 'failed':
            color = 'error';
            icon = <CloseCircleOutlined />;
            break;
          default:
            break;
        }
        
        return (
          <Tag color={color} icon={icon}>
            {status === 'pending' ? '等待中' : 
             status === 'running' ? '运行中' : 
             status === 'completed' ? '已完成' : 
             status === 'failed' ? '失败' : status}
          </Tag>
        );
      }
    },
    {
      title: '创建时间',
      dataIndex: 'created_at',
      key: 'created_at',
      render: (date: string) => new Date(date).toLocaleString('zh-CN')
    },
    {
      title: '操作',
      key: 'action',
      render: (text: string, record: any) => (
        <Space size="small">
          <Button 
            type="primary" 
            size="small" 
            icon={<EyeOutlined />}
            onClick={() => navigate(`/calculations/${record.id}`)}
          >
            查看
          </Button>
          {record.status === 'failed' && (
            <Button 
              size="small" 
              icon={<ReloadOutlined />}
              onClick={() => handleRestartCalculation(record.id)}
            >
              重新计算
            </Button>
          )}
        </Space>
      )
    }
  ];

  // 计算模块列表
  const calculationModules = [
    {
      title: '电池材料计算',
      description: '提供电极材料、电解液和全电池系统的计算功能，帮助设计和优化高性能电池材料',
      icon: <ThunderboltOutlined style={{ fontSize: 48, color: '#1C64F2' }} />,
      path: '/battery'
    },
    {
      title: '催化材料计算',
      description: '支持表面催化、电催化和光催化反应机理与性能评估，加速催化剂筛选和设计',
      icon: <ExperimentOutlined style={{ fontSize: 48, color: '#1C64F2' }} />,
      path: '/catalysis'
    },
    {
      title: '半导体材料计算',
      description: '计算半导体材料的电子结构、光学性质和载流子动力学，辅助设计新型电子和光电材料',
      icon: <DotChartOutlined style={{ fontSize: 48, color: '#1C64F2' }} />,
      path: '/semiconductor'
    },
    {
      title: '能源材料计算',
      description: '专注于太阳能、热电和储能材料的性能预测和优化，促进可再生能源技术发展',
      icon: <LineChartOutlined style={{ fontSize: 48, color: '#1C64F2' }} />,
      path: '/energy'
    },
    {
      title: '结构材料计算',
      description: '分析结构材料的力学性能、热稳定性和失效机制，支持先进结构材料的开发',
      icon: <ApartmentOutlined style={{ fontSize: 48, color: '#1C64F2' }} />,
      path: '/structural'
    }
  ];

  // 重启失败的计算
  const handleRestartCalculation = async (id: number) => {
    try {
      message.loading('正在重启计算...');
      await electrolyteService.restartCalculation(id);
      message.success('计算任务已重启');
      // 重新获取计算任务列表
      fetchCalculations();
    } catch (error) {
      console.error('重启计算失败:', error);
      message.error('重启计算失败，请稍后重试');
    }
  };

  return (
    <div>
      <Tabs activeKey={activeTab} onChange={setActiveTab}>
        <TabPane tab="计算任务列表" key="tasks">
          <div style={{ marginBottom: 20 }}>
            <Space>
              <Button 
                type="primary" 
                icon={<ReloadOutlined />} 
                onClick={fetchCalculations}
                loading={loading}
              >
                刷新
              </Button>
              <Button onClick={() => navigate('/battery/electrolyte')}>
                新建电解液计算
              </Button>
            </Space>
          </div>
          
          {loading ? (
            <div style={{ textAlign: 'center', padding: '30px 0' }}>
              <Spin size="large" />
            </div>
          ) : (
            calculations.length === 0 ? (
              <Empty 
                description="暂无计算任务" 
                image={Empty.PRESENTED_IMAGE_SIMPLE}
              />
            ) : (
              <Table 
                dataSource={calculations} 
                columns={columns} 
                rowKey="id"
                pagination={{ pageSize: 10 }}
              />
            )
          )}
        </TabPane>
        
        <TabPane tab="计算模块" key="modules">
          <div style={{ marginBottom: 40 }}>
            <Title level={2}>选择计算类型</Title>
            <Paragraph>
              选择您要进行的计算类型，平台将为您提供专业的计算工具和分析能力。
              我们的计算模块涵盖材料科学的多个领域，帮助您高效开展研究工作。
            </Paragraph>
          </div>

          <Row gutter={[24, 24]}>
            {calculationModules.map((module, index) => (
              <Col xs={24} sm={12} md={8} key={index}>
                <Card 
                  hoverable 
                  style={{ height: '100%', display: 'flex', flexDirection: 'column' }}
                  bodyStyle={{ flex: 1, display: 'flex', flexDirection: 'column' }}
                >
                  <div style={{ textAlign: 'center', marginBottom: 20 }}>
                    {module.icon}
                  </div>
                  <Title level={3} style={{ textAlign: 'center', fontSize: 20 }}>
                    {module.title}
                  </Title>
                  <Paragraph style={{ flexGrow: 1 }}>
                    {module.description}
                  </Paragraph>
                  <Button 
                    type="primary" 
                    block 
                    size="large"
                    onClick={() => navigate(module.path)}
                  >
                    进入模块
                  </Button>
                </Card>
              </Col>
            ))}
          </Row>
        </TabPane>
      </Tabs>

      <div style={{ marginTop: 40, textAlign: 'center' }}>
        <Space>
          <Button onClick={() => navigate('/user/dashboard')}>返回控制台</Button>
          <Button type="primary" onClick={() => navigate('/documentation')}>查看使用文档</Button>
        </Space>
      </div>
    </div>
  );
};

export default CalculationsIndex; 