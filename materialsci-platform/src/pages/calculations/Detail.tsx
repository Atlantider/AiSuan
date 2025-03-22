import React, { useState, useEffect } from 'react';
import { useParams, useNavigate } from 'react-router-dom';
import { Typography, Card, Descriptions, Button, Spin, message, Tabs, Space, Tag, Result, Timeline, Divider } from 'antd';
import { 
  CheckCircleOutlined, 
  ClockCircleOutlined, 
  SyncOutlined, 
  CloseCircleOutlined,
  ReloadOutlined,
  ArrowLeftOutlined,
  DownloadOutlined,
  BarChartOutlined,
  FileTextOutlined
} from '@ant-design/icons';
import * as electrolyteService from '../../services/electrolyteService';

const { Title, Paragraph } = Typography;

// 创建Tab项类型
type TabItem = {
  key: string;
  label: React.ReactNode;
  children: React.ReactNode;
};

// 定义计算任务类型
interface ICalculation {
  id: number;
  name: string;
  description: string;
  status: string;
  created_at: string;
  started_at?: string;
  finished_at?: string;
  formulation: number;
  slurm_job_id?: string;
  error_message?: string | Record<string, any>;
}

// 定义计算结果类型
interface ICalculationResult {
  ionic_conductivity?: number;
  density?: number;
  viscosity?: number;
  molar_conductivity?: number;
  diffusion_coefficients?: {
    cation?: number;
    anion?: number;
    overall?: number;
  };
  transport_number?: number;
}

const CalculationDetail: React.FC = () => {
  const { id } = useParams<{ id: string }>();
  const navigate = useNavigate();
  const [calculation, setCalculation] = useState<ICalculation | null>(null);
  const [calculationResult, setCalculationResult] = useState<ICalculationResult | null>(null);
  const [loading, setLoading] = useState(true);
  const [resultLoading, setResultLoading] = useState(false);
  const [activeTab, setActiveTab] = useState('info');
  const [pollingTimer, setPollingTimer] = useState<NodeJS.Timeout | null>(null);

  // 获取计算任务详情
  const fetchCalculation = async () => {
    if (!id) return;
    
    setLoading(true);
    try {
      const response = await electrolyteService.getCalculation(Number(id));
      setCalculation(response.data);
      
      // 如果计算已完成，获取结果
      if (response.data.status === 'completed') {
        fetchCalculationResult();
      }
      // 如果计算正在运行中，启动轮询
      else if (response.data.status === 'running' || response.data.status === 'submitted' || response.data.status === 'queued') {
        startPolling();
      }
    } catch (error) {
      console.error('获取计算任务详情失败:', error);
      message.error('获取计算任务详情失败');
    } finally {
      setLoading(false);
    }
  };

  // 停止轮询
  const stopPolling = () => {
    if (pollingTimer) {
      clearInterval(pollingTimer);
      setPollingTimer(null);
    }
  };

  // 开始轮询
  const startPolling = () => {
    // 清除已有的轮询
    stopPolling();
    // 启动新的轮询
    pollTaskStatus();
  };

  // 轮询计算状态
  const pollTaskStatus = () => {
    const timer = setInterval(async () => {
      if (!id) return;
      
      try {
        const response = await electrolyteService.getCalculationStatus(id);
        const newStatus = response.data?.status;
        
        // 更新计算状态
        setCalculation((prev: ICalculation | null) => prev ? ({ ...prev, status: newStatus }) : null);
        
        // 如果计算已完成或失败，停止轮询并获取结果
        if (newStatus === 'completed' || newStatus === 'failed') {
          stopPolling();
          if (newStatus === 'completed') {
            fetchCalculationResult();
          }
        }
      } catch (error) {
        console.error('轮询计算状态失败:', error);
        stopPolling();
      }
    }, 5000); // 每5秒轮询一次
    
    setPollingTimer(timer);
  };

  // 获取计算结果
  const fetchCalculationResult = async () => {
    if (!id) return;
    
    setResultLoading(true);
    try {
      const response = await electrolyteService.getCalculationResults(id);
      setCalculationResult((prev: ICalculationResult | null) => {
        if (prev && JSON.stringify(prev) === JSON.stringify(response.data)) {
          return prev;
        }
        return response.data;
      });
    } catch (error) {
      console.error('获取计算结果出错:', error);
      message.error('获取计算结果失败，请重试');
    } finally {
      setResultLoading(false);
    }
  };

  // 重启计算
  const handleRestart = async () => {
    if (!id || !calculation) return;
    
    try {
      const response = await electrolyteService.restartCalculation(Number(id));
      message.success('重启计算成功');
      // 更新计算状态
      if (response.data && response.data.status) {
        const newStatus = response.data.status;
        setCalculation((prev: ICalculation | null) => prev ? ({ ...prev, status: newStatus }) : null);
      }
      // 重启成功后导航到任务列表页面
      setTimeout(() => {
        navigate('/user/tasks');
      }, 1500); // 延迟1.5秒跳转，让用户先看到成功消息
    } catch (error) {
      console.error('重启计算失败:', error);
      message.error('重启计算失败，请稍后重试');
    }
  };

  // 手动提交计算
  const handleSubmit = async () => {
    if (!id || !calculation) return;
    
    try {
      await electrolyteService.submitCalculation(Number(id));
      message.success('提交计算成功');
      
      // 更新计算状态为pending
      setCalculation((prev: ICalculation | null) => prev ? ({ ...prev, status: 'pending' }) : null);
      
      // 开始轮询状态
      startPolling();
    } catch (error) {
      console.error('提交计算失败:', error);
      message.error('提交计算失败，请稍后重试');
    }
  };

  // 组件挂载时获取计算任务详情
  useEffect(() => {
    fetchCalculation();
    
    // 组件卸载时清除轮询
    return () => {
      if (pollingTimer) {
        clearInterval(pollingTimer);
      }
    };
  }, [id]);

  // 渲染计算状态标签
  const renderStatusTag = (status: string) => {
    switch (status) {
      case 'pending':
        return <Tag color="blue" icon={<ClockCircleOutlined />}>等待提交</Tag>;
      case 'submitted':
        return <Tag color="orange" icon={<ClockCircleOutlined />}>已提交</Tag>;
      case 'queued':
        return <Tag color="gold" icon={<ClockCircleOutlined />}>排队中</Tag>;
      case 'running':
        return <Tag color="processing" icon={<SyncOutlined spin />}>运行中</Tag>;
      case 'completed':
        return <Tag color="success" icon={<CheckCircleOutlined />}>已完成</Tag>;
      case 'failed':
        return <Tag color="error" icon={<CloseCircleOutlined />}>失败</Tag>;
      case 'cancelled':
        return <Tag color="default" icon={<CloseCircleOutlined />}>已取消</Tag>;
      default:
        return <Tag color="default">{status}</Tag>;
    }
  };

  // 获取计算结果Tab内容
  const getResultTabContent = () => {
    if (!calculation) {
      return <div>加载中...</div>;
    }
    
    if (calculation.status !== 'completed') {
      return (
        <Result
          icon={calculation.status === 'failed' ? <CloseCircleOutlined /> : <ClockCircleOutlined />}
          title={calculation.status === 'failed' ? '计算失败' : '计算未完成'}
          subTitle={calculation.status === 'failed' ? '您可以尝试重新启动计算' : '请等待计算完成后查看结果'}
          extra={calculation.status === 'failed' && (
            <Button type="primary" onClick={handleRestart}>
              重启计算
            </Button>
          )}
        />
      );
    }
    
    if (resultLoading) {
      return (
        <div style={{ textAlign: 'center', padding: '30px 0' }}>
          <Spin size="large" tip="加载计算结果..." />
        </div>
      );
    }
    
    if (!calculationResult) {
      return (
        <Result
          status="warning"
          title="未找到计算结果"
          subTitle="系统找不到此计算任务的结果数据"
          extra={
            <Button type="primary" onClick={fetchCalculationResult}>
              重新获取
            </Button>
          }
        />
      );
    }
    
    // 安全访问扩散系数，避免undefined
    const diffCoef = calculationResult.diffusion_coefficients || {};
    
    return (
      <div>
        <Card title="电解液性能参数" style={{ marginBottom: 16 }}>
          <Descriptions bordered column={2}>
            <Descriptions.Item label="离子电导率 (mS/cm)">
              {calculationResult.ionic_conductivity?.toFixed(4) || '暂无数据'}
            </Descriptions.Item>
            <Descriptions.Item label="密度 (g/cm³)">
              {calculationResult.density?.toFixed(4) || '暂无数据'}
            </Descriptions.Item>
            <Descriptions.Item label="粘度 (mPa·s)">
              {calculationResult.viscosity?.toFixed(4) || '暂无数据'}
            </Descriptions.Item>
            <Descriptions.Item label="摩尔电导率 (S·cm²/mol)">
              {calculationResult.molar_conductivity?.toFixed(4) || '暂无数据'}
            </Descriptions.Item>
          </Descriptions>
        </Card>
        
        <Card title="扩散系数" style={{ marginBottom: 16 }}>
          <Descriptions bordered column={2}>
            <Descriptions.Item label="阳离子扩散系数 (10⁻¹⁰ m²/s)">
              {(diffCoef.cation !== undefined ? (diffCoef.cation * 1e10).toFixed(4) : '暂无数据')}
            </Descriptions.Item>
            <Descriptions.Item label="阴离子扩散系数 (10⁻¹⁰ m²/s)">
              {(diffCoef.anion !== undefined ? (diffCoef.anion * 1e10).toFixed(4) : '暂无数据')}
            </Descriptions.Item>
            <Descriptions.Item label="总体扩散系数 (10⁻¹⁰ m²/s)">
              {(diffCoef.overall !== undefined ? (diffCoef.overall * 1e10).toFixed(4) : '暂无数据')}
            </Descriptions.Item>
            <Descriptions.Item label="迁移数">
              {calculationResult.transport_number?.toFixed(4) || '暂无数据'}
            </Descriptions.Item>
          </Descriptions>
        </Card>
        
        <div style={{ marginTop: 16 }}>
          <Space>
            <Button 
              type="primary" 
              icon={<BarChartOutlined />}
              onClick={() => setActiveTab('visualization')}
            >
              查看可视化结果
            </Button>
            <Button 
              icon={<FileTextOutlined />}
              onClick={() => setActiveTab('log')}
            >
              查看计算日志
            </Button>
            <Button 
              icon={<DownloadOutlined />}
              onClick={() => message.info('下载功能即将上线')}
            >
              下载计算数据
            </Button>
          </Space>
        </div>
      </div>
    );
  };

  // 生成Tab项
  const getTabItems = (): TabItem[] => {
    if (!calculation) {
      return [];
    }
    
    const items: TabItem[] = [
      {
        key: 'info',
        label: '基本信息',
        children: (
          <>
            <Descriptions bordered column={2}>
              <Descriptions.Item label="ID">{calculation.id}</Descriptions.Item>
              <Descriptions.Item label="状态">{renderStatusTag(calculation.status)}</Descriptions.Item>
              <Descriptions.Item label="名称">{calculation.name || '-'}</Descriptions.Item>
              <Descriptions.Item label="描述">{calculation.description || '无'}</Descriptions.Item>
              <Descriptions.Item label="创建时间">
                {calculation.created_at ? new Date(calculation.created_at).toLocaleString('zh-CN') : '-'}
              </Descriptions.Item>
              <Descriptions.Item label="配方ID">
                {typeof calculation.formulation === 'object' 
                  ? JSON.stringify(calculation.formulation)
                  : (calculation.formulation || '-')}
              </Descriptions.Item>
              {calculation.started_at && (
                <Descriptions.Item label="开始计算时间">
                  {new Date(calculation.started_at).toLocaleString('zh-CN')}
                </Descriptions.Item>
              )}
              {calculation.finished_at && (
                <Descriptions.Item label="计算完成时间">
                  {new Date(calculation.finished_at).toLocaleString('zh-CN')}
                </Descriptions.Item>
              )}
              {calculation.error_message && (
                <Descriptions.Item label="错误信息" span={2}>
                  {typeof calculation.error_message === 'object' 
                    ? JSON.stringify(calculation.error_message, null, 2)
                    : calculation.error_message}
                </Descriptions.Item>
              )}
            </Descriptions>
            
            <Divider>计算过程</Divider>
            
            <Timeline mode="left">
              <Timeline.Item color="blue">
                创建计算任务 - {calculation.created_at ? new Date(calculation.created_at).toLocaleString('zh-CN') : '-'}
              </Timeline.Item>
              {calculation.started_at && (
                <Timeline.Item color="blue">
                  开始计算 - {new Date(calculation.started_at).toLocaleString('zh-CN')}
                </Timeline.Item>
              )}
              {calculation.status === 'running' && (
                <Timeline.Item color="blue" dot={<SyncOutlined spin />}>
                  正在计算中...
                </Timeline.Item>
              )}
              {calculation.finished_at && (
                <Timeline.Item 
                  color={calculation.status === 'completed' ? 'green' : 'red'}
                  dot={calculation.status === 'completed' ? <CheckCircleOutlined /> : <CloseCircleOutlined />}
                >
                  计算{calculation.status === 'completed' ? '完成' : '失败'} - {new Date(calculation.finished_at).toLocaleString('zh-CN')}
                </Timeline.Item>
              )}
            </Timeline>
          </>
        )
      },
      {
        key: 'result',
        label: '计算结果',
        children: getResultTabContent()
      },
      {
        key: 'visualization',
        label: '可视化结果',
        children: (
          <Result
            icon={<BarChartOutlined />}
            title="可视化功能正在开发中"
            subTitle="请期待后续版本更新"
          />
        )
      },
      {
        key: 'log',
        label: '计算日志',
        children: (
          <Result
            icon={<FileTextOutlined />}
            title="计算日志功能正在开发中"
            subTitle="请期待后续版本更新"
          />
        )
      }
    ];
    
    return items;
  };

  // 渲染页面
  return (
    <div style={{ padding: 24, backgroundColor: '#f5f7fa', minHeight: 'calc(100vh - 64px)' }}>
      {loading ? (
        <div style={{ textAlign: 'center', padding: '50px 0' }}>
          <Spin size="large" tip="加载计算任务详情..." />
        </div>
      ) : !calculation ? (
        <Result
          status="error"
          title="未找到计算任务"
          subTitle={`无法找到ID为 ${id} 的计算任务`}
          extra={
            <Button type="primary" onClick={() => navigate('/user/tasks')}>
              返回计算任务列表
            </Button>
          }
        />
      ) : (
        <>
          <div style={{ marginBottom: 16 }}>
            <Button 
              type="primary" 
              icon={<ArrowLeftOutlined />} 
              onClick={() => navigate('/user/tasks')}
            >
              返回计算任务列表
            </Button>
          </div>
          
          <Card className="detail-card" style={{ boxShadow: '0 1px 4px rgba(0,0,0,0.1)', borderRadius: '8px' }}>
            <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 24 }}>
              <div>
                <Title level={3} style={{ marginBottom: 8, color: '#1890ff' }}>{calculation.name || '-'}</Title>
                <Space>
                  {renderStatusTag(calculation.status)}
                  <span>创建时间: {calculation.created_at ? new Date(calculation.created_at).toLocaleString('zh-CN') : '-'}</span>
                </Space>
              </div>
              
              <Space>
                {calculation.status === 'failed' && (
                  <Button 
                    type="primary" 
                    danger
                    icon={<ReloadOutlined />} 
                    onClick={handleRestart}
                  >
                    重启计算
                  </Button>
                )}
              </Space>
            </div>
            
            <Tabs 
              activeKey={activeTab} 
              onChange={setActiveTab}
              items={getTabItems()}
              className="detail-tabs"
              type="card"
            />
          </Card>
        </>
      )}
    </div>
  );
};

export default CalculationDetail; 