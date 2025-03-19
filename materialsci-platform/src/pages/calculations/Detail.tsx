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
const { TabPane } = Tabs;

const CalculationDetail: React.FC = () => {
  const { id } = useParams<{ id: string }>();
  const navigate = useNavigate();
  const [calculation, setCalculation] = useState<any>(null);
  const [calculationResult, setCalculationResult] = useState<any>(null);
  const [loading, setLoading] = useState(true);
  const [resultLoading, setResultLoading] = useState(false);
  const [activeTab, setActiveTab] = useState('info');
  const [pollingTimer, setPollingTimer] = useState<NodeJS.Timeout | null>(null);

  // 获取计算任务详情
  const fetchCalculation = async () => {
    if (!id) return;
    
    try {
      const response = await electrolyteService.getCalculation(Number(id));
      setCalculation(response.data);
      
      // 如果计算任务状态为running或pending，开始轮询
      if (response.data?.status === 'running' || response.data?.status === 'pending') {
        startPolling();
      } else if (response.data?.status === 'completed') {
        // 如果计算已完成，获取结果
        fetchCalculationResult();
      }
    } catch (error) {
      console.error('获取计算任务详情失败:', error);
      message.error('获取计算任务详情失败，请稍后重试');
    } finally {
      setLoading(false);
    }
  };

  // 获取计算结果
  const fetchCalculationResult = async () => {
    if (!id) return;
    
    try {
      setResultLoading(true);
      const response = await electrolyteService.getCalculationResults(id);
      setCalculationResult(response.data);
    } catch (error) {
      console.error('获取计算结果失败:', error);
      message.error('获取计算结果失败，请稍后重试');
    } finally {
      setResultLoading(false);
    }
  };

  // 轮询计算状态
  const startPolling = () => {
    if (pollingTimer) {
      clearInterval(pollingTimer);
    }
    
    const timer = setInterval(async () => {
      try {
        const response = await electrolyteService.getCalculationStatus(id!);
        const newStatus = response.data?.status;
        
        // 更新计算状态
        setCalculation(prev => ({ ...prev, status: newStatus }));
        
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

  // 停止轮询
  const stopPolling = () => {
    if (pollingTimer) {
      clearInterval(pollingTimer);
      setPollingTimer(null);
    }
  };

  // 重启计算
  const handleRestart = async () => {
    if (!id) return;
    
    try {
      message.loading('正在重启计算...');
      await electrolyteService.restartCalculation(Number(id));
      message.success('计算任务已重启');
      
      // 更新计算状态为pending
      setCalculation(prev => ({ ...prev, status: 'pending' }));
      
      // 开始轮询状态
      startPolling();
    } catch (error) {
      console.error('重启计算失败:', error);
      message.error('重启计算失败，请稍后重试');
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
    let color = 'default';
    let icon = null;
    let text = status;
    
    switch (status) {
      case 'pending':
        color = 'warning';
        icon = <ClockCircleOutlined />;
        text = '等待中';
        break;
      case 'running':
        color = 'processing';
        icon = <SyncOutlined spin />;
        text = '运行中';
        break;
      case 'completed':
        color = 'success';
        icon = <CheckCircleOutlined />;
        text = '已完成';
        break;
      case 'failed':
        color = 'error';
        icon = <CloseCircleOutlined />;
        text = '失败';
        break;
      default:
        break;
    }
    
    return <Tag color={color} icon={icon}>{text}</Tag>;
  };

  if (loading) {
    return (
      <div style={{ display: 'flex', justifyContent: 'center', alignItems: 'center', height: '50vh' }}>
        <Spin size="large" tip="加载中..." />
      </div>
    );
  }

  if (!calculation) {
    return (
      <Result
        status="404"
        title="计算任务不存在"
        subTitle="您请求的计算任务不存在或已被删除"
        extra={
          <Button type="primary" onClick={() => navigate('/calculations')}>
            返回计算任务列表
          </Button>
        }
      />
    );
  }

  return (
    <div style={{ padding: '24px 0' }}>
      <div style={{ marginBottom: 16 }}>
        <Button 
          icon={<ArrowLeftOutlined />} 
          onClick={() => navigate('/calculations')}
        >
          返回计算任务列表
        </Button>
      </div>
      
      <Card>
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 24 }}>
          <div>
            <Title level={3} style={{ marginBottom: 8 }}>{calculation.name}</Title>
            <Space>
              {renderStatusTag(calculation.status)}
              <span>创建时间: {new Date(calculation.created_at).toLocaleString('zh-CN')}</span>
            </Space>
          </div>
          
          <Space>
            {calculation.status === 'failed' && (
              <Button 
                type="primary" 
                icon={<ReloadOutlined />} 
                onClick={handleRestart}
              >
                重启计算
              </Button>
            )}
          </Space>
        </div>
        
        <Tabs activeKey={activeTab} onChange={setActiveTab}>
          <TabPane tab="基本信息" key="info">
            <Descriptions bordered column={2}>
              <Descriptions.Item label="ID">{calculation.id}</Descriptions.Item>
              <Descriptions.Item label="状态">{renderStatusTag(calculation.status)}</Descriptions.Item>
              <Descriptions.Item label="名称">{calculation.name}</Descriptions.Item>
              <Descriptions.Item label="描述">{calculation.description || '无'}</Descriptions.Item>
              <Descriptions.Item label="创建时间">{new Date(calculation.created_at).toLocaleString('zh-CN')}</Descriptions.Item>
              <Descriptions.Item label="配方ID">{calculation.formulation}</Descriptions.Item>
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
                  {calculation.error_message}
                </Descriptions.Item>
              )}
            </Descriptions>
            
            <Divider>计算过程</Divider>
            
            <Timeline mode="left">
              <Timeline.Item color="blue">
                创建计算任务 - {new Date(calculation.created_at).toLocaleString('zh-CN')}
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
          </TabPane>
          
          <TabPane tab="计算结果" key="result">
            {calculation.status !== 'completed' ? (
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
            ) : resultLoading ? (
              <div style={{ textAlign: 'center', padding: '30px 0' }}>
                <Spin size="large" tip="加载计算结果..." />
              </div>
            ) : !calculationResult ? (
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
            ) : (
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
                      {(calculationResult.diffusion_coefficients?.cation * 1e10)?.toFixed(4) || '暂无数据'}
                    </Descriptions.Item>
                    <Descriptions.Item label="阴离子扩散系数 (10⁻¹⁰ m²/s)">
                      {(calculationResult.diffusion_coefficients?.anion * 1e10)?.toFixed(4) || '暂无数据'}
                    </Descriptions.Item>
                    <Descriptions.Item label="总体扩散系数 (10⁻¹⁰ m²/s)">
                      {(calculationResult.diffusion_coefficients?.overall * 1e10)?.toFixed(4) || '暂无数据'}
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
            )}
          </TabPane>
          
          <TabPane tab="可视化结果" key="visualization">
            <Result
              icon={<BarChartOutlined />}
              title="可视化功能正在开发中"
              subTitle="请期待后续版本更新"
            />
          </TabPane>
          
          <TabPane tab="计算日志" key="log">
            <Result
              icon={<FileTextOutlined />}
              title="计算日志功能正在开发中"
              subTitle="请期待后续版本更新"
            />
          </TabPane>
        </Tabs>
      </Card>
    </div>
  );
};

export default CalculationDetail; 