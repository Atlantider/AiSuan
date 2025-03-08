import React from 'react';
import { Form, Input, Button, Checkbox, Typography, Card, Row, Col } from 'antd';
import { UserOutlined, LockOutlined } from '@ant-design/icons';
import { useNavigate } from 'react-router-dom';
import { useDispatch } from 'react-redux';
import { setLoginState, setUserInfo, setToken } from '../store/slices/userSlice';

const { Title, Paragraph } = Typography;

const Login: React.FC = () => {
  const navigate = useNavigate();
  const dispatch = useDispatch();

  const onFinish = (values: any) => {
    console.log('登录信息:', values);
    
    // 假设这是成功的登录响应
    const mockUserInfo = {
      id: '12345',
      username: values.username,
      email: 'user@example.com',
      role: 'regular' as 'admin' | 'advanced' | 'regular',
      avatar: 'https://joeschmoe.io/api/v1/random',
    };
    
    const mockToken = 'mock_jwt_token_' + Math.random().toString(36).substring(2);
    
    // 存储到本地存储，以便在页面刷新后保持登录状态
    localStorage.setItem('token', mockToken);
    localStorage.setItem('userInfo', JSON.stringify(mockUserInfo));
    
    // 更新Redux状态
    dispatch(setToken(mockToken));
    dispatch(setUserInfo(mockUserInfo));
    dispatch(setLoginState(true));
    
    // 导航到首页
    navigate('/');
  };

  return (
    <div style={{ display: 'flex', justifyContent: 'center', alignItems: 'center', minHeight: '100vh', padding: '20px' }}>
      <Card style={{ width: '100%', maxWidth: 800, boxShadow: '0 4px 12px rgba(0,0,0,0.1)' }}>
        <Row gutter={24}>
          <Col xs={24} md={12} style={{ 
            background: 'linear-gradient(135deg, #1c64f2 0%, #3c83f6 100%)',
            borderRadius: '8px 0 0 8px',
            padding: '40px',
            color: 'white',
            display: 'flex',
            flexDirection: 'column',
            justifyContent: 'center'
          }}>
            <Title level={2} style={{ color: 'white' }}>欢迎回来</Title>
            <Paragraph style={{ color: 'white', fontSize: 16 }}>
              登录计算材料科学平台，开始您的材料计算之旅，探索前沿材料设计与性能预测。
            </Paragraph>
            <div style={{ marginTop: 40 }}>
              <Title level={4} style={{ color: 'white' }}>还没有账号？</Title>
              <Button ghost size="large" onClick={() => navigate('/register')}>
                立即注册
              </Button>
            </div>
          </Col>
          
          <Col xs={24} md={12} style={{ display: 'flex', flexDirection: 'column', justifyContent: 'center' }}>
            <div style={{ padding: '20px 40px' }}>
              <Title level={2}>用户登录</Title>
              <Form
                name="login"
                initialValues={{ remember: true }}
                onFinish={onFinish}
                size="large"
              >
                <Form.Item
                  name="username"
                  rules={[{ required: true, message: '请输入用户名!' }]}
                >
                  <Input prefix={<UserOutlined />} placeholder="用户名" />
                </Form.Item>

                <Form.Item
                  name="password"
                  rules={[{ required: true, message: '请输入密码!' }]}
                >
                  <Input.Password prefix={<LockOutlined />} placeholder="密码" />
                </Form.Item>

                <Form.Item>
                  <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
                    <Form.Item name="remember" valuePropName="checked" noStyle>
                      <Checkbox>记住我</Checkbox>
                    </Form.Item>
                    <a href="#">忘记密码?</a>
                  </div>
                </Form.Item>

                <Form.Item>
                  <Button type="primary" htmlType="submit" block>
                    登录
                  </Button>
                </Form.Item>
                
                <div style={{ textAlign: 'center', marginTop: 20 }}>
                  <Button type="link" onClick={() => navigate('/')}>
                    返回首页
                  </Button>
                </div>
              </Form>
            </div>
          </Col>
        </Row>
      </Card>
    </div>
  );
};

export default Login; 