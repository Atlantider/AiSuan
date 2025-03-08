import React from 'react';
import { Form, Input, Button, Typography, Card, Row, Col, Checkbox } from 'antd';
import { useNavigate } from 'react-router-dom';
import { UserOutlined, LockOutlined, MailOutlined, TeamOutlined } from '@ant-design/icons';

const { Title, Paragraph } = Typography;

const Register: React.FC = () => {
  const navigate = useNavigate();
  
  const onFinish = (values: any) => {
    console.log('注册信息:', values);
    // 这里应该发送注册请求到后端
    // 临时模拟注册成功
    alert('注册成功！请登录');
    navigate('/login');
  };

  return (
    <div style={{ display: 'flex', justifyContent: 'center', alignItems: 'center', minHeight: '100vh', padding: '20px' }}>
      <Card style={{ width: '100%', maxWidth: 800, boxShadow: '0 4px 12px rgba(0,0,0,0.1)' }}>
        <Row gutter={24}>
          <Col xs={24} md={12} style={{ display: 'flex', flexDirection: 'column', justifyContent: 'center' }}>
            <div style={{ padding: '20px' }}>
              <Title level={2}>创建账号</Title>
              <Paragraph>
                欢迎加入计算材料科学平台，注册账号后即可使用我们的高级计算服务。
              </Paragraph>
              
              <Form
                name="register"
                onFinish={onFinish}
                layout="vertical"
              >
                <Form.Item
                  name="username"
                  rules={[{ required: true, message: '请输入用户名!' }]}
                >
                  <Input prefix={<UserOutlined />} placeholder="用户名" size="large" />
                </Form.Item>
                
                <Form.Item
                  name="email"
                  rules={[
                    { required: true, message: '请输入邮箱!' },
                    { type: 'email', message: '请输入有效的邮箱地址!' }
                  ]}
                >
                  <Input prefix={<MailOutlined />} placeholder="邮箱" size="large" />
                </Form.Item>
                
                <Form.Item
                  name="organization"
                >
                  <Input prefix={<TeamOutlined />} placeholder="所属组织/机构（选填）" size="large" />
                </Form.Item>
                
                <Form.Item
                  name="password"
                  rules={[
                    { required: true, message: '请输入密码!' },
                    { min: 8, message: '密码长度至少为8个字符!' }
                  ]}
                >
                  <Input.Password prefix={<LockOutlined />} placeholder="密码" size="large" />
                </Form.Item>
                
                <Form.Item
                  name="confirm"
                  dependencies={['password']}
                  rules={[
                    { required: true, message: '请确认密码!' },
                    ({ getFieldValue }) => ({
                      validator(_, value) {
                        if (!value || getFieldValue('password') === value) {
                          return Promise.resolve();
                        }
                        return Promise.reject(new Error('两次输入的密码不一致!'));
                      },
                    }),
                  ]}
                >
                  <Input.Password prefix={<LockOutlined />} placeholder="确认密码" size="large" />
                </Form.Item>
                
                <Form.Item
                  name="agreement"
                  valuePropName="checked"
                  rules={[
                    { validator: (_, value) => value ? Promise.resolve() : Promise.reject(new Error('请阅读并同意用户协议')) },
                  ]}
                >
                  <Checkbox>我已阅读并同意<a href="#">用户协议</a>和<a href="#">隐私政策</a></Checkbox>
                </Form.Item>
                
                <Form.Item>
                  <Button type="primary" htmlType="submit" size="large" block>
                    注册
                  </Button>
                </Form.Item>
                
                <div style={{ textAlign: 'center' }}>
                  已有账号？ <a onClick={() => navigate('/login')}>立即登录</a>
                </div>
              </Form>
            </div>
          </Col>
          
          <Col xs={24} md={12} style={{ 
            background: 'linear-gradient(135deg, #1c64f2 0%, #3c83f6 100%)',
            borderRadius: '0 8px 8px 0',
            padding: '40px',
            color: 'white',
            display: 'flex',
            flexDirection: 'column',
            justifyContent: 'center'
          }}>
            <Title level={2} style={{ color: 'white' }}>计算材料科学平台优势</Title>
            <ul style={{ fontSize: '16px', lineHeight: '2', paddingLeft: '20px' }}>
              <li>高效的材料计算和分析工具</li>
              <li>先进的可视化和数据处理功能</li>
              <li>庞大的材料数据库和参考资料</li>
              <li>专业的技术支持和社区资源</li>
              <li>安全可靠的数据存储和共享机制</li>
            </ul>
            <div style={{ marginTop: '30px' }}>
              <Button ghost size="large" onClick={() => navigate('/')}>了解更多</Button>
            </div>
          </Col>
        </Row>
      </Card>
    </div>
  );
};

export default Register; 