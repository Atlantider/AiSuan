import React from 'react';
import { Typography, Collapse, Divider, Card, Tabs, Anchor } from 'antd';

const { Title, Paragraph, Text } = Typography;
const { Panel } = Collapse;
const { TabPane } = Tabs;
const { Link } = Anchor;

const Documentation: React.FC = () => {
  return (
    <div>
      <Title level={2}>文档中心</Title>
      <Paragraph style={{ fontSize: 16 }}>
        欢迎使用计算材料科学平台文档！这里包含了平台的使用指南、计算方法说明和科学参考资料。
      </Paragraph>

      <Divider />

      <div style={{ display: 'flex', gap: 20 }}>
        <div style={{ width: 200 }}>
          <Anchor>
            <Link href="#intro" title="平台简介" />
            <Link href="#tutorial" title="使用教程" />
            <Link href="#technical" title="计算方法" />
            <Link href="#faq" title="常见问题" />
            <Link href="#references" title="参考文献" />
          </Anchor>
        </div>

        <div style={{ flex: 1 }}>
          <section id="intro" style={{ marginBottom: 32 }}>
            <Title level={3}>平台简介</Title>
            <Paragraph>
              计算材料科学平台是一个面向材料科学研究人员的在线计算服务平台，
              提供电池材料、催化材料等多种计算服务。该平台旨在简化复杂的材料计算流程，
              使研究人员能够轻松地进行各种材料性能预测和设计。
            </Paragraph>
            <Card title="平台核心功能" style={{ marginTop: 16 }}>
              <ul>
                <li>多种材料计算模块（电池材料、催化材料等）</li>
                <li>交互式计算工作台</li>
                <li>计算模板自动生成</li>
                <li>计算结果可视化</li>
                <li>材料数据库集成</li>
              </ul>
            </Card>
          </section>

          <section id="tutorial" style={{ marginBottom: 32 }}>
            <Title level={3}>使用教程</Title>
            <Tabs defaultActiveKey="1">
              <TabPane tab="入门指南" key="1">
                <Paragraph>
                  <ol>
                    <li>注册账号并登录系统</li>
                    <li>选择计算模块（电池/催化等）</li>
                    <li>选择具体计算类型（电极材料/电解液等）</li>
                    <li>在工作台选择计算内容（可多选）</li>
                    <li>系统生成计算模板</li>
                    <li>填写模板并提交计算</li>
                    <li>等待系统执行计算并返回结果</li>
                    <li>查看和导出结果</li>
                  </ol>
                </Paragraph>
              </TabPane>
              <TabPane tab="数据导入导出" key="2">
                <Paragraph>
                  计算材料科学平台支持多种格式的材料结构数据导入，包括CIF、POSCAR、XYZ等格式。
                  您可以从本地上传文件，也可以从材料数据库中直接导入结构。
                </Paragraph>
                <Paragraph>
                  计算结果可以导出为多种格式，包括CSV、JSON、图片等，方便用户进行后续分析和发表。
                </Paragraph>
              </TabPane>
              <TabPane tab="高级功能" key="3">
                <Paragraph>
                  <ul>
                    <li><Text strong>参数优化</Text>：系统可以自动优化计算参数，提高计算效率和精度。</li>
                    <li><Text strong>工作流自定义</Text>：高级用户可以自定义计算工作流，实现更复杂的计算任务。</li>
                    <li><Text strong>数据共享</Text>：用户可以选择共享计算结果，促进科研合作。</li>
                  </ul>
                </Paragraph>
              </TabPane>
            </Tabs>
          </section>

          <section id="technical" style={{ marginBottom: 32 }}>
            <Title level={3}>计算方法</Title>
            <Collapse defaultActiveKey={['1']}>
              <Panel header="密度泛函理论计算" key="1">
                <Paragraph>
                  平台集成了常用的DFT计算软件，如VASP、Quantum ESPRESSO等，支持多种交换关联泛函，
                  如PBE、HSE06、B3LYP等。DFT计算主要用于获取材料的电子结构、能带、态密度等信息。
                </Paragraph>
              </Panel>
              <Panel header="分子动力学计算" key="2">
                <Paragraph>
                  平台集成了LAMMPS、GROMACS等分子动力学软件，用于模拟材料在不同温度、压力下的动态行为，
                  可以计算扩散系数、热导率等动力学性质。
                </Paragraph>
              </Panel>
              <Panel header="机器学习辅助计算" key="3">
                <Paragraph>
                  平台集成了多种机器学习方法，可以加速材料性能预测和筛选，包括回归模型、分类模型、
                  主动学习等方法，大幅提高计算效率。
                </Paragraph>
              </Panel>
            </Collapse>
          </section>

          <section id="faq" style={{ marginBottom: 32 }}>
            <Title level={3}>常见问题</Title>
            <Collapse>
              <Panel header="计算需要多长时间？" key="1">
                <Paragraph>
                  计算时间取决于计算类型和系统复杂度。简单的单点能计算可能只需几分钟，
                  而复杂的反应路径计算可能需要几小时甚至几天。系统会在提交计算时给出估计时间。
                </Paragraph>
              </Panel>
              <Panel header="计算结果的精度如何？" key="2">
                <Paragraph>
                  计算精度取决于选择的理论方法和计算参数。平台默认设置旨在平衡计算精度和效率，
                  高级用户可以调整参数以获得更高精度的结果。
                </Paragraph>
              </Panel>
              <Panel header="如何获取技术支持？" key="3">
                <Paragraph>
                  您可以通过"联系我们"页面提交问题，或者查看论坛中的讨论。
                  平台团队会在工作日内回复您的技术问题。
                </Paragraph>
              </Panel>
            </Collapse>
          </section>

          <section id="references" style={{ marginBottom: 32 }}>
            <Title level={3}>参考文献</Title>
            <Paragraph>
              <ol>
                <li>Kohn, W., & Sham, L. J. (1965). Self-consistent equations including exchange and correlation effects. Physical Review, 140(4A), A1133.</li>
                <li>Perdew, J. P., Burke, K., & Ernzerhof, M. (1996). Generalized gradient approximation made simple. Physical Review Letters, 77(18), 3865.</li>
                <li>Henkelman, G., Uberuaga, B. P., & Jónsson, H. (2000). A climbing image nudged elastic band method for finding saddle points and minimum energy paths. The Journal of Chemical Physics, 113(22), 9901-9904.</li>
                <li>Jain, A., Ong, S. P., Hautier, G., Chen, W., Richards, W. D., Dacek, S., ... & Persson, K. A. (2013). Commentary: The Materials Project: A materials genome approach to accelerating materials innovation. APL Materials, 1(1), 011002.</li>
              </ol>
            </Paragraph>
          </section>
        </div>
      </div>
    </div>
  );
};

export default Documentation; 