import React, { useEffect, useRef, useState } from 'react';
import { Spin, Button, Space, Radio, Tooltip } from 'antd';
import './CrystalViewer.css';

// 引入3Dmol脚本
// 注意：需要在public/index.html中添加3Dmol.js的CDN链接，因为npm版本可能不完整
// <script src="https://3Dmol.org/build/3Dmol-min.js"></script>

// 修正Window接口的类型定义
declare global {
  interface Window {
    $3Dmol: typeof $3Dmol;
  }
}

interface AtomPosition {
  element: string;
  x: number;
  y: number;
  z: number;
  charge?: string;
}

interface CrystalViewerProps {
  atoms?: AtomPosition[];
  width?: number;
  height?: number;
  style?: 'stick' | 'sphere' | 'cartoon' | 'line' | 'polyhedra';
  colorScheme?: 'element' | 'spectrum' | 'charge';
  showUnitCell?: boolean;
  showLabels?: boolean;
  backgroundColor?: string;
  rotationSpeed?: number;
  latticeParams?: {
    a: number;
    b: number;
    c: number;
    alpha: number;
    beta: number;
    gamma: number;
  };
}

// 元素颜色映射 - 标准配色，用于生成颜色定义
const elementColors: Record<string, string> = {
  'H': '#FFFFFF',
  'Li': '#8F40D4',
  'C': '#909090',
  'N': '#3050F8',
  'O': '#FF0D0D',
  'Na': '#AB5CF2',
  'Mg': '#8AFF00',
  'Al': '#BFA6A6',
  'Si': '#F0C8A0',
  'P': '#FF8000',
  'S': '#FFFF30',
  'K': '#8F40D4',
  'Ca': '#3DFF00',
  'Fe': '#E06633',
  'Co': '#0000FF',
  'Ni': '#50D050',
  'Cu': '#C88033',
  'Zn': '#7D80B0',
};

const CrystalViewer: React.FC<CrystalViewerProps> = ({
  atoms = [],
  width = 600,
  height = 400,
  style = 'polyhedra',
  colorScheme = 'element',
  showUnitCell = true,
  showLabels = true,
  backgroundColor = '#ffffff',
  rotationSpeed = 0.5,
  latticeParams = { a: 5.4307, b: 5.4307, c: 5.4307, alpha: 90, beta: 90, gamma: 90 }
}) => {
  const containerRef = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<any>(null);
  const [loading, setLoading] = useState(true);
  const [scriptLoaded, setScriptLoaded] = useState(false);
  const [uniqueElements, setUniqueElements] = useState<{element: string, color: string, charge?: string}[]>([]);
  const [currentStyle, setCurrentStyle] = useState<string>(style);
  const [isRotating, setIsRotating] = useState<boolean>(rotationSpeed > 0);

  // 加载3Dmol.js脚本
  useEffect(() => {
    if (window.$3Dmol) {
      setScriptLoaded(true);
      return;
    }

    const script = document.createElement('script');
    script.src = 'https://3Dmol.org/build/3Dmol-min.js';
    script.async = true;
    script.onload = () => {
      console.log('3Dmol.js loaded successfully');
      setScriptLoaded(true);
    };
    script.onerror = () => {
      console.error('Failed to load 3Dmol.js');
    };
    document.body.appendChild(script);

    return () => {
      document.body.removeChild(script);
    };
  }, []);

  // 收集唯一元素并确保图例颜色与显示匹配
  useEffect(() => {
    if (atoms.length > 0) {
      const elements = new Set<string>();
      const elementCharges = new Map<string, string>();
      
      atoms.forEach(atom => {
        const elementWithCharge = atom.charge ? 
          `${atom.element}${atom.charge.replace(/[^0-9+-]+/g, '')}` : 
          atom.element;
          
        elements.add(elementWithCharge);
        if (atom.charge) {
          elementCharges.set(elementWithCharge, atom.charge);
        }
      });
      
      const legendItems: {element: string, color: string, charge?: string}[] = [];
      elements.forEach(element => {
        const baseElement = element.replace(/[0-9+-]+/g, '');
        legendItems.push({
          element,
          color: elementColors[baseElement] || '#808080'
        });
      });
      
      setUniqueElements(legendItems);
    }
  }, [atoms]);

  // 当脚本加载完成和原子数据可用时初始化查看器
  useEffect(() => {
    if (!scriptLoaded || !containerRef.current || atoms.length === 0) return;
    
    // 清除现有内容
    const molContainer = containerRef.current.querySelector('.mol-container');
    if (molContainer) {
      containerRef.current.removeChild(molContainer);
    }
    
    // 创建容器元素
    const viewerDiv = document.createElement('div');
    viewerDiv.className = 'mol-container';
    viewerDiv.style.width = `${width}px`;
    viewerDiv.style.height = `${height}px`;
    containerRef.current.appendChild(viewerDiv);
    
    setLoading(true);
    
    // 初始化查看器
    try {
      const viewer = window.$3Dmol.createViewer(viewerDiv, {
        backgroundColor: backgroundColor
      });
      viewerRef.current = viewer;
      
      // 生成PDB格式数据
      const pdbData = generatePDB(atoms, latticeParams);
      
      // 添加模型
      const model = viewer.addModel(pdbData, 'pdb');
      
      // 应用样式
      applyStyle(model, currentStyle, colorScheme);
      
      // 显示晶胞
      if (showUnitCell) {
        addUnitCell(viewer, latticeParams);
      }
      
      // 显示标签
      if (showLabels) {
        addLabels(viewer, atoms);
      }
      
      // 缩放适应
      viewer.zoomTo();
      
      // 设置自动旋转
      if (isRotating) {
        viewer.spin(rotationSpeed < 1 ? rotationSpeed : 1);
      }
      
      // 渲染
      viewer.render();
      setLoading(false);
    } catch (error) {
      console.error('Error initializing 3Dmol viewer:', error);
      setLoading(false);
    }
    
    // 清理函数
    return () => {
      if (viewerRef.current) {
        try {
          viewerRef.current.spin(false);
          viewerRef.current = null;
        } catch (error) {
          console.error('Error cleaning up viewer:', error);
        }
      }
    };
  }, [scriptLoaded, atoms, currentStyle, colorScheme, showUnitCell, showLabels, backgroundColor, rotationSpeed, latticeParams, width, height, isRotating]);

  // 生成PDB格式的结构数据
  const generatePDB = (atoms: AtomPosition[], lattice: any) => {
    let pdb = 'HEADER Created by AiSuan Materials Science Platform\n';
    
    // 添加晶格信息
    const { a, b, c, alpha, beta, gamma } = lattice;
    pdb += `CRYST1${formatPDBFloat(a, 9)}${formatPDBFloat(b, 9)}${formatPDBFloat(c, 9)}${formatPDBFloat(alpha, 7)}${formatPDBFloat(beta, 7)}${formatPDBFloat(gamma, 7)} P 1\n`;
    
    // 添加原子信息
    atoms.forEach((atom, idx) => {
      const i = idx + 1;
      const element = atom.element.padEnd(2);
      
      pdb += `ATOM  ${i.toString().padStart(5)} ${element} MAT A   1    ${formatPDBFloat(atom.x, 8)}${formatPDBFloat(atom.y, 8)}${formatPDBFloat(atom.z, 8)}  1.00  0.00           ${element.trim()}\n`;
    });
    
    pdb += 'END\n';
    return pdb;
  };
  
  // 格式化PDB浮点数
  const formatPDBFloat = (num: number, width: number) => {
    const str = num.toFixed(3);
    return str.padStart(width);
  };
  
  // 应用样式到模型
  const applyStyle = (model: any, styleType: string, colorType: string) => {
    let styleObj: any = {};
    
    // 根据样式类型设置
    switch (styleType) {
      case 'stick':
        styleObj.stick = { radius: 0.2 };
        break;
      case 'sphere':
        styleObj.sphere = { scale: 0.5 };
        break;
      case 'cartoon':
        styleObj.cartoon = {};
        break;
      case 'line':
        styleObj.line = {};
        break;
      case 'polyhedra':
        // 多面体样式，结合球体和棍棒
        styleObj.sphere = { scale: 0.5 };
        styleObj.stick = { radius: 0.15 };
        
        // 尝试添加多面体效果
        try {
          if (model.addPolyhedra) {
            // 识别中心原子（如过渡金属）和配位原子（如氧）
            const metalElements = ['Li', 'Na', 'K', 'Mg', 'Ca', 'Fe', 'Co', 'Ni', 'Mn', 'Cr', 'V', 'Ti', 'Zn', 'Cu'];
            const ligandElements = ['O', 'S', 'N', 'F', 'Cl'];
            
            // 为金属原子周围添加多面体
            model.addPolyhedra({
              center: {elem: metalElements},
              neighbors: {elem: ligandElements},
              color: 'rgba(255,165,0,0.5)' // 半透明橙色多面体
            });
          }
        } catch (error) {
          console.warn('多面体渲染失败，回退到默认样式', error);
        }
        break;
      default:
        styleObj.stick = { radius: 0.2 };
        styleObj.sphere = { scale: 0.4 };
    }
    
    // 根据颜色方案设置
    if (colorType === 'element') {
      styleObj.colors = elementColors;
    } else if (colorType === 'spectrum') {
      // 使用光谱颜色
      model.setStyle({}, styleObj);
      model.setColorByProperty({}, 'index', new window.$3Dmol.Gradient.RWB(0, atoms.length));
      return;
    } else if (colorType === 'charge') {
      // 根据电荷着色
      styleObj.colorScheme = 'Jmol';
    }
    
    model.setStyle({}, styleObj);
  };
  
  // 添加单元晶胞
  const addUnitCell = (viewer: any, lattice: any) => {
    const { a, b, c, alpha, beta, gamma } = lattice;
    
    // 转换角度为弧度
    const alphaRad = (alpha * Math.PI) / 180;
    const betaRad = (beta * Math.PI) / 180;
    const gammaRad = (gamma * Math.PI) / 180;
    
    // 计算晶格矢量
    const ax = a;
    const ay = 0;
    const az = 0;
    
    const bx = b * Math.cos(gammaRad);
    const by = b * Math.sin(gammaRad);
    const bz = 0;
    
    const cx = c * Math.cos(betaRad);
    const cy = (c * Math.cos(alphaRad) - cx * Math.cos(gammaRad)) / Math.sin(gammaRad);
    const cz = Math.sqrt(c*c - cx*cx - cy*cy);
    
    // 定义晶胞顶点
    const v = [
      {x: 0, y: 0, z: 0},
      {x: ax, y: ay, z: az},
      {x: bx, y: by, z: bz},
      {x: ax+bx, y: ay+by, z: az+bz},
      
      {x: cx, y: cy, z: cz},
      {x: ax+cx, y: ay+cy, z: az+cz},
      {x: bx+cx, y: by+cy, z: bz+cz},
      {x: ax+bx+cx, y: ay+by+cy, z: az+bz+cz}
    ];
    
    // 添加晶胞边线
    const unitCellLines = [
      [0,1], [1,3], [3,2], [2,0],
      [4,5], [5,7], [7,6], [6,4],
      [0,4], [1,5], [2,6], [3,7]
    ];
    
    unitCellLines.forEach(([i, j]) => {
      viewer.addLine({
        start: v[i], 
        end: v[j],
        color: '0x000000',
        linewidth: 1.5
      });
    });
    
    // 添加晶轴标签
    if (showLabels) {
      // 使用箭头表示晶轴方向
      viewer.addArrow({
        start: v[0],
        end: {x: ax*1.1, y: 0, z: 0},
        color: 'red',
        radius: 0.1
      });
      viewer.addArrow({
        start: v[0],
        end: {x: bx*1.1, y: by*1.1, z: 0},
        color: 'green',
        radius: 0.1
      });
      viewer.addArrow({
        start: v[0],
        end: {x: cx*1.1, y: cy*1.1, z: cz*1.1},
        color: 'blue',
        radius: 0.1
      });
      
      // 添加文字标签
      viewer.addLabel("a", {position: {x: ax*1.15, y: 0, z: 0}, backgroundColor: "transparent", fontColor: "red", fontSize: 14});
      viewer.addLabel("b", {position: {x: bx*1.15, y: by*1.15, z: 0}, backgroundColor: "transparent", fontColor: "green", fontSize: 14});
      viewer.addLabel("c", {position: {x: cx*1.15, y: cy*1.15, z: cz*1.15}, backgroundColor: "transparent", fontColor: "blue", fontSize: 14});
    }
  };
  
  // 添加原子标签
  const addLabels = (viewer: any, atoms: AtomPosition[]) => {
    atoms.forEach((atom, i) => {
      // 只为少量原子添加标签，避免过度拥挤
      if (i % 5 === 0) {
        const label = atom.charge ? 
          `${atom.element}${atom.charge.replace(/[^0-9+-]+/g, '')}` : 
          atom.element;
          
        viewer.addLabel(label, {
          position: {x: atom.x, y: atom.y, z: atom.z},
          backgroundColor: "rgba(255, 255, 255, 0.5)",
          fontColor: "black",
          fontSize: 12
        });
      }
    });
  };

  // 更改显示方式
  const changeStyle = (newStyle: string) => {
    setCurrentStyle(newStyle);
    
    if (viewerRef.current) {
      const model = viewerRef.current.getModel();
      if (model) {
        applyStyle(model, newStyle, colorScheme);
        viewerRef.current.render();
      }
    }
  };

  // 切换旋转状态
  const toggleRotation = () => {
    setIsRotating(!isRotating);
    
    if (viewerRef.current) {
      if (!isRotating) {
        viewerRef.current.spin(rotationSpeed < 1 ? rotationSpeed : 1);
      } else {
        viewerRef.current.spin(false);
      }
    }
  };

  // 没有原子数据时显示提示
  if (atoms.length === 0) {
    return (
      <div className="crystal-viewer-container" style={{ width, height }} ref={containerRef}>
        <div className="loading-container">
          <p style={{ color: '#888', fontStyle: 'italic' }}>等待结构数据...</p>
        </div>
      </div>
    );
  }

  return (
    <div className="crystal-viewer-container" style={{ width, height }} ref={containerRef}>
      {loading && (
        <div className="loading-container">
          <Spin tip="加载结构..." />
        </div>
      )}
      
      {/* 控制面板 */}
      {!loading && (
        <div className="viewer-controls">
          <Space>
            <Radio.Group 
              value={currentStyle} 
              onChange={(e) => changeStyle(e.target.value)}
              size="small"
              optionType="button"
              buttonStyle="solid"
            >
              <Tooltip title="多面体">
                <Radio.Button value="polyhedra">多面体</Radio.Button>
              </Tooltip>
              <Tooltip title="球棍式">
                <Radio.Button value="stick">球棍</Radio.Button>
              </Tooltip>
              <Tooltip title="球式">
                <Radio.Button value="sphere">球式</Radio.Button>
              </Tooltip>
              <Tooltip title="线框">
                <Radio.Button value="line">线框</Radio.Button>
              </Tooltip>
            </Radio.Group>
            <Tooltip title={isRotating ? "停止旋转" : "开始旋转"}>
              <Button 
                size="small" 
                type={isRotating ? "primary" : "default"}
                onClick={toggleRotation}
              >
                {isRotating ? "停止旋转" : "开始旋转"}
              </Button>
            </Tooltip>
          </Space>
        </div>
      )}
      
      {/* 元素图例 */}
      {!loading && uniqueElements.length > 0 && (
        <div className="element-legend">
          {uniqueElements.map((item, index) => {
            // 简化显示文本 - 如果是 "CoCo3+" 这种格式就显示为 "Co3+"
            const displayText = item.element;
            
            return (
              <div key={index} className="legend-item">
                <div 
                  className="color-dot"
                  style={{ background: item.color }}
                />
                <span className="element-label">
                  {displayText}
                </span>
              </div>
            );
          })}
        </div>
      )}
    </div>
  );
};

export default CrystalViewer; 