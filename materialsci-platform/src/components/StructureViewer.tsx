import React, { useEffect, useRef, useState } from 'react';
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

interface AtomPosition {
  element: string;
  x: number;
  y: number;
  z: number;
  charge?: string; // 可选的电荷状态，例如 "Li+"、"Co3+"
}

interface StructureViewerProps {
  atoms?: AtomPosition[];
  width?: number;
  height?: number;
  showBonds?: boolean;  // 是否显示原子键
  showUnitCell?: boolean; // 是否显示晶胞
  showLabels?: boolean; // 是否显示元素标签
  backgroundColor?: string;
  rotationSpeed?: number;
  quality?: 'low' | 'medium' | 'high'; // 渲染质量
  latticeParams?: {
    a: number;
    b: number;
    c: number;
    alpha: number;
    beta: number;
    gamma: number;
  };
}

// 元素颜色映射 - 标准配色
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

// 原子半径
const atomicRadii: Record<string, number> = {
  'H': 0.46,
  'Li': 1.82,
  'C': 0.91,
  'N': 0.92,
  'O': 0.88,
  'Na': 2.27,
  'Mg': 1.73,
  'Al': 1.84,
  'Si': 1.52,
  'P': 1.48,
  'S': 1.45,
  'K': 2.75,
  'Ca': 2.31,
  'Fe': 1.56,
  'Co': 1.53,
  'Ni': 1.63,
  'Cu': 1.4,
  'Zn': 1.39,
};

// 键长阈值
const bondingDistanceThreshold: Record<string, number> = {
  'Li-O': 2.4,
  'Co-O': 2.2,
  'O-O': 2.6,
  'Si-Si': 2.8,
  'Si-O': 1.8,
  'default': 2.2
};

const StructureViewer: React.FC<StructureViewerProps> = ({ 
  atoms = [], 
  width = 400, 
  height = 300,
  showBonds = true,
  showUnitCell = true,
  showLabels = true,
  backgroundColor = '#f7f7f7',
  rotationSpeed = 0.003,
  quality = 'medium',
  latticeParams = { a: 5.4307, b: 5.4307, c: 5.4307, alpha: 90, beta: 90, gamma: 90 }
}) => {
  const containerRef = useRef<HTMLDivElement>(null);
  const [legendItems, setLegendItems] = useState<{element: string, color: string, charge?: string}[]>([]);
  
  useEffect(() => {
    if (!containerRef.current) return;

    // 清除可能存在的旧内容
    while (containerRef.current.firstChild) {
      containerRef.current.removeChild(containerRef.current.firstChild);
    }
    
    // 设置基本的渲染质量
    const segments = {
      'low': { sphere: 16, cylinder: 8 },
      'medium': { sphere: 24, cylinder: 12 },
      'high': { sphere: 32, cylinder: 16 }
    }[quality];
    
    // 初始化场景
    const scene = new THREE.Scene();
    scene.background = new THREE.Color(backgroundColor);
    
    // 基本相机设置
    const camera = new THREE.PerspectiveCamera(50, width / height, 0.1, 1000);
    camera.position.z = 15;
    camera.position.y = 5;
    camera.position.x = 5;
    
    // 基本渲染器
    const renderer = new THREE.WebGLRenderer({ 
      antialias: true,
      alpha: true
    });
    renderer.setSize(width, height);
    renderer.setPixelRatio(window.devicePixelRatio);
    containerRef.current.appendChild(renderer.domElement);
    
    // 控制器
    const controls = new OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true;
    controls.dampingFactor = 0.1;
    controls.rotateSpeed = 0.7;
    controls.enableZoom = true;
    controls.autoRotate = true;
    controls.autoRotateSpeed = rotationSpeed * 10;
    
    // 基础光照
    const ambientLight = new THREE.AmbientLight(0xffffff, 0.6);
    scene.add(ambientLight);
    
    const directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
    directionalLight.position.set(1, 1, 1);
    scene.add(directionalLight);
    
    // 计算结构中心
    let centerX = 0, centerY = 0, centerZ = 0;
    
    if (atoms.length > 0) {
      atoms.forEach(atom => {
        centerX += atom.x;
        centerY += atom.y;
        centerZ += atom.z;
      });
      centerX /= atoms.length;
      centerY /= atoms.length;
      centerZ /= atoms.length;
    }
    
    // 收集唯一元素
    const uniqueElements = new Set<string>();
    const elementCharges = new Map<string, string>();
    
    // 创建原子组
    const atomsGroup = new THREE.Group();
    scene.add(atomsGroup);
    
    // 添加原子
    atoms.forEach(atom => {
      uniqueElements.add(atom.element);
      if (atom.charge) {
        elementCharges.set(atom.element, atom.charge);
      }
      
      const color = elementColors[atom.element] || '#808080';
      const radius = atomicRadii[atom.element] || 1.0;
      
      const geometry = new THREE.SphereGeometry(radius * 0.4, segments.sphere, segments.sphere);
      const material = new THREE.MeshPhongMaterial({ 
        color: color,
        shininess: 60,
        specular: 0x444444
      });
      
      const sphere = new THREE.Mesh(geometry, material);
      
      sphere.position.set(
        atom.x - centerX, 
        atom.y - centerY, 
        atom.z - centerZ
      );
      
      atomsGroup.add(sphere);
    });
    
    // 添加键连接
    if (showBonds && atoms.length > 1) {
      const bondsGroup = new THREE.Group();
      scene.add(bondsGroup);
      
      for (let i = 0; i < atoms.length; i++) {
        for (let j = i + 1; j < atoms.length; j++) {
          const atom1 = atoms[i];
          const atom2 = atoms[j];
          
          const distance = Math.sqrt(
            Math.pow(atom1.x - atom2.x, 2) +
            Math.pow(atom1.y - atom2.y, 2) +
            Math.pow(atom1.z - atom2.z, 2)
          );
          
          const bondKey1 = `${atom1.element}-${atom2.element}`;
          const bondKey2 = `${atom2.element}-${atom1.element}`;
          const threshold = bondingDistanceThreshold[bondKey1] || 
                           bondingDistanceThreshold[bondKey2] || 
                           bondingDistanceThreshold.default;
          
          if (distance < threshold) {
            const bondMaterial = new THREE.MeshPhongMaterial({ 
              color: 0xCCCCCC,
              shininess: 30,
              transparent: true,
              opacity: 0.85
            });
            
            const direction = new THREE.Vector3(
              atom2.x - atom1.x,
              atom2.y - atom1.y,
              atom2.z - atom1.z
            );
            
            const bondGeometry = new THREE.CylinderGeometry(
              0.15, 0.15, distance, segments.cylinder, 1
            );
            bondGeometry.translate(0, distance/2, 0);
            bondGeometry.rotateX(Math.PI/2);
            
            const bond = new THREE.Mesh(bondGeometry, bondMaterial);
            
            const axis = new THREE.Vector3(0, 1, 0);
            bond.quaternion.setFromUnitVectors(
              axis, 
              direction.clone().normalize()
            );
            
            bond.position.set(
              atom1.x - centerX,
              atom1.y - centerY,
              atom1.z - centerZ
            );
            
            bondsGroup.add(bond);
          }
        }
      }
    }
    
    // 添加晶胞
    if (showUnitCell && latticeParams) {
      const unitCellGroup = new THREE.Group();
      scene.add(unitCellGroup);
      
      const { a, b, c, alpha, beta, gamma } = latticeParams;
      
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
      
      const lineMaterial = new THREE.LineBasicMaterial({ 
        color: 0x000000, 
        linewidth: 1,
        opacity: 0.7,
        transparent: true
      });
      
      // 创建晶胞顶点
      const vertices = [
        new THREE.Vector3(0, 0, 0),
        new THREE.Vector3(ax, ay, az),
        new THREE.Vector3(ax + bx, ay + by, az + bz),
        new THREE.Vector3(bx, by, bz),
        
        new THREE.Vector3(cx, cy, cz),
        new THREE.Vector3(ax + cx, ay + cy, az + cz),
        new THREE.Vector3(ax + bx + cx, ay + by + cy, az + bz + cz),
        new THREE.Vector3(bx + cx, by + cy, bz + cz)
      ];
      
      // 添加底面边线
      const bottomEdges = new THREE.BufferGeometry();
      bottomEdges.setFromPoints([
        vertices[0], vertices[1], vertices[2], vertices[3], vertices[0]
      ]);
      unitCellGroup.add(new THREE.Line(bottomEdges, lineMaterial));
      
      // 添加顶面边线
      const topEdges = new THREE.BufferGeometry();
      topEdges.setFromPoints([
        vertices[4], vertices[5], vertices[6], vertices[7], vertices[4]
      ]);
      unitCellGroup.add(new THREE.Line(topEdges, lineMaterial));
      
      // 添加侧面连接线
      for (let i = 0; i < 4; i++) {
        const sideEdge = new THREE.BufferGeometry();
        sideEdge.setFromPoints([vertices[i], vertices[i + 4]]);
        unitCellGroup.add(new THREE.Line(sideEdge, lineMaterial));
      }
      
      // 添加顶点小球
      const vertexMaterial = new THREE.MeshBasicMaterial({ color: 0x000000 });
      const vertexGeometry = new THREE.SphereGeometry(0.08, 8, 8);
      
      vertices.forEach(position => {
        const vertex = new THREE.Mesh(vertexGeometry, vertexMaterial);
        vertex.position.copy(position);
        unitCellGroup.add(vertex);
      });
      
      // 添加晶轴标签
      if (showLabels) {
        const createSimpleAxisLabel = (position: THREE.Vector3, text: string) => {
          const sprite = new THREE.Sprite(
            new THREE.SpriteMaterial({
              map: createTextTexture(text),
              transparent: true
            })
          );
          sprite.position.copy(position);
          sprite.scale.set(1, 0.5, 1);
          unitCellGroup.add(sprite);
        };
        
        createSimpleAxisLabel(new THREE.Vector3(ax * 1.1, ay, az), 'a');
        createSimpleAxisLabel(new THREE.Vector3(bx, by * 1.1, bz), 'b');
        createSimpleAxisLabel(new THREE.Vector3(cx, cy, cz * 1.1), 'c');
      }
      
      // 晶胞定位到中心
      unitCellGroup.position.set(-centerX, -centerY, -centerZ);
    }
    
    // 准备图例数据
    if (uniqueElements.size > 0) {
      const legend: {element: string, color: string, charge?: string}[] = [];
      uniqueElements.forEach(element => {
        legend.push({
          element,
          color: elementColors[element] || '#808080',
          charge: elementCharges.get(element)
        });
      });
      setLegendItems(legend);
    }
    
    // 自动调整相机视角
    const box = new THREE.Box3().setFromObject(scene);
    const size = box.getSize(new THREE.Vector3());
    const center = box.getCenter(new THREE.Vector3());
    
    const maxDim = Math.max(size.x, size.y, size.z);
    const fov = camera.fov * (Math.PI / 180);
    let cameraDistance = (maxDim / 2) / Math.tan(fov / 2);
    cameraDistance *= 1.5;
    
    camera.position.x = center.x + cameraDistance * 0.8;
    camera.position.y = center.y + cameraDistance * 0.5;
    camera.position.z = center.z + cameraDistance * 0.8;
    
    controls.target.set(center.x, center.y, center.z);
    controls.update();
    
    // 渲染循环
    const animate = () => {
      requestAnimationFrame(animate);
      controls.update();
      renderer.render(scene, camera);
    };
    
    animate();
    
    return () => {
      if (containerRef.current) {
        if (containerRef.current.contains(renderer.domElement)) {
          containerRef.current.removeChild(renderer.domElement);
        }
      }
      renderer.dispose();
    };
  }, [atoms, width, height, showBonds, showUnitCell, showLabels, backgroundColor, rotationSpeed, quality, latticeParams]);
  
  // 创建文本纹理
  const createTextTexture = (text: string): THREE.Texture => {
    const canvas = document.createElement('canvas');
    canvas.width = 64;
    canvas.height = 32;
    const context = canvas.getContext('2d');
    if (context) {
      context.fillStyle = 'white';
      context.fillRect(0, 0, 64, 32);
      context.font = 'Bold 20px Arial';
      context.fillStyle = 'black';
      context.textAlign = 'center';
      context.textBaseline = 'middle';
      context.fillText(text, 32, 16);
    }
    const texture = new THREE.CanvasTexture(canvas);
    return texture;
  };
  
  // 当没有原子数据时显示占位符
  if (atoms.length === 0) {
    return (
      <div 
        ref={containerRef}
        style={{ 
          width, 
          height, 
          background: backgroundColor, 
          display: 'flex', 
          justifyContent: 'center', 
          alignItems: 'center',
          borderRadius: '8px',
          boxShadow: 'inset 0 0 5px rgba(0,0,0,0.1)'
        }}
      >
        <span style={{ 
          color: '#888', 
          fontStyle: 'italic',
          fontSize: '14px'
        }}>等待结构数据...</span>
      </div>
    );
  }
  
  return (
    <div style={{ position: 'relative' }}>
      <div 
        ref={containerRef} 
        style={{ 
          width, 
          height, 
          borderRadius: '8px',
          overflow: 'hidden',
          boxShadow: '0 4px 8px rgba(0,0,0,0.08)'
        }} 
      />
      
      {/* 元素图例 */}
      {legendItems.length > 0 && (
        <div 
          style={{ 
            position: 'absolute',
            bottom: '10px',
            left: '50%',
            transform: 'translateX(-50%)',
            display: 'flex',
            backgroundColor: 'rgba(255, 255, 255, 0.8)',
            borderRadius: '20px',
            padding: '5px 15px',
            boxShadow: '0 2px 6px rgba(0, 0, 0, 0.15)'
          }}
        >
          {legendItems.map((item, index) => (
            <div 
              key={index}
              style={{ 
                margin: '0 8px', 
                display: 'flex', 
                alignItems: 'center',
                flexDirection: 'column',
                padding: '5px'
              }}
            >
              <div 
                style={{ 
                  width: '24px', 
                  height: '24px', 
                  borderRadius: '50%', 
                  background: item.color,
                  boxShadow: '0 2px 4px rgba(0, 0, 0, 0.2)',
                  marginBottom: '4px'
                }}
              />
              <span style={{ fontSize: '12px', fontWeight: 'bold' }}>
                {item.element}{item.charge || ''}
              </span>
            </div>
          ))}
        </div>
      )}
    </div>
  );
};

export default StructureViewer; 