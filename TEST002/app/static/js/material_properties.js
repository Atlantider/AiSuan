$(document).ready(function() {
    console.log("材料属性页面已加载");
    
    // 全局变量
    let currentStructureData = null;
    let currentVisualizationMode = '3dmol';
    let viewer = null;
    let threeJsRenderer = null;
    let threeJsScene = null;
    let threeJsCamera = null;
    let threeJsControls = null;
    
    // 检查3Dmol.js是否正确加载
    const is3DmolLoaded = typeof $3Dmol !== 'undefined';
    if (!is3DmolLoaded) {
        console.error("错误: 3Dmol.js 库未正确加载");
        $('#structure-debug').removeClass('d-none').addClass('alert-danger');
        $('#structure-data-debug').text("错误: 3Dmol.js 库未正确加载，将使用备用可视化方法");
        
        // 默认切换到Three.js或表格视图
        currentVisualizationMode = 'table';
        $('#view-3dmol').prop('disabled', true).removeClass('active');
        $('#view-table').addClass('active');
    } else {
        console.log("3Dmol.js 库已正确加载");
    }
    
    // 检查Three.js是否正确加载
    const isThreeJsLoaded = typeof THREE !== 'undefined';
    if (!isThreeJsLoaded) {
        console.error("错误: Three.js 库未正确加载");
        $('#view-threejs').prop('disabled', true);
        
        if (!is3DmolLoaded) {
            // 如果3Dmol和Three.js都未加载，只能使用表格视图
            $('#structure-debug').removeClass('d-none').addClass('alert-danger');
            $('#structure-data-debug').text("错误: 3Dmol.js 和 Three.js 库均未正确加载，将使用表格视图");
        }
    } else {
        console.log("Three.js 库已正确加载");
    }
    
    // 可视化方式切换
    $('#view-3dmol').on('click', function() {
        if (is3DmolLoaded) {
            switchVisualizationMode('3dmol');
        }
    });
    
    $('#view-threejs').on('click', function() {
        if (isThreeJsLoaded) {
            switchVisualizationMode('threejs');
        }
    });
    
    $('#view-table').on('click', function() {
        switchVisualizationMode('table');
    });
    
    // 切换可视化模式
    function switchVisualizationMode(mode) {
        // 更新按钮状态
        $('.btn-group button').removeClass('active');
        $(`#view-${mode}`).addClass('active');
        
        // 隐藏所有查看器
        $('#structure-viewer, #threejs-viewer, #table-viewer').addClass('d-none');
        
        // 显示选中的查看器
        currentVisualizationMode = mode;
        
        if (mode === '3dmol') {
            $('#structure-viewer').removeClass('d-none');
            if (currentStructureData && is3DmolLoaded) {
                displayStructureWith3Dmol(currentStructureData);
            }
        } else if (mode === 'threejs') {
            $('#threejs-viewer').removeClass('d-none');
            if (currentStructureData && isThreeJsLoaded) {
                displayStructureWithThreeJs(currentStructureData);
            }
        } else if (mode === 'table') {
            $('#table-viewer').removeClass('d-none');
            if (currentStructureData) {
                displayStructureAsTable(currentStructureData);
            }
        }
    }
    
    // 表单提交处理
    $('#properties-form').on('submit', function(e) {
        e.preventDefault();
        
        // 显示加载状态
        $('#loading').removeClass('d-none');
        $('#results-container').html('');
        
        // 获取选中的属性
        const selectedProperties = [];
        $('input[name="properties"]:checked').each(function() {
            selectedProperties.push($(this).val());
        });
        
        if (selectedProperties.length === 0) {
            alert('请至少选择一个要计算的属性');
            $('#loading').addClass('d-none');
            return;
        }
        
        // 创建FormData对象
        const formData = new FormData();
        const structureFile = $('#structure-file')[0].files[0];
        
        if (!structureFile) {
            alert('请上传材料结构文件');
            $('#loading').addClass('d-none');
            return;
        }
        
        // 添加文件和属性到FormData
        formData.append('structure_file', structureFile);
        formData.append('properties', JSON.stringify(selectedProperties));
        
        console.log('提交数据：', {
            文件名: structureFile.name,
            文件大小: structureFile.size,
            文件类型: structureFile.type,
            选择的属性: selectedProperties
        });
        
        // 发送AJAX请求
        $.ajax({
            url: '/calculate_properties',
            type: 'POST',
            data: formData,
            processData: false,  // 不处理数据
            contentType: false,  // 不设置内容类型
            success: function(response) {
                $('#loading').addClass('d-none');
                console.log('成功响应：', response);
                
                if (response.success) {
                    displayResults(response.result);
                    
                    // 如果包含结构数据，显示3D结构
                    if (response.result.structure_data) {
                        console.log("发现结构数据，长度:", response.result.structure_data.length);
                        currentStructureData = response.result.structure_data;
                        
                        // 根据当前选择的可视化模式显示结构
                        switchVisualizationMode(currentVisualizationMode);
                    } else {
                        console.warn("响应中没有结构数据");
                        $('#structure-debug').removeClass('d-none').addClass('alert-warning');
                        $('#structure-data-debug').text("服务器响应中没有包含结构数据");
                        
                        // 显示错误信息在结构查看器中
                        $('#structure-viewer').html(`
                            <div class="alert alert-warning">
                                <h5>无法显示结构</h5>
                                <p>服务器响应中没有包含结构数据。请确保选择了"晶体结构"属性，或者联系管理员。</p>
                            </div>
                        `);
                    }
                    
                    // 如果包含能带数据，显示能带图
                    if (response.result.bandstructure_data) {
                        displayBandstructure(response.result.bandstructure_data);
                    }
                    
                    // 如果包含态密度数据，显示态密度图
                    if (response.result.dos_data) {
                        displayDOS(response.result.dos_data);
                    }
                } else {
                    $('#results-container').html(`
                        <div class="alert alert-danger">
                            <h5>计算出错</h5>
                            <p>${response.error}</p>
                        </div>
                    `);
                }
            },
            error: function(xhr, status, error) {
                $('#loading').addClass('d-none');
                console.log('错误响应：', xhr.responseJSON || xhr.responseText);
                
                let errorMessage = error;
                if (xhr.responseJSON && xhr.responseJSON.error) {
                    errorMessage = xhr.responseJSON.error;
                }
                
                $('#results-container').html(`
                    <div class="alert alert-danger">
                        <h5>请求出错</h5>
                        <p>服务器返回错误: ${errorMessage}</p>
                    </div>
                `);
            }
        });
    });
    
    // 显示计算结果
    function displayResults(result) {
        let html = '<div class="results-summary">';
        
        // 显示基本属性
        if (result.formula) {
            html += `<div class="mb-3">
                <h5>化学式</h5>
                <p class="lead">${result.formula}</p>
            </div>`;
        }
        
        if (result.density) {
            html += `<div class="mb-3">
                <h5>密度</h5>
                <p>${result.density !== null && result.density !== undefined ? result.density.toFixed(4) + ' g/cm³' : '计算失败'}</p>
            </div>`;
        }
        
        // 显示体积信息
        if (result.volume) {
            const vol = result.volume;
            html += `<div class="mb-3">
                <h5>体积信息</h5>
                <table class="table table-sm">
                    <tr>
                        <td>晶胞体积</td>
                        <td>${vol.cell_volume !== null && vol.cell_volume !== undefined ? vol.cell_volume.toFixed(4) + ' Å³' : '计算失败'}</td>
                    </tr>
                    <tr>
                        <td>每原子体积</td>
                        <td>${vol.volume_per_atom !== null && vol.volume_per_atom !== undefined ? vol.volume_per_atom.toFixed(4) + ' Å³/atom' : '计算失败'}</td>
                    </tr>
                    <tr>
                        <td>每公式单位体积</td>
                        <td>${vol.volume_per_formula !== null && vol.volume_per_formula !== undefined ? vol.volume_per_formula.toFixed(4) + ' Å³/formula' : '计算失败'}</td>
                    </tr>
                    <tr>
                        <td>公式单位数</td>
                        <td>${vol.formula_units !== null && vol.formula_units !== undefined ? vol.formula_units : '计算失败'}</td>
                    </tr>
                </table>
            </div>`;
        }
        
        // 显示最短原子距离和配位数
        if (result.min_distance) {
            const md = result.min_distance;
            html += `<div class="mb-3">
                <h5>原子距离和配位</h5>
                <table class="table table-sm">
                    <tr>
                        <td>最短原子距离</td>
                        <td>${md.min_distance !== null && md.min_distance !== undefined ? md.min_distance.toFixed(4) + ' Å' : '计算失败'}</td>
                    </tr>
                    <tr>
                        <td>平均配位数</td>
                        <td>${md.avg_coordination_number !== null && md.avg_coordination_number !== undefined ? md.avg_coordination_number.toFixed(2) : '计算失败'}</td>
                    </tr>
                </table>`;
                
            // 显示配位细节
            if (md.coordination_details && md.coordination_details.length > 0) {
                html += `<div class="mt-2">
                    <button class="btn btn-sm btn-outline-primary" type="button" data-bs-toggle="collapse" data-bs-target="#collapseCoordination">
                        显示配位细节
                    </button>
                    <div class="collapse mt-2" id="collapseCoordination">
                        <div class="card card-body">
                            <h6>配位环境 (前${md.coordination_details.length}个原子)</h6>
                            <table class="table table-sm table-striped">
                                <thead>
                                    <tr>
                                        <th>原子</th>
                                        <th>元素</th>
                                        <th>配位数</th>
                                        <th>近邻</th>
                                    </tr>
                                </thead>
                                <tbody>`;
                
                md.coordination_details.forEach((site, index) => {
                    html += `<tr>
                        <td>${site.site_index}</td>
                        <td>${site.element}</td>
                        <td>${site.coordination_number}</td>
                        <td>`;
                    
                    if (site.neighbors && site.neighbors.length > 0) {
                        site.neighbors.forEach(neighbor => {
                            html += `${neighbor.element} (${neighbor.distance !== null && neighbor.distance !== undefined ? neighbor.distance.toFixed(2) : '?'} Å), `;
                        });
                        html = html.slice(0, -2); // 移除最后的逗号和空格
                    }
                    
                    html += `</td></tr>`;
                });
                
                html += `</tbody>
                            </table>
                        </div>
                    </div>
                </div>`;
            }
            
            html += `</div>`;
        }
        
        // 显示表面能估算
        if (result.surface_energy) {
            const se = result.surface_energy;
            
            if (se.error) {
                html += `<div class="mb-3">
                    <h5>表面能估算</h5>
                    <div class="alert alert-warning">
                        ${se.note || '无法估算表面能'}
                    </div>
                </div>`;
            } else if (se.surfaces && se.surfaces.length > 0) {
                html += `<div class="mb-3">
                    <h5>表面能估算</h5>
                    <p class="text-muted small">${se.note}</p>
                    <table class="table table-sm table-striped">
                        <thead>
                            <tr>
                                <th>晶面</th>
                                <th>表面积 (Å²)</th>
                                <th>厚度 (Å)</th>
                                <th>原子数</th>
                                <th>表面原子数</th>
                            </tr>
                        </thead>
                        <tbody>`;
                
                se.surfaces.forEach(surface => {
                    const millerStr = `(${surface.miller_index.join(',')})`;
                    html += `<tr>
                        <td>${millerStr}</td>
                        <td>${surface.surface_area !== null && surface.surface_area !== undefined ? surface.surface_area.toFixed(2) : '?'}</td>
                        <td>${surface.thickness !== null && surface.thickness !== undefined ? surface.thickness.toFixed(2) : '?'}</td>
                        <td>${surface.n_atoms}</td>
                        <td>${surface.n_surface_atoms}</td>
                    </tr>`;
                });
                
                html += `</tbody>
                    </table>
                </div>`;
            }
        }
        
        // 显示VASP输入文件
        if (result.vasp_input) {
            const vi = result.vasp_input;
            
            if (vi.error) {
                html += `<div class="mb-3">
                    <h5>VASP输入文件</h5>
                    <div class="alert alert-warning">
                        ${vi.note || '无法生成VASP输入文件'}
                    </div>
                </div>`;
            } else {
                html += `<div class="mb-3">
                    <h5>VASP输入文件</h5>
                    <p class="text-muted small">${vi.note || ''}</p>
                    
                    <div class="accordion" id="accordionVasp">
                        <!-- POSCAR -->
                        <div class="accordion-item">
                            <h2 class="accordion-header" id="headingPoscar">
                                <button class="accordion-button collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#collapsePoscar" aria-expanded="false" aria-controls="collapsePoscar">
                                    POSCAR
                                </button>
                            </h2>
                            <div id="collapsePoscar" class="accordion-collapse collapse" aria-labelledby="headingPoscar" data-bs-parent="#accordionVasp">
                                <div class="accordion-body">
                                    <pre class="border p-2 bg-light">${vi.POSCAR}</pre>
                                    <button class="btn btn-sm btn-outline-primary copy-btn" data-content="poscar">复制</button>
                                    <a href="#" class="btn btn-sm btn-outline-secondary download-btn" data-content="poscar" data-filename="POSCAR">下载</a>
                                </div>
                            </div>
                        </div>
                        
                        <!-- INCAR -->
                        <div class="accordion-item">
                            <h2 class="accordion-header" id="headingIncar">
                                <button class="accordion-button collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#collapseIncar" aria-expanded="false" aria-controls="collapseIncar">
                                    INCAR
                                </button>
                            </h2>
                            <div id="collapseIncar" class="accordion-collapse collapse" aria-labelledby="headingIncar" data-bs-parent="#accordionVasp">
                                <div class="accordion-body">
                                    <pre class="border p-2 bg-light">${vi.INCAR}</pre>
                                    <button class="btn btn-sm btn-outline-primary copy-btn" data-content="incar">复制</button>
                                    <a href="#" class="btn btn-sm btn-outline-secondary download-btn" data-content="incar" data-filename="INCAR">下载</a>
                                </div>
                            </div>
                        </div>
                        
                        <!-- KPOINTS -->
                        <div class="accordion-item">
                            <h2 class="accordion-header" id="headingKpoints">
                                <button class="accordion-button collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#collapseKpoints" aria-expanded="false" aria-controls="collapseKpoints">
                                    KPOINTS
                                </button>
                            </h2>
                            <div id="collapseKpoints" class="accordion-collapse collapse" aria-labelledby="headingKpoints" data-bs-parent="#accordionVasp">
                                <div class="accordion-body">
                                    <pre class="border p-2 bg-light">${vi.KPOINTS}</pre>
                                    <button class="btn btn-sm btn-outline-primary copy-btn" data-content="kpoints">复制</button>
                                    <a href="#" class="btn btn-sm btn-outline-secondary download-btn" data-content="kpoints" data-filename="KPOINTS">下载</a>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>`;
                
                // 存储VASP文件内容用于复制和下载
                window.vaspFiles = {
                    poscar: vi.POSCAR,
                    incar: vi.INCAR,
                    kpoints: vi.KPOINTS
                };
            }
        }
        
        if (result.lattice_parameters) {
            const lp = result.lattice_parameters;
            html += `<div class="mb-3">
                <h5>晶格参数</h5>
                <table class="table table-sm">
                    <tr>
                        <td>a</td>
                        <td>${lp.a.toFixed(4)} Å</td>
                    </tr>
                    <tr>
                        <td>b</td>
                        <td>${lp.b.toFixed(4)} Å</td>
                    </tr>
                    <tr>
                        <td>c</td>
                        <td>${lp.c.toFixed(4)} Å</td>
                    </tr>
                    <tr>
                        <td>α</td>
                        <td>${lp.alpha.toFixed(2)}°</td>
                    </tr>
                    <tr>
                        <td>β</td>
                        <td>${lp.beta.toFixed(2)}°</td>
                    </tr>
                    <tr>
                        <td>γ</td>
                        <td>${lp.gamma.toFixed(2)}°</td>
                    </tr>
                </table>
            </div>`;
        }
        
        if (result.space_group) {
            html += `<div class="mb-3">
                <h5>空间群</h5>
                <p>${result.space_group}</p>
            </div>`;
        }
        
        html += '</div>';
        $('#results-container').html(html);
        
        // 为复制按钮添加事件处理
        $('.copy-btn').on('click', function() {
            const contentType = $(this).data('content');
            const content = window.vaspFiles[contentType];
            
            if (content) {
                // 创建临时文本区域
                const textarea = document.createElement('textarea');
                textarea.value = content;
                document.body.appendChild(textarea);
                textarea.select();
                
                // 复制文本
                try {
                    document.execCommand('copy');
                    $(this).text('已复制!');
                    setTimeout(() => {
                        $(this).text('复制');
                    }, 2000);
                } catch (err) {
                    console.error('复制失败:', err);
                }
                
                // 移除临时文本区域
                document.body.removeChild(textarea);
            }
        });
        
        // 为下载按钮添加事件处理
        $('.download-btn').on('click', function(e) {
            e.preventDefault();
            
            const contentType = $(this).data('content');
            const filename = $(this).data('filename');
            const content = window.vaspFiles[contentType];
            
            if (content) {
                // 创建Blob对象
                const blob = new Blob([content], { type: 'text/plain' });
                
                // 创建下载链接
                const url = URL.createObjectURL(blob);
                const a = document.createElement('a');
                a.href = url;
                a.download = filename;
                
                // 触发下载
                document.body.appendChild(a);
                a.click();
                
                // 清理
                setTimeout(() => {
                    document.body.removeChild(a);
                    URL.revokeObjectURL(url);
                }, 0);
            }
        });
    }
    
    // 使用3Dmol.js显示3D结构
    function displayStructureWith3Dmol(structureData) {
        console.log("使用3Dmol.js显示3D结构");
        
        // 检查结构数据
        if (!structureData || structureData.trim() === '') {
            console.error("结构数据为空");
            $('#structure-debug').removeClass('d-none').addClass('alert-danger');
            $('#structure-data-debug').text("错误: 结构数据为空");
            return;
        }
        
        // 显示调试信息
        $('#structure-debug').removeClass('d-none').removeClass('alert-danger').addClass('alert-info');
        $('#structure-data-debug').text(structureData.length > 500 ? 
            structureData.substring(0, 500) + "..." : 
            structureData);
        
        // 确保结构查看器容器是空的
        $('#structure-viewer').empty();
        
        try {
            console.log("初始化3Dmol查看器");
            
            // 初始化3Dmol.js查看器
            let element = $('#structure-viewer');
            let config = {
                backgroundColor: 'white',
                antialias: true
            };
            
            // 创建查看器
            viewer = $3Dmol.createViewer(element, config);
            
            if (!viewer) {
                throw new Error("无法创建3Dmol查看器");
            }
            
            console.log("查看器创建成功，添加模型");
            
            // 加载结构数据
            let model = viewer.addModel(structureData, "xyz");
            
            if (!model) {
                throw new Error("无法加载结构模型");
            }
            
            console.log("模型添加成功，设置样式");
            
            // 设置样式
            viewer.setStyle({}, {stick:{radius: 0.15}, sphere:{scale: 0.3}});
            
            // 缩放和渲染
            console.log("缩放和渲染模型");
            viewer.zoomTo();
            viewer.render();
            
            console.log("3D结构显示完成");
        } catch (e) {
            console.error('显示3D结构时出错:', e);
            $('#structure-debug').removeClass('d-none').addClass('alert-danger');
            $('#structure-data-debug').text(`显示3D结构时出错: ${e.message}\n${e.stack}`);
            
            $('#structure-viewer').html(`
                <div class="alert alert-warning">
                    <p>显示3D结构时出错: ${e.message}</p>
                    <p>请尝试使用其他可视化方式</p>
                </div>
            `);
            
            // 自动切换到备用可视化方式
            if (isThreeJsLoaded) {
                switchVisualizationMode('threejs');
            } else {
                switchVisualizationMode('table');
            }
        }
    }
    
    // 使用Three.js显示3D结构
    function displayStructureWithThreeJs(structureData) {
        console.log("使用Three.js显示3D结构");
        
        // 检查结构数据
        if (!structureData || structureData.trim() === '') {
            console.error("结构数据为空");
            return;
        }
        
        // 解析XYZ数据
        const atoms = parseXyzData(structureData);
        if (!atoms || atoms.length === 0) {
            console.error("无法解析XYZ数据");
            $('#threejs-viewer').html(`
                <div class="alert alert-warning">
                    <p>无法解析结构数据</p>
                </div>
            `);
            return;
        }
        
        // 清空容器
        $('#threejs-viewer').empty();
        
        try {
            // 初始化Three.js
            initThreeJs();
            
            // 添加原子到场景
            addAtomsToScene(atoms);
            
            // 开始渲染
            animateThreeJs();
            
            console.log("Three.js 3D结构显示完成");
        } catch (e) {
            console.error('使用Three.js显示3D结构时出错:', e);
            $('#threejs-viewer').html(`
                <div class="alert alert-warning">
                    <p>显示3D结构时出错: ${e.message}</p>
                    <p>请尝试使用表格视图</p>
                </div>
            `);
            
            // 自动切换到表格视图
            switchVisualizationMode('table');
        }
    }
    
    // 初始化Three.js
    function initThreeJs() {
        // 创建场景
        threeJsScene = new THREE.Scene();
        threeJsScene.background = new THREE.Color(0xffffff);
        
        // 创建相机
        threeJsCamera = new THREE.PerspectiveCamera(75, $('#threejs-viewer').width() / $('#threejs-viewer').height(), 0.1, 1000);
        threeJsCamera.position.z = 5;
        
        // 创建渲染器
        threeJsRenderer = new THREE.WebGLRenderer({ antialias: true });
        threeJsRenderer.setSize($('#threejs-viewer').width(), $('#threejs-viewer').height());
        $('#threejs-viewer').append(threeJsRenderer.domElement);
        
        // 添加轨道控制
        threeJsControls = new THREE.OrbitControls(threeJsCamera, threeJsRenderer.domElement);
        threeJsControls.enableDamping = true;
        threeJsControls.dampingFactor = 0.25;
        
        // 添加光源
        const ambientLight = new THREE.AmbientLight(0x404040);
        threeJsScene.add(ambientLight);
        
        const directionalLight = new THREE.DirectionalLight(0xffffff, 0.5);
        directionalLight.position.set(1, 1, 1);
        threeJsScene.add(directionalLight);
    }
    
    // 添加原子到场景
    function addAtomsToScene(atoms) {
        // 元素颜色映射
        const elementColors = {
            'H': 0xFFFFFF,  // 白色
            'C': 0x808080,  // 灰色
            'N': 0x0000FF,  // 蓝色
            'O': 0xFF0000,  // 红色
            'F': 0x00FF00,  // 绿色
            'Cl': 0x00FF00, // 绿色
            'Br': 0xA52A2A, // 棕色
            'I': 0x800080,  // 紫色
            'He': 0xFFC0CB, // 粉色
            'Ne': 0xFFC0CB, // 粉色
            'Ar': 0xFFC0CB, // 粉色
            'Kr': 0xFFC0CB, // 粉色
            'Xe': 0xFFC0CB, // 粉色
            'P': 0xFFA500,  // 橙色
            'S': 0xFFFF00,  // 黄色
            'B': 0xFFC0CB,  // 粉色
            'Li': 0xFF0000, // 红色
            'Na': 0xFF0000, // 红色
            'K': 0xFF0000,  // 红色
            'Rb': 0xFF0000, // 红色
            'Cs': 0xFF0000, // 红色
            'Fr': 0xFF0000, // 红色
            'Be': 0x00FF00, // 绿色
            'Mg': 0x00FF00, // 绿色
            'Ca': 0x00FF00, // 绿色
            'Sr': 0x00FF00, // 绿色
            'Ba': 0x00FF00, // 绿色
            'Ra': 0x00FF00, // 绿色
            'Ti': 0x808080, // 灰色
            'Fe': 0xFFA500, // 橙色
            'Co': 0xFFA500, // 橙色
            'Ni': 0xFFA500, // 橙色
            'Cu': 0xFFA500, // 橙色
            'Zn': 0xFFA500, // 橙色
            'Ag': 0xC0C0C0, // 银色
            'Au': 0xFFD700, // 金色
            'Pt': 0xC0C0C0  // 银色
        };
        
        // 元素半径映射（单位：埃）
        const elementRadii = {
            'H': 0.25,
            'C': 0.7,
            'N': 0.65,
            'O': 0.6,
            'F': 0.5,
            'Cl': 1.0,
            'Br': 1.15,
            'I': 1.4,
            'He': 0.3,
            'Ne': 0.7,
            'Ar': 0.9,
            'Kr': 1.1,
            'Xe': 1.3,
            'P': 1.0,
            'S': 1.0,
            'B': 0.8,
            'Li': 1.45,
            'Na': 1.8,
            'K': 2.2,
            'Rb': 2.35,
            'Cs': 2.6,
            'Fr': 2.7,
            'Be': 1.05,
            'Mg': 1.5,
            'Ca': 1.8,
            'Sr': 2.0,
            'Ba': 2.15,
            'Ra': 2.15,
            'Ti': 1.4,
            'Fe': 1.4,
            'Co': 1.35,
            'Ni': 1.35,
            'Cu': 1.35,
            'Zn': 1.35,
            'Ag': 1.6,
            'Au': 1.35,
            'Pt': 1.35
        };
        
        // 计算分子中心
        let center = { x: 0, y: 0, z: 0 };
        atoms.forEach(atom => {
            center.x += atom.x;
            center.y += atom.y;
            center.z += atom.z;
        });
        center.x /= atoms.length;
        center.y /= atoms.length;
        center.z /= atoms.length;
        
        // 添加原子
        atoms.forEach(atom => {
            const element = atom.element;
            const color = elementColors[element] || 0x808080; // 默认灰色
            const radius = elementRadii[element] || 0.5; // 默认半径
            
            const geometry = new THREE.SphereGeometry(radius, 32, 32);
            const material = new THREE.MeshPhongMaterial({ color: color });
            const sphere = new THREE.Mesh(geometry, material);
            
            // 设置位置（相对于分子中心）
            sphere.position.set(
                atom.x - center.x,
                atom.y - center.y,
                atom.z - center.z
            );
            
            threeJsScene.add(sphere);
        });
        
        // 添加坐标轴辅助
        const axesHelper = new THREE.AxesHelper(5);
        threeJsScene.add(axesHelper);
    }
    
    // Three.js动画循环
    function animateThreeJs() {
        if (!threeJsRenderer) return;
        
        requestAnimationFrame(animateThreeJs);
        
        if (threeJsControls) threeJsControls.update();
        
        threeJsRenderer.render(threeJsScene, threeJsCamera);
    }
    
    // 以表格形式显示结构
    function displayStructureAsTable(structureData) {
        console.log("以表格形式显示结构");
        
        // 检查结构数据
        if (!structureData || structureData.trim() === '') {
            console.error("结构数据为空");
            $('#table-viewer').html(`
                <div class="alert alert-warning">
                    <p>结构数据为空</p>
                </div>
            `);
            return;
        }
        
        // 解析XYZ数据
        const atoms = parseXyzData(structureData);
        if (!atoms || atoms.length === 0) {
            console.error("无法解析XYZ数据");
            $('#table-viewer').html(`
                <div class="alert alert-warning">
                    <p>无法解析结构数据</p>
                </div>
            `);
            return;
        }
        
        // 清空表格
        $('#atoms-table tbody').empty();
        
        // 填充表格
        atoms.forEach((atom, index) => {
            $('#atoms-table tbody').append(`
                <tr>
                    <td>${index + 1}</td>
                    <td>${atom.element}</td>
                    <td>${atom.x.toFixed(4)}</td>
                    <td>${atom.y.toFixed(4)}</td>
                    <td>${atom.z.toFixed(4)}</td>
                </tr>
            `);
        });
    }
    
    // 解析XYZ格式数据
    function parseXyzData(xyzData) {
        try {
            const lines = xyzData.trim().split('\n');
            
            // 第一行是原子数量
            const numAtoms = parseInt(lines[0].trim());
            
            // 检查数据格式
            if (isNaN(numAtoms) || lines.length < numAtoms + 2) {
                console.error("XYZ数据格式不正确");
                return null;
            }
            
            // 解析原子数据（从第三行开始）
            const atoms = [];
            for (let i = 2; i < numAtoms + 2; i++) {
                const parts = lines[i].trim().split(/\s+/);
                if (parts.length >= 4) {
                    atoms.push({
                        element: parts[0],
                        x: parseFloat(parts[1]),
                        y: parseFloat(parts[2]),
                        z: parseFloat(parts[3])
                    });
                }
            }
            
            return atoms;
        } catch (e) {
            console.error("解析XYZ数据时出错:", e);
            return null;
        }
    }
    
    // 显示3D结构（主函数，根据当前模式选择合适的显示方法）
    function displayStructure(structureData) {
        console.log("开始显示3D结构，当前模式:", currentVisualizationMode);
        
        // 保存结构数据
        currentStructureData = structureData;
        
        // 根据当前模式选择显示方法
        if (currentVisualizationMode === '3dmol' && is3DmolLoaded) {
            displayStructureWith3Dmol(structureData);
        } else if (currentVisualizationMode === 'threejs' && isThreeJsLoaded) {
            displayStructureWithThreeJs(structureData);
        } else {
            displayStructureAsTable(structureData);
        }
    }
    
    // 显示能带结构
    function displayBandstructure(bandData) {
        try {
            // 处理能带数据格式
            const traces = [];
            
            // 如果是多条能带
            if (Array.isArray(bandData.energies[0])) {
                for (let i = 0; i < bandData.energies.length; i++) {
                    traces.push({
                        x: bandData.kpoints,
                        y: bandData.energies[i],
                        mode: 'lines',
                        line: {color: 'blue', width: 2},
                        name: `能带 ${i+1}`
                    });
                }
            } else {
                // 单条能带
                traces.push({
                    x: bandData.kpoints,
                    y: bandData.energies,
                    mode: 'lines',
                    line: {color: 'blue', width: 2},
                    name: '能带'
                });
            }
            
            const layout = {
                title: '能带结构',
                xaxis: {
                    title: 'k路径',
                    tickvals: bandData.tick_vals,
                    ticktext: bandData.tick_labels
                },
                yaxis: {
                    title: '能量 (eV)'
                },
                showlegend: false,
                margin: {
                    l: 50,
                    r: 30,
                    t: 50,
                    b: 50
                },
                autosize: true
            };
            
            const config = {
                responsive: true
            };
            
            Plotly.newPlot('bandstructure-plot', traces, layout, config);
        } catch (e) {
            console.error('显示能带结构时出错:', e);
            $('#bandstructure-plot').html(`
                <div class="alert alert-warning">
                    <p>显示能带结构时出错: ${e.message}</p>
                </div>
            `);
        }
    }
    
    // 显示态密度
    function displayDOS(dosData) {
        try {
            const data = [{
                x: dosData.energies,
                y: dosData.densities,
                mode: 'lines',
                line: {color: 'red', width: 2},
                name: '态密度'
            }];
            
            const layout = {
                title: '态密度',
                xaxis: {
                    title: '能量 (eV)'
                },
                yaxis: {
                    title: '态密度'
                },
                showlegend: false,
                margin: {
                    l: 50,
                    r: 30,
                    t: 50,
                    b: 50
                },
                autosize: true
            };
            
            const config = {
                responsive: true
            };
            
            Plotly.newPlot('dos-plot', data, layout, config);
        } catch (e) {
            console.error('显示态密度时出错:', e);
            $('#dos-plot').html(`
                <div class="alert alert-warning">
                    <p>显示态密度时出错: ${e.message}</p>
                </div>
            `);
        }
    }
    
    // 窗口大小改变时调整图表大小
    $(window).on('resize', function() {
        if (viewer) {
            viewer.resize();
        }
        
        if (threeJsRenderer) {
            threeJsCamera.aspect = $('#threejs-viewer').width() / $('#threejs-viewer').height();
            threeJsCamera.updateProjectionMatrix();
            threeJsRenderer.setSize($('#threejs-viewer').width(), $('#threejs-viewer').height());
        }
        
        // 重新绘制Plotly图表 - 修复：只对已经初始化的图表调用relayout
        $('.plot-container').each(function() {
            const plotId = $(this).attr('id');
            // 检查元素是否存在且已经被Plotly初始化
            if (plotId && document.getElementById(plotId) && 
                $(this).html() !== '' && 
                !$(this).find('.alert').length && 
                document.getElementById(plotId).data && 
                document.getElementById(plotId).data.length > 0) {
                try {
                    Plotly.relayout(plotId, {
                        autosize: true
                    });
                } catch (e) {
                    console.warn(`无法重新布局图表 ${plotId}:`, e);
                }
            }
        });
    });
}); 