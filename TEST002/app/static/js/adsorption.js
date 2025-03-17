document.addEventListener('DOMContentLoaded', function() {
    // 获取DOM元素
    const structureFileInput = document.getElementById('structure-file');
    const millerIndexInput = document.getElementById('miller-index');
    const slabThicknessInput = document.getElementById('slab-thickness');
    const vacuumThicknessInput = document.getElementById('vacuum-thickness');
    const adsorbateSpeciesSelect = document.getElementById('adsorbate-species');
    const adsorptionSiteSelect = document.getElementById('adsorption-site');
    const calculateBtn = document.getElementById('calculate-btn');
    
    const viewABtn = document.getElementById('view-a');
    const viewBBtn = document.getElementById('view-b');
    const viewCBtn = document.getElementById('view-c');
    
    const tabSitesBtn = document.getElementById('tab-sites');
    const tabStructureBtn = document.getElementById('tab-structure');
    
    const sitesView = document.getElementById('sites-view');
    const structureView = document.getElementById('structure-view');
    
    const sitesImage = document.getElementById('sites-image');
    const structureImage = document.getElementById('structure-image');
    
    const resultsContainer = document.getElementById('results-container');
    
    // 当前视图方向
    let currentViewDirection = 'c';
    
    // 上传的文件路径
    let uploadedFilePath = null;
    
    // 计算结果
    let calculationResult = null;
    
    // 初始化加载状态
    const sitesLoadingOverlay = sitesView.querySelector('.loading-overlay');
    const structureLoadingOverlay = structureView.querySelector('.loading-overlay');
    sitesLoadingOverlay.style.display = 'none';
    structureLoadingOverlay.style.display = 'none';
    
    // 处理文件上传
    structureFileInput.addEventListener('change', function(event) {
        const file = event.target.files[0];
        if (!file) return;
        
        // 创建FormData对象
        const formData = new FormData();
        formData.append('structure_file', file);
        
        // 显示加载状态
        showMessage('正在上传文件...');
        
        // 发送请求
        fetch('/upload_structure', {
            method: 'POST',
            body: formData
        })
        .then(response => response.json())
        .then(data => {
            if (data.error) {
                showError(data.error);
                return;
            }
            
            // 保存上传的文件路径
            uploadedFilePath = data.file_path;
            showMessage('文件上传成功: ' + file.name);
        })
        .catch(error => {
            showError('上传文件时出错: ' + error.message);
        });
    });
    
    // 处理计算按钮点击
    calculateBtn.addEventListener('click', function() {
        // 检查是否上传了文件
        if (!uploadedFilePath) {
            showError('请先上传结构文件');
            return;
        }
        
        // 获取输入参数
        const millerIndex = millerIndexInput.value;
        const slabThickness = slabThicknessInput.value;
        const vacuumThickness = vacuumThicknessInput.value;
        const species = adsorbateSpeciesSelect.value;
        const adsorptionSite = adsorptionSiteSelect.value;
        
        // 验证输入
        if (!millerIndex || !slabThickness || !vacuumThickness) {
            showError('请填写所有必填字段');
            return;
        }
        
        // 准备请求数据
        const requestData = {
            structure_file_path: uploadedFilePath,
            miller_index: millerIndex,
            slab_thickness: slabThickness,
            vacuum_thickness: vacuumThickness,
            species: species,
            adsorption_site: adsorptionSite
        };
        
        // 显示加载状态
        showMessage('正在计算吸附能...');
        showMessage('请求数据: ' + JSON.stringify(requestData));
        sitesLoadingOverlay.style.display = 'flex';
        structureLoadingOverlay.style.display = 'flex';
        
        // 发送请求
        fetch('/calculate_adsorption', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify(requestData)
        })
        .then(response => {
            if (!response.ok) {
                return response.json().then(data => {
                    throw new Error('计算失败: ' + (data.error || response.statusText));
                });
            }
            return response.json();
        })
        .then(data => {
            if (data.error) {
                showError(data.error);
                sitesLoadingOverlay.style.display = 'none';
                structureLoadingOverlay.style.display = 'none';
                return;
            }
            
            // 保存计算结果
            calculationResult = data;
            
            // 显示计算结果
            displayResults(data);
            
            // 加载图像
            loadImages();
        })
        .catch(error => {
            showError('计算吸附能时出错: ' + error.message);
            sitesLoadingOverlay.style.display = 'none';
            structureLoadingOverlay.style.display = 'none';
        });
    });
    
    // 处理视图方向按钮点击
    viewABtn.addEventListener('click', function() {
        setViewDirection('a');
    });
    
    viewBBtn.addEventListener('click', function() {
        setViewDirection('b');
    });
    
    viewCBtn.addEventListener('click', function() {
        setViewDirection('c');
    });
    
    // 处理标签页按钮点击
    tabSitesBtn.addEventListener('click', function() {
        setActiveTab('sites');
    });
    
    tabStructureBtn.addEventListener('click', function() {
        setActiveTab('structure');
    });
    
    // 设置视图方向
    function setViewDirection(direction) {
        // 更新当前视图方向
        currentViewDirection = direction;
        
        // 更新按钮状态
        viewABtn.classList.remove('active');
        viewBBtn.classList.remove('active');
        viewCBtn.classList.remove('active');
        
        if (direction === 'a') {
            viewABtn.classList.add('active');
        } else if (direction === 'b') {
            viewBBtn.classList.add('active');
        } else {
            viewCBtn.classList.add('active');
        }
        
        // 如果有计算结果，重新加载图像
        if (calculationResult) {
            loadImages();
        }
    }
    
    // 设置活动标签页
    function setActiveTab(tab) {
        // 更新按钮状态
        tabSitesBtn.classList.remove('active');
        tabStructureBtn.classList.remove('active');
        
        // 更新标签页显示
        sitesView.classList.remove('active');
        structureView.classList.remove('active');
        
        if (tab === 'sites') {
            tabSitesBtn.classList.add('active');
            sitesView.classList.add('active');
        } else {
            tabStructureBtn.classList.add('active');
            structureView.classList.add('active');
        }
    }
    
    // 加载图像
    function loadImages() {
        // 显示加载状态
        sitesLoadingOverlay.style.display = 'flex';
        structureLoadingOverlay.style.display = 'flex';
        
        // 加载吸附位点图像
        fetch('/get_adsorption_sites_image', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({ view_direction: currentViewDirection })
        })
        .then(response => {
            if (!response.ok) {
                throw new Error('加载吸附位点图像失败');
            }
            return response.blob();
        })
        .then(blob => {
            const imageUrl = URL.createObjectURL(blob);
            sitesImage.src = imageUrl;
            sitesImage.onload = function() {
                sitesLoadingOverlay.style.display = 'none';
            };
        })
        .catch(error => {
            showError('加载吸附位点图像时出错: ' + error.message);
            sitesLoadingOverlay.style.display = 'none';
        });
        
        // 加载吸附结构图像
        fetch('/get_adsorption_structure_image', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({ view_direction: currentViewDirection })
        })
        .then(response => {
            if (!response.ok) {
                throw new Error('加载吸附结构图像失败');
            }
            return response.blob();
        })
        .then(blob => {
            const imageUrl = URL.createObjectURL(blob);
            structureImage.src = imageUrl;
            structureImage.onload = function() {
                structureLoadingOverlay.style.display = 'none';
            };
        })
        .catch(error => {
            showError('加载吸附结构图像时出错: ' + error.message);
            structureLoadingOverlay.style.display = 'none';
        });
    }
    
    // 显示计算结果
    function displayResults(data) {
        // 创建结果HTML
        let html = `
            <div class="result-item">
                <span class="result-label">表面:</span>
                <span class="result-value">${data.surface}</span>
            </div>
            <div class="result-item">
                <span class="result-label">吸附物种:</span>
                <span class="result-value">${data.species}</span>
            </div>
            <div class="result-item">
                <span class="result-label">吸附位点类型:</span>
                <span class="result-value">${data.site_type}</span>
            </div>
            <div class="result-item">
                <span class="result-label">吸附位点坐标:</span>
                <span class="result-value">[${data.site_position.join(', ')}]</span>
            </div>
            <div class="result-item highlight">
                <span class="result-label">吸附能:</span>
                <span class="result-value">${data.adsorption_energy} eV</span>
            </div>
        `;
        
        // 更新结果容器
        resultsContainer.innerHTML = html;
    }
    
    // 显示错误消息
    function showError(message) {
        resultsContainer.innerHTML = `<p class="error">${message}</p>`;
    }
    
    // 显示一般消息
    function showMessage(message) {
        resultsContainer.innerHTML = `<p>${message}</p>`;
    }
}); 