// 等待DOM加载完成
document.addEventListener('DOMContentLoaded', function() {
    // 获取DOM元素
    const calculationForm = document.getElementById('calculation-form');
    const speciesSelect = document.getElementById('species');
    const initialMessage = document.getElementById('initial-message');
    const loadingIndicator = document.getElementById('loading-indicator');
    const resultsDisplay = document.getElementById('results-display');
    const errorDisplay = document.getElementById('error-display');
    const errorMessage = document.getElementById('error-message');
    
    // 结果显示元素
    const adsorptionEnergyElement = document.getElementById('adsorption-energy');
    const resultSurfaceElement = document.getElementById('result-surface');
    const resultSpeciesElement = document.getElementById('result-species');
    const resultSiteTypeElement = document.getElementById('result-site-type');
    const resultSitePositionElement = document.getElementById('result-site-position');
    
    // 加载吸附物种列表
    loadSpeciesList();
    
    // 监听表单提交事件
    calculationForm.addEventListener('submit', function(event) {
        event.preventDefault();
        
        // 验证表单
        if (!validateForm()) {
            return;
        }
        
        // 显示加载指示器
        showElement(loadingIndicator);
        hideElement(initialMessage);
        hideElement(resultsDisplay);
        hideElement(errorDisplay);
        
        // 提交表单数据
        const formData = new FormData(calculationForm);
        
        // 发送AJAX请求
        fetch('/calculate', {
            method: 'POST',
            body: formData
        })
        .then(response => response.json())
        .then(data => {
            // 隐藏加载指示器
            hideElement(loadingIndicator);
            
            if (data.success) {
                // 显示计算结果
                displayResults(data.result);
            } else {
                // 显示错误信息
                displayError(data.error);
            }
        })
        .catch(error => {
            // 隐藏加载指示器
            hideElement(loadingIndicator);
            
            // 显示错误信息
            displayError('请求失败: ' + error.message);
        });
    });
    
    // 加载吸附物种列表
    function loadSpeciesList() {
        fetch('/species')
            .then(response => response.json())
            .then(species => {
                // 清空现有选项
                speciesSelect.innerHTML = '<option value="" selected disabled>选择吸附物种</option>';
                
                // 添加物种选项
                species.forEach(item => {
                    const option = document.createElement('option');
                    option.value = item.id;
                    option.textContent = item.name;
                    speciesSelect.appendChild(option);
                });
            })
            .catch(error => {
                console.error('加载吸附物种列表失败:', error);
            });
    }
    
    // 验证表单
    function validateForm() {
        // 检查结构文件
        const structureFile = document.getElementById('structure-file').files[0];
        if (!structureFile) {
            displayError('请上传材料结构文件');
            return false;
        }
        
        // 检查吸附物种
        const species = document.getElementById('species').value;
        if (!species) {
            displayError('请选择吸附物种');
            return false;
        }
        
        // 检查Miller指数
        const millerIndex = document.getElementById('miller-index').value;
        if (!millerIndex || !/^\d+,\d+,\d+$/.test(millerIndex)) {
            displayError('请输入有效的Miller指数 (例如: 1,1,1)');
            return false;
        }
        
        return true;
    }
    
    // 显示计算结果
    function displayResults(result) {
        // 检查是否有错误
        if (result.error) {
            displayError(result.error);
            return;
        }
        
        // 填充结果数据
        adsorptionEnergyElement.textContent = result.adsorption_energy_formatted;
        resultSurfaceElement.textContent = result.surface;
        resultSpeciesElement.textContent = result.species;
        
        // 转换位点类型为中文
        let siteTypeText = '';
        switch (result.site_type) {
            case 'top':
                siteTypeText = '顶位 (Top)';
                break;
            case 'bridge':
                siteTypeText = '桥位 (Bridge)';
                break;
            case 'hollow':
                siteTypeText = '空位 (Hollow)';
                break;
            default:
                siteTypeText = result.site_type;
        }
        resultSiteTypeElement.textContent = siteTypeText;
        
        // 格式化位点坐标
        const position = result.site_position;
        resultSitePositionElement.textContent = `(${position[0].toFixed(3)}, ${position[1].toFixed(3)}, ${position[2].toFixed(3)})`;
        
        // 显示结果区域
        showElement(resultsDisplay);
        hideElement(initialMessage);
        hideElement(errorDisplay);
    }
    
    // 显示错误信息
    function displayError(message) {
        errorMessage.textContent = message;
        showElement(errorDisplay);
        hideElement(resultsDisplay);
        hideElement(initialMessage);
    }
    
    // 辅助函数：显示元素
    function showElement(element) {
        element.classList.remove('d-none');
    }
    
    // 辅助函数：隐藏元素
    function hideElement(element) {
        element.classList.add('d-none');
    }
}); 