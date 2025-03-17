let currentStructure = null;
let currentView = 'a'; // 默认视图方向：a轴

// 页面加载完成后初始化
document.addEventListener('DOMContentLoaded', function() {
    // 设置对话框拖动功能
    makeDialogsDraggable();
});

// 使对话框可拖动
function makeDialogsDraggable() {
    const dialogs = document.querySelectorAll('.dialog');
    
    dialogs.forEach(dialog => {
        const header = dialog.querySelector('h3');
        if (header) {
            header.style.cursor = 'move';
            header.addEventListener('mousedown', function(e) {
                const dialog = this.parentElement;
                const startX = e.clientX;
                const startY = e.clientY;
                const startLeft = parseInt(dialog.style.left) || 50;
                const startTop = parseInt(dialog.style.top) || 50;
                
                function moveDialog(e) {
                    dialog.style.left = `${startLeft + e.clientX - startX}px`;
                    dialog.style.top = `${startTop + e.clientY - startY}px`;
                    dialog.style.transform = 'none';
                }
                
                function stopMoving() {
                    document.removeEventListener('mousemove', moveDialog);
                    document.removeEventListener('mouseup', stopMoving);
                }
                
                document.addEventListener('mousemove', moveDialog);
                document.addEventListener('mouseup', stopMoving);
            });
        }
    });
}

// 文件上传处理
document.getElementById('structure-file').addEventListener('change', async (event) => {
    const file = event.target.files[0];
    if (!file) return;

    const formData = new FormData();
    formData.append('file', file);

    try {
        const response = await fetch('/upload_structure', {
            method: 'POST',
            body: formData
        });
        const data = await response.json();
        if (data.success) {
            currentStructure = data.structure;
            updateStructureViewer();
        }
    } catch (error) {
        console.error('上传文件失败:', error);
        alert('上传文件失败，请重试');
    }
});

// 创建超胞
function createSupercell() {
    const dialog = document.getElementById('supercell-dialog');
    dialog.style.display = 'block';
}

async function confirmSupercell() {
    const x = parseInt(document.getElementById('supercell-x').value) || 1;
    const y = parseInt(document.getElementById('supercell-y').value) || 1;
    const z = parseInt(document.getElementById('supercell-z').value) || 1;

    if (x < 1 || y < 1 || z < 1) {
        alert('超胞参数必须大于等于1');
        return;
    }

    try {
        const response = await fetch('/create_supercell', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({x, y, z})
        });
        const data = await response.json();
        if (data.success) {
            currentStructure = data.structure;
            updateStructureViewer();
        }
    } catch (error) {
        console.error('创建超胞失败:', error);
        alert('创建超胞失败，请重试');
    }

    document.getElementById('supercell-dialog').style.display = 'none';
}

// 添加原子
function addAtom() {
    const dialog = document.getElementById('atom-dialog');
    dialog.querySelector('h3').textContent = '添加原子';
    dialog.style.display = 'block';
    
    // 清空输入框
    document.getElementById('atom-symbol').value = '';
    document.getElementById('position-x').value = '';
    document.getElementById('position-y').value = '';
    document.getElementById('position-z').value = '';
    
    // 设置确认按钮的操作
    const confirmButton = dialog.querySelector('.dialog-buttons button:last-child');
    confirmButton.onclick = confirmAddAtom;
}

async function confirmAddAtom() {
    const symbol = document.getElementById('atom-symbol').value.trim();
    const x = parseFloat(document.getElementById('position-x').value);
    const y = parseFloat(document.getElementById('position-y').value);
    const z = parseFloat(document.getElementById('position-z').value);

    if (!symbol) {
        alert('请输入有效的元素符号');
        return;
    }

    if (isNaN(x) || isNaN(y) || isNaN(z)) {
        alert('请输入有效的坐标值');
        return;
    }

    try {
        const response = await fetch('/add_atom', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({
                symbol, 
                x, 
                y, 
                z, 
                coords_are_cartesian: false
            })
        });
        const data = await response.json();
        if (data.success) {
            currentStructure = data.structure;
            updateStructureViewer();
        } else {
            alert(data.error || '添加原子失败');
        }
    } catch (error) {
        console.error('添加原子失败:', error);
        alert('添加原子失败，请重试');
    }

    document.getElementById('atom-dialog').style.display = 'none';
}

// 选择原子
async function selectAtom() {
    if (!currentStructure || !currentStructure.species) {
        alert('请先加载或创建一个结构');
        return null;
    }
    
    // 创建原子选择对话框
    const overlay = document.createElement('div');
    overlay.className = 'overlay';
    
    const container = document.createElement('div');
    container.className = 'modal-dialog';
    
    container.innerHTML = `
        <h3 style="margin-top:0;margin-bottom:15px;">选择原子</h3>
        <div style="max-height: 400px; overflow-y: auto;">
            <table class="atom-table">
                <thead>
                    <tr>
                        <th>索引</th>
                        <th>元素</th>
                        <th>分数坐标</th>
                        <th>操作</th>
                    </tr>
                </thead>
                <tbody>
                    ${currentStructure.species.map((species, index) => `
                        <tr>
                            <td>${index}</td>
                            <td>${species}</td>
                            <td>(${currentStructure.coords[index].map(c => c.toFixed(4)).join(', ')})</td>
                            <td><button class="select-atom-btn" data-index="${index}">选择</button></td>
                        </tr>
                    `).join('')}
                </tbody>
            </table>
        </div>
        <div style="text-align: right; margin-top: 15px;">
            <button id="cancel-select-atom">取消</button>
        </div>
    `;
    
    overlay.appendChild(container);
    document.body.appendChild(overlay);
    
    return new Promise((resolve) => {
        // 为每个选择按钮添加事件监听器
        container.querySelectorAll('.select-atom-btn').forEach(btn => {
            btn.addEventListener('click', () => {
                const index = parseInt(btn.getAttribute('data-index'));
                document.body.removeChild(overlay);
                resolve(index);
            });
        });
        
        // 取消按钮
        document.getElementById('cancel-select-atom').addEventListener('click', () => {
            document.body.removeChild(overlay);
            resolve(null);
        });
    });
}

// 删除原子
async function deleteAtom() {
    // 让用户选择要删除的原子
    const selectedAtom = await selectAtom();
    if (selectedAtom === null) return;

    try {
        const response = await fetch('/delete_atom', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({atom_index: selectedAtom})
        });
        const data = await response.json();
        if (data.success) {
            currentStructure = data.structure;
            updateStructureViewer();
        } else {
            alert(data.error || '删除原子失败');
        }
    } catch (error) {
        console.error('删除原子失败:', error);
        alert('删除原子失败，请重试');
    }
}

// 修改原子位置
async function modifyPosition() {
    // 让用户选择要修改的原子
    const selectedAtom = await selectAtom();
    if (selectedAtom === null) return;
    
    const dialog = document.getElementById('atom-dialog');
    dialog.querySelector('h3').textContent = '修改原子位置';
    
    // 预填充当前坐标
    document.getElementById('atom-symbol').value = currentStructure.species[selectedAtom];
    document.getElementById('position-x').value = currentStructure.coords[selectedAtom][0];
    document.getElementById('position-y').value = currentStructure.coords[selectedAtom][1];
    document.getElementById('position-z').value = currentStructure.coords[selectedAtom][2];
    
    dialog.style.display = 'block';
    
    // 设置确认按钮的操作
    const confirmButton = dialog.querySelector('.dialog-buttons button:last-child');
    confirmButton.onclick = () => confirmModifyPosition(selectedAtom);
}

async function confirmModifyPosition(atomIndex) {
    const x = parseFloat(document.getElementById('position-x').value);
    const y = parseFloat(document.getElementById('position-y').value);
    const z = parseFloat(document.getElementById('position-z').value);

    if (isNaN(x) || isNaN(y) || isNaN(z)) {
        alert('请输入有效的坐标值');
        return;
    }

    try {
        const response = await fetch('/modify_position', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({
                atom_index: atomIndex,
                x, y, z,
                coords_are_cartesian: false
            })
        });
        const data = await response.json();
        if (data.success) {
            currentStructure = data.structure;
            updateStructureViewer();
        } else {
            alert(data.error || '修改原子位置失败');
        }
    } catch (error) {
        console.error('修改原子位置失败:', error);
        alert('修改原子位置失败，请重试');
    }

    document.getElementById('atom-dialog').style.display = 'none';
}

// 创建不可变结构
async function createImmutableStructure() {
    if (!currentStructure) {
        alert('请先加载或创建一个结构');
        return;
    }
    
    try {
        const response = await fetch('/create_immutable_structure', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify(currentStructure)
        });
        const data = await response.json();
        if (data.success) {
            alert('成功创建不可变结构！');
        } else {
            alert(data.error || '创建不可变结构失败');
        }
    } catch (error) {
        console.error('创建不可变结构失败:', error);
        alert('创建不可变结构失败，请重试');
    }
}

// 重置结构
function resetStructure() {
    currentStructure = null;
    updateStructureViewer();
    document.getElementById('structureForm').reset();
}

// 切换视图方向
function changeView(direction) {
    if (!currentStructure) return;
    
    currentView = direction;
    updateStructureImage();
    
    // 更新视图按钮状态
    const buttons = document.querySelectorAll('.view-controls button');
    buttons.forEach(btn => {
        if (btn.dataset.view === direction) {
            btn.classList.add('active');
        } else {
            btn.classList.remove('active');
        }
    });
}

// 获取结构图像
async function updateStructureImage() {
    if (!currentStructure || !currentStructure.species) return;
    
    try {
        const response = await fetch('/get_structure_image', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({
                structure: currentStructure,
                view_direction: currentView
            })
        });
        
        if (response.ok) {
            const blob = await response.blob();
            const imageUrl = URL.createObjectURL(blob);
            
            const imageContainer = document.querySelector('.structure-image');
            if (imageContainer) {
                const img = imageContainer.querySelector('img');
                if (img) {
                    img.src = imageUrl;
                } else {
                    const newImg = document.createElement('img');
                    newImg.src = imageUrl;
                    imageContainer.appendChild(newImg);
                }
            }
        } else {
            console.error('获取结构图像失败');
        }
    } catch (error) {
        console.error('获取结构图像失败:', error);
    }
}

// 更新结构显示
function updateStructureViewer() {
    const container = document.getElementById('structure-container');
    
    if (!currentStructure || !currentStructure.species) {
        container.innerHTML = '<div style="text-align:center;padding:20px;">无结构数据</div>';
        return;
    }
    
    // 显示结构信息
    const infoDiv = document.createElement('div');
    infoDiv.className = 'structure-info';
    infoDiv.innerHTML = `
        <h3>结构信息</h3>
        <p>化学式: ${currentStructure.formula || '未知'}</p>
        <p>原子数: ${currentStructure.num_sites || 0}</p>
        <p>晶格参数:</p>
        <ul>
            ${currentStructure.lattice ? `
                <li>a = ${Math.sqrt(currentStructure.lattice[0][0]**2 + currentStructure.lattice[0][1]**2 + currentStructure.lattice[0][2]**2).toFixed(4)} Å</li>
                <li>b = ${Math.sqrt(currentStructure.lattice[1][0]**2 + currentStructure.lattice[1][1]**2 + currentStructure.lattice[1][2]**2).toFixed(4)} Å</li>
                <li>c = ${Math.sqrt(currentStructure.lattice[2][0]**2 + currentStructure.lattice[2][1]**2 + currentStructure.lattice[2][2]**2).toFixed(4)} Å</li>
            ` : '<li>无晶格信息</li>'}
        </ul>
    `;
    
    // 创建视图控制按钮
    const viewControls = document.createElement('div');
    viewControls.className = 'view-controls';
    viewControls.innerHTML = `
        <button data-view="a" class="active" onclick="changeView('a')">a轴视图</button>
        <button data-view="b" onclick="changeView('b')">b轴视图</button>
        <button data-view="c" onclick="changeView('c')">c轴视图</button>
    `;
    
    // 创建图像容器
    const imageContainer = document.createElement('div');
    imageContainer.className = 'structure-image';
    
    // 清空容器并添加新元素
    container.innerHTML = '';
    container.appendChild(infoDiv);
    container.appendChild(viewControls);
    container.appendChild(imageContainer);
    
    // 获取结构图像
    updateStructureImage();
} 