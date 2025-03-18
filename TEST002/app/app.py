from flask import Flask, render_template, request, jsonify, send_from_directory
import os
import json
import traceback
from models.adsorption_calculator import AdsorptionCalculator
from models.material_properties import MaterialPropertiesCalculator

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = 'uploads'
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024  # 16MB限制

# 确保上传目录存在
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)

@app.route('/')
def index():
    """渲染主页"""
    return render_template('index.html')

@app.route('/material_properties')
def material_properties():
    """渲染材料属性页面"""
    return render_template('material_properties.html')

@app.route('/calculate', methods=['POST'])
def calculate_adsorption():
    """处理吸附能计算请求"""
    try:
        # 获取表单数据
        data = request.form.to_dict()
        
        # 处理上传的结构文件
        structure_file = request.files.get('structure_file')
        if structure_file:
            filename = os.path.join(app.config['UPLOAD_FOLDER'], structure_file.filename)
            structure_file.save(filename)
            data['structure_file_path'] = filename
        
        # 创建计算器实例
        calculator = AdsorptionCalculator()
        
        # 执行计算
        result = calculator.calculate(data)
        
        return jsonify({
            'success': True,
            'result': result
        })
    
    except Exception as e:
        print(f"计算吸附能时出错: {str(e)}")
        print(traceback.format_exc())
        return jsonify({
            'success': False,
            'error': str(e)
        }), 400

@app.route('/calculate_properties', methods=['POST'])
def calculate_properties():
    """处理材料属性计算请求"""
    try:
        print("接收到材料属性计算请求")
        
        # 处理上传的结构文件
        if 'structure_file' not in request.files:
            print("未找到结构文件")
            return jsonify({
                'success': False,
                'error': '未提供结构文件'
            }), 400
            
        structure_file = request.files['structure_file']
        if structure_file.filename == '':
            print("文件名为空")
            return jsonify({
                'success': False,
                'error': '文件名为空'
            }), 400
        
        # 保存上传的文件
        filename = os.path.join(app.config['UPLOAD_FOLDER'], structure_file.filename)
        structure_file.save(filename)
        print(f"文件已保存到: {filename}")
        
        # 获取要计算的属性列表
        properties_json = request.form.get('properties', '[]')
        try:
            properties = json.loads(properties_json)
            print(f"要计算的属性: {properties}")
        except json.JSONDecodeError as e:
            print(f"JSON解析错误: {str(e)}")
            return jsonify({
                'success': False,
                'error': f'属性列表格式错误: {str(e)}'
            }), 400
        
        # 创建计算器实例
        calculator = MaterialPropertiesCalculator()
        
        # 执行计算
        result = calculator.calculate({
            'structure_file_path': filename,
            'properties': properties
        })
        
        print("计算完成，返回结果")
        return jsonify({
            'success': True,
            'result': result
        })
    
    except Exception as e:
        error_trace = traceback.format_exc()
        print(f"计算材料属性时出错: {str(e)}")
        print(error_trace)
        return jsonify({
            'success': False,
            'error': str(e),
            'trace': error_trace
        }), 400

@app.route('/species')
def get_species():
    """获取可用的吸附物种列表"""
    species = [
        {'id': 'H', 'name': '氢 (H)'},
        {'id': 'O', 'name': '氧 (O)'},
        {'id': 'N', 'name': '氮 (N)'},
        {'id': 'C', 'name': '碳 (C)'},
        {'id': 'CO', 'name': '一氧化碳 (CO)'},
        {'id': 'CO2', 'name': '二氧化碳 (CO2)'},
        {'id': 'OH', 'name': '羟基 (OH)'},
        {'id': 'H2O', 'name': '水 (H2O)'}
    ]
    return jsonify(species)

@app.route('/examples/<filename>')
def download_example(filename):
    """下载示例结构文件"""
    return send_from_directory('static/examples', filename)

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000) 