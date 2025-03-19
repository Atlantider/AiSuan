from rest_framework import serializers
from django.contrib.auth.models import User
from django.db import transaction
from .models import (
    Solvent, 
    Salt, 
    ElectrolyteFormulation, 
    FormulationComponent, 
    SimulationParameters,
    Calculation,
    CalculationResult,
    InputFile
)

class UserSerializer(serializers.ModelSerializer):
    class Meta:
        model = User
        fields = ['id', 'username', 'email', 'first_name', 'last_name']
        read_only_fields = ['id']

class SolventSerializer(serializers.ModelSerializer):
    class Meta:
        model = Solvent
        fields = '__all__'

class SaltSerializer(serializers.ModelSerializer):
    class Meta:
        model = Salt
        fields = '__all__'

class FormulationComponentSerializer(serializers.ModelSerializer):
    salt_name = serializers.SerializerMethodField(read_only=True)
    solvent_name = serializers.SerializerMethodField(read_only=True)
    
    class Meta:
        model = FormulationComponent
        fields = ['id', 'component_type', 'salt', 'solvent', 'concentration', 'salt_name', 'solvent_name']
    
    def get_salt_name(self, obj):
        if obj.salt:
            return obj.salt.name
        return None
    
    def get_solvent_name(self, obj):
        if obj.solvent:
            return obj.solvent.name
        return None

class SimulationParametersSerializer(serializers.ModelSerializer):
    class Meta:
        model = SimulationParameters
        exclude = ['formulation']

class ElectrolyteFormulationSerializer(serializers.ModelSerializer):
    components = FormulationComponentSerializer(many=True, read_only=True)
    parameters = SimulationParametersSerializer(read_only=True)
    user = UserSerializer(read_only=True)
    
    class Meta:
        model = ElectrolyteFormulation
        fields = '__all__'
        read_only_fields = ['id', 'created_at', 'updated_at', 'user']
    
    def create(self, validated_data):
        user = self.context['request'].user
        validated_data['user'] = user
        return super().create(validated_data)

class CalculationResultSerializer(serializers.ModelSerializer):
    class Meta:
        model = CalculationResult
        exclude = ['calculation']

class CalculationSerializer(serializers.ModelSerializer):
    formulation = ElectrolyteFormulationSerializer(read_only=True)
    formulation_id = serializers.PrimaryKeyRelatedField(
        write_only=True,
        queryset=ElectrolyteFormulation.objects.all(),
        source='formulation'
    )
    result = CalculationResultSerializer(read_only=True)
    user = UserSerializer(read_only=True)
    
    class Meta:
        model = Calculation
        fields = '__all__'
        read_only_fields = [
            'id', 'created_at', 'updated_at', 'user', 'status', 
            'task_id', 'input_file', 'output_file', 'log_file', 
            'error_message', 'started_at', 'finished_at'
        ]
    
    def create(self, validated_data):
        user = self.context['request'].user
        validated_data['user'] = user
        return super().create(validated_data)

# 表单序列化器
class CreateFormulationSerializer(serializers.Serializer):
    name = serializers.CharField(max_length=255)
    description = serializers.CharField(required=False, allow_blank=True)
    
    # 盐组件
    salts = serializers.ListField(
        child=serializers.DictField(
            child=serializers.CharField(),
            allow_empty=False
        )
    )
    
    # 溶剂组件
    solvents = serializers.ListField(
        child=serializers.DictField(
            child=serializers.CharField(),
            allow_empty=False
        )
    )
    
    # 模拟参数
    temperature = serializers.FloatField()
    pressure = serializers.FloatField(required=False, default=1.0)
    time_step = serializers.FloatField(required=False, default=1.0)
    equilibration_steps = serializers.IntegerField(required=False, default=10000)
    production_steps = serializers.IntegerField(required=False, default=50000)
    cutoff = serializers.FloatField(required=False, default=12.0)
    additional_params = serializers.JSONField(required=False, allow_null=True)
    
    def create(self, validated_data):
        user = self.context['request'].user
        
        # 检查用户是否是匿名用户，如果是，则使用ID为1的系统用户
        if user.is_anonymous:
            # 获取会话ID，用于后续识别当前用户创建的配方
            request = self.context['request']
            session_id = request.session.session_key
            if not session_id:
                request.session.create()
                session_id = request.session.session_key
                
            # 将会话ID添加到配方名称中，便于后续筛选
            name = validated_data.get('name', '')
            validated_data['name'] = f"{name}_{session_id}"
            
            try:
                # 尝试获取ID为1的用户作为系统用户
                user = User.objects.get(id=1)
            except User.DoesNotExist:
                # 尝试获取用户名为system_user的用户
                try:
                    user = User.objects.get(username="system_user")
                except User.DoesNotExist:
                    # 如果不存在，则创建一个系统用户
                    user = User.objects.create_user(
                        username="system_user",
                        email="system@aisuan.com",
                        password="defaultpassword"
                    )
                
        # 创建电解质配方
        formulation = ElectrolyteFormulation.objects.create(
            name=validated_data['name'],
            description=validated_data.get('description', ''),
            user=user
        )
        
        # 添加盐组件
        for salt_data in validated_data['salts']:
            salt_id = salt_data.get('id')
            salt = Salt.objects.get(id=salt_id)
            concentration = float(salt_data.get('concentration', 0))
            
            FormulationComponent.objects.create(
                formulation=formulation,
                component_type='salt',
                salt=salt,
                concentration=concentration
            )
        
        # 添加溶剂组件
        for solvent_data in validated_data['solvents']:
            solvent_id = solvent_data.get('id')
            solvent = Solvent.objects.get(id=solvent_id)
            concentration = float(solvent_data.get('concentration', 0))
            
            FormulationComponent.objects.create(
                formulation=formulation,
                component_type='solvent',
                solvent=solvent,
                concentration=concentration
            )
        
        # 创建模拟参数
        SimulationParameters.objects.create(
            formulation=formulation,
            temperature=validated_data['temperature'],
            pressure=validated_data.get('pressure', 1.0),
            time_step=validated_data.get('time_step', 1.0),
            equilibration_steps=validated_data.get('equilibration_steps', 10000),
            production_steps=validated_data.get('production_steps', 50000),
            cutoff=validated_data.get('cutoff', 12.0),
            additional_params=validated_data.get('additional_params')
        )
        
        return formulation 