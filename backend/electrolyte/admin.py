from django.contrib import admin
from .models import (
    Solvent, 
    Salt, 
    ElectrolyteFormulation, 
    FormulationComponent, 
    SimulationParameters,
    Calculation,
    CalculationResult
)

class FormulationComponentInline(admin.TabularInline):
    model = FormulationComponent
    extra = 1

class SimulationParametersInline(admin.StackedInline):
    model = SimulationParameters
    can_delete = False
    max_num = 1

class CalculationResultInline(admin.StackedInline):
    model = CalculationResult
    can_delete = False
    max_num = 1
    readonly_fields = [
        'ionic_conductivity', 'diffusion_coefficients', 
        'radial_distribution', 'density', 'viscosity', 
        'additional_results'
    ]

@admin.register(Solvent)
class SolventAdmin(admin.ModelAdmin):
    list_display = ['name', 'smile', 'created_at', 'updated_at']
    search_fields = ['name', 'smile', 'description']
    list_filter = ['created_at', 'updated_at']

@admin.register(Salt)
class SaltAdmin(admin.ModelAdmin):
    list_display = ['name', 'cation', 'anion', 'created_at', 'updated_at']
    search_fields = ['name', 'cation', 'anion', 'description']
    list_filter = ['created_at', 'updated_at']

@admin.register(ElectrolyteFormulation)
class ElectrolyteFormulationAdmin(admin.ModelAdmin):
    list_display = ['name', 'user', 'created_at', 'updated_at']
    search_fields = ['name', 'description', 'user__username']
    list_filter = ['created_at', 'updated_at', 'user']
    inlines = [FormulationComponentInline, SimulationParametersInline]

@admin.register(Calculation)
class CalculationAdmin(admin.ModelAdmin):
    list_display = ['name', 'user', 'formulation', 'status', 'created_at', 'started_at', 'finished_at']
    search_fields = ['name', 'description', 'user__username', 'formulation__name']
    list_filter = ['status', 'created_at', 'started_at', 'finished_at', 'user']
    readonly_fields = ['task_id', 'input_file', 'output_file', 'log_file', 'error_message', 'started_at', 'finished_at']
    inlines = [CalculationResultInline]
