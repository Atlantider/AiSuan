from django.core.management.base import BaseCommand
from django.core.management import call_command
import os
from django.conf import settings

class Command(BaseCommand):
    help = '加载电解质应用的初始数据'

    def handle(self, *args, **options):
        self.stdout.write(self.style.SUCCESS('开始加载初始数据...'))
        
        # 加载初始数据
        fixture_path = os.path.join('electrolyte', 'fixtures', 'initial_data.json')
        call_command('loaddata', fixture_path, verbosity=1)
        
        self.stdout.write(self.style.SUCCESS('初始数据加载完成!'))
        self.stdout.write(self.style.SUCCESS('已添加示例溶剂和盐数据。')) 