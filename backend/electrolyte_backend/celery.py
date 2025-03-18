import os
from celery import Celery

# 设置Django环境
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'electrolyte_backend.settings')

app = Celery('electrolyte_backend')

# 使用字符串表示，这样解析器不会将解析对象导入到worker中
app.config_from_object('django.conf:settings', namespace='CELERY')

# 从所有已注册的Django应用配置中加载任务
app.autodiscover_tasks()

@app.task(bind=True, ignore_result=True)
def debug_task(self):
    print(f'Request: {self.request!r}') 