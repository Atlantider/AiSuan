"""
用户管理模块

负责用户认证、权限管理和用户数据管理。
"""

import os
import json
import hashlib
import logging
from typing import Dict, List, Optional, Any, Union
from datetime import datetime
import shutil

class UserManager:
    """
    用户管理类
    
    处理用户认证、权限管理和用户数据操作
    """
    
    # 用户角色定义
    ROLE_ADMIN = "admin"       # 管理员
    ROLE_RESEARCHER = "researcher"  # 研究人员
    ROLE_VIEWER = "viewer"     # 只读用户
    
    # 权限定义
    PERMISSION_READ = "read"       # 读权限
    PERMISSION_WRITE = "write"     # 写权限
    PERMISSION_EXECUTE = "execute" # 执行权限
    PERMISSION_ADMIN = "admin"     # 管理权限
    
    def __init__(self, workspace_root: Optional[str] = None):
        """
        初始化用户管理器
        
        Args:
            workspace_root: 工作空间根目录，默认为AiSuan/workspace
        """
        # 设置工作空间根目录
        if workspace_root is None:
            # 获取AiSuan根目录
            aisuan_root = os.environ.get('AISUAN_ROOT', '/public/home/xiaoji/AiSuan')
            workspace_root = os.path.join(aisuan_root, 'workspace')
        
        self.workspace_root = workspace_root
        
        # 创建系统配置目录
        self.system_dir = os.path.join(workspace_root, 'system')
        self.config_dir = os.path.join(self.system_dir, 'configs')
        os.makedirs(self.config_dir, exist_ok=True)
        
        # 用户数据文件
        self.users_file = os.path.join(self.config_dir, 'users.json')
        
        # 初始化日志
        self.logger = logging.getLogger(__name__)
        
        # 加载用户数据
        self.users = self._load_users()
        
        # 如果用户文件不存在或为空，初始化默认用户
        if not self.users:
            self._initialize_default_users()
    
    def _load_users(self) -> Dict[str, Any]:
        """
        加载用户数据
        
        Returns:
            Dict: 用户数据字典
        """
        if os.path.exists(self.users_file):
            try:
                with open(self.users_file, 'r', encoding='utf-8') as f:
                    return json.load(f)
            except Exception as e:
                self.logger.error(f"加载用户数据失败: {str(e)}")
        
        return {}
    
    def _save_users(self) -> None:
        """保存用户数据到文件"""
        try:
            with open(self.users_file, 'w', encoding='utf-8') as f:
                json.dump(self.users, f, indent=2, ensure_ascii=False)
            
            # 设置文件权限
            os.chmod(self.users_file, 0o600)  # 只有所有者可读写
            
            self.logger.info(f"用户数据已保存到: {self.users_file}")
        except Exception as e:
            self.logger.error(f"保存用户数据失败: {str(e)}")
    
    def _initialize_default_users(self) -> None:
        """初始化默认用户数据"""
        # 创建admin用户
        admin_user = {
            "email": "admin@aisuan.com",
            "username": "admin",
            "password_hash": self._hash_password("admin123"),
            "role": self.ROLE_ADMIN,
            "permissions": [self.PERMISSION_READ, self.PERMISSION_WRITE, self.PERMISSION_EXECUTE, self.PERMISSION_ADMIN],
            "created_at": datetime.now().isoformat(),
            "last_login": None,
            "projects": []
        }
        
        self.users = {"admin@aisuan.com": admin_user}
        self._save_users()
        
        # 创建admin用户目录
        admin_dir = os.path.join(self.workspace_root, 'users', 'admin')
        os.makedirs(admin_dir, exist_ok=True)
        
        # 确保admin用户目录下有必要的子目录
        for subdir in ['projects', 'molecules', 'shared']:
            os.makedirs(os.path.join(admin_dir, subdir), exist_ok=True)
        
        self.logger.info("已初始化默认用户数据和目录")
    
    def _hash_password(self, password: str) -> str:
        """
        对密码进行哈希处理
        
        Args:
            password: 原始密码
            
        Returns:
            str: 密码哈希值
        """
        return hashlib.sha256(password.encode()).hexdigest()
    
    def authenticate(self, email: str, password: str) -> Optional[Dict[str, Any]]:
        """
        验证用户身份
        
        Args:
            email: 用户邮箱
            password: 用户密码
            
        Returns:
            Optional[Dict]: 用户信息，如果验证失败则返回None
        """
        if email in self.users:
            user = self.users[email]
            if user["password_hash"] == self._hash_password(password):
                # 更新最后登录时间
                user["last_login"] = datetime.now().isoformat()
                self._save_users()
                
                # 返回用户信息（不包含密码哈希）
                user_info = user.copy()
                user_info.pop("password_hash", None)
                return user_info
        
        return None
    
    def create_user(self, email: str, username: str, password: str, role: str = "researcher") -> bool:
        """
        创建新用户
        
        Args:
            email: 用户邮箱（用作唯一标识）
            username: 用户名
            password: 用户密码
            role: 用户角色，默认为researcher
            
        Returns:
            bool: 创建成功返回True，否则返回False
        """
        if email in self.users:
            self.logger.warning(f"用户已存在: {email}")
            return False
        
        # 检查角色有效性
        if role not in [self.ROLE_ADMIN, self.ROLE_RESEARCHER, self.ROLE_VIEWER]:
            self.logger.warning(f"无效的角色: {role}")
            return False
        
        # 根据角色设置权限
        permissions = [self.PERMISSION_READ]
        if role == self.ROLE_RESEARCHER:
            permissions.extend([self.PERMISSION_WRITE, self.PERMISSION_EXECUTE])
        elif role == self.ROLE_ADMIN:
            permissions.extend([self.PERMISSION_WRITE, self.PERMISSION_EXECUTE, self.PERMISSION_ADMIN])
        
        # 创建用户记录
        user = {
            "email": email,
            "username": username,
            "password_hash": self._hash_password(password),
            "role": role,
            "permissions": permissions,
            "created_at": datetime.now().isoformat(),
            "last_login": None,
            "projects": []
        }
        
        self.users[email] = user
        self._save_users()
        
        # 创建用户目录
        user_dir = os.path.join(self.workspace_root, 'users', username)
        os.makedirs(user_dir, exist_ok=True)
        
        # 创建用户子目录
        for subdir in ['projects', 'molecules', 'shared']:
            os.makedirs(os.path.join(user_dir, subdir), exist_ok=True)
        
        self.logger.info(f"成功创建用户: {email} ({username})")
        return True
    
    def update_user(self, email: str, data: Dict[str, Any]) -> bool:
        """
        更新用户信息
        
        Args:
            email: 用户邮箱
            data: 要更新的用户数据
            
        Returns:
            bool: 更新成功返回True，否则返回False
        """
        if email not in self.users:
            self.logger.warning(f"用户不存在: {email}")
            return False
        
        user = self.users[email]
        
        # 更新用户数据
        for key, value in data.items():
            if key == 'password':
                user['password_hash'] = self._hash_password(value)
            elif key != 'password_hash' and key in user:  # 不直接设置password_hash
                user[key] = value
        
        self._save_users()
        self.logger.info(f"成功更新用户信息: {email}")
        return True
    
    def delete_user(self, email: str) -> bool:
        """
        删除用户
        
        Args:
            email: 用户邮箱
            
        Returns:
            bool: 删除成功返回True，否则返回False
        """
        if email not in self.users:
            self.logger.warning(f"用户不存在: {email}")
            return False
        
        # 不允许删除最后一个管理员账户
        if self.users[email]['role'] == self.ROLE_ADMIN:
            admin_count = sum(1 for u in self.users.values() if u['role'] == self.ROLE_ADMIN)
            if admin_count <= 1:
                self.logger.warning("不能删除唯一的管理员账户")
                return False
        
        # 获取用户名，用于删除用户目录
        username = self.users[email]['username']
        
        # 删除用户记录
        del self.users[email]
        self._save_users()
        
        # 删除或归档用户目录
        user_dir = os.path.join(self.workspace_root, 'users', username)
        if os.path.exists(user_dir):
            # 可以选择删除或归档
            # 这里选择重命名为归档
            archive_dir = os.path.join(self.workspace_root, 'users', f"{username}_archived_{int(datetime.now().timestamp())}")
            try:
                os.rename(user_dir, archive_dir)
                self.logger.info(f"用户目录已归档: {archive_dir}")
            except Exception as e:
                self.logger.error(f"归档用户目录失败: {str(e)}")
        
        self.logger.info(f"成功删除用户: {email}")
        return True
    
    def check_permission(self, user_email: str, permission: str) -> bool:
        """
        检查用户是否拥有指定权限
        
        Args:
            user_email: 用户邮箱
            permission: 要检查的权限
            
        Returns:
            bool: 用户拥有权限返回True，否则返回False
        """
        if user_email not in self.users:
            return False
        
        user = self.users[user_email]
        return permission in user['permissions']
    
    def get_user_info(self, email: str) -> Optional[Dict[str, Any]]:
        """
        获取用户信息
        
        Args:
            email: 用户邮箱
            
        Returns:
            Optional[Dict]: 用户信息（不包含密码哈希），用户不存在则返回None
        """
        if email in self.users:
            user_info = self.users[email].copy()
            user_info.pop("password_hash", None)
            return user_info
        
        return None
    
    def list_users(self) -> List[Dict[str, Any]]:
        """
        获取所有用户列表
        
        Returns:
            List[Dict]: 用户信息列表（不包含密码哈希）
        """
        return [
            {k: v for k, v in user.items() if k != "password_hash"}
            for user in self.users.values()
        ]
    
    def create_project(self, user_email: str, project_name: str) -> Optional[str]:
        """
        为用户创建项目
        
        Args:
            user_email: 用户邮箱
            project_name: 项目名称
            
        Returns:
            Optional[str]: 项目ID，如果创建失败则返回None
        """
        if user_email not in self.users:
            self.logger.warning(f"用户不存在: {user_email}")
            return None
        
        user = self.users[user_email]
        
        # 检查写权限
        if self.PERMISSION_WRITE not in user['permissions']:
            self.logger.warning(f"用户 {user_email} 没有创建项目的权限")
            return None
        
        # 生成项目ID
        timestamp = datetime.now().strftime('%Y%m%d%H%M%S')
        project_id = f"proj_{timestamp}_{project_name.replace(' ', '_')}"
        
        # 创建项目目录
        user_dir = os.path.join(self.workspace_root, 'users', user['username'])
        project_dir = os.path.join(user_dir, 'projects', project_id)
        os.makedirs(project_dir, exist_ok=True)
        
        # 创建项目子目录
        for subdir in ['raw', 'molecules', 'systems', 'simulations', 'analysis', 'reports', 'logs']:
            os.makedirs(os.path.join(project_dir, subdir), exist_ok=True)
        
        # 更新用户项目列表
        project_info = {
            "id": project_id,
            "name": project_name,
            "created_at": datetime.now().isoformat(),
            "path": project_dir
        }
        
        user['projects'].append(project_info)
        self._save_users()
        
        self.logger.info(f"为用户 {user_email} 创建项目: {project_id}")
        return project_id
    
    def get_user_projects(self, user_email: str) -> List[Dict[str, Any]]:
        """
        获取用户的项目列表
        
        Args:
            user_email: 用户邮箱
            
        Returns:
            List[Dict]: 项目信息列表
        """
        if user_email in self.users:
            return self.users[user_email].get('projects', [])
        
        return []
    
    def get_project_path(self, user_email: str, project_id: str) -> Optional[str]:
        """
        获取项目路径
        
        Args:
            user_email: 用户邮箱
            project_id: 项目ID
            
        Returns:
            Optional[str]: 项目路径，如果项目不存在则返回None
        """
        projects = self.get_user_projects(user_email)
        for project in projects:
            if project['id'] == project_id:
                return project['path']
        
        return None
    
    def set_directory_permissions(self, directory: str, is_admin: bool = False) -> None:
        """
        设置目录权限
        
        Args:
            directory: 目录路径
            is_admin: 是否为管理员目录
        """
        if not os.path.exists(directory):
            os.makedirs(directory, exist_ok=True)
        
        if os.name == 'posix':  # 仅在Unix/Linux系统上设置
            if is_admin:
                # 管理员目录: 只有所有者有完全权限 (700)
                os.chmod(directory, 0o700)
            else:
                # 普通目录: 所有者有完全权限，组用户可读可执行 (750)
                os.chmod(directory, 0o750)

# 创建命令行工具函数
def create_admin_user(email: str, username: str, password: str) -> bool:
    """
    创建管理员用户
    
    Args:
        email: 管理员邮箱
        username: 管理员用户名
        password: 管理员密码
        
    Returns:
        bool: 创建成功返回True，否则返回False
    """
    try:
        # 初始化用户管理器
        user_manager = UserManager()
        
        # 检查用户是否已存在
        if email in user_manager.users:
            user = user_manager.users[email]
            if user['role'] == UserManager.ROLE_ADMIN:
                print(f"管理员 {email} 已存在")
                return True
            else:
                # 将现有用户升级为管理员
                update_data = {
                    "role": UserManager.ROLE_ADMIN,
                    "permissions": [
                        UserManager.PERMISSION_READ,
                        UserManager.PERMISSION_WRITE,
                        UserManager.PERMISSION_EXECUTE,
                        UserManager.PERMISSION_ADMIN
                    ]
                }
                if user_manager.update_user(email, update_data):
                    print(f"已将用户 {email} 升级为管理员")
                    return True
                else:
                    print(f"升级用户 {email} 为管理员失败")
                    return False
        
        # 创建新管理员用户
        if user_manager.create_user(email, username, password, role=UserManager.ROLE_ADMIN):
            print(f"成功创建管理员用户: {email} ({username})")
            return True
        else:
            print(f"创建管理员用户失败: {email}")
            return False
    
    except Exception as e:
        print(f"创建管理员用户时出错: {str(e)}")
        return False

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) == 4:
        email = sys.argv[1]
        username = sys.argv[2]
        password = sys.argv[3]
        
        create_admin_user(email, username, password)
    else:
        print("用法: python user_manager.py <email> <username> <password>") 