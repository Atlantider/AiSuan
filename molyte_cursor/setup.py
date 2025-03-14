from setuptools import setup, find_packages

setup(
    name="molyte_cursor",
    version="0.1.0",
    description="重构版的Molyte分子动力学模拟工具",
    author="原作者",
    packages=find_packages(),
    package_data={
        "molyte_cursor": ["config/*/*.yaml"],
    },
    install_requires=[
        "openpyxl",
        "numpy",
        "pyyaml",
    ],
    entry_points={
        "console_scripts": [
            "molyte=molyte_cursor.src.core.main:main",
        ],
    },
) 