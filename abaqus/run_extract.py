# -*- coding: utf-8 -*-
"""
Runner script to invoke extract.py via Abaqus Python interpreter.
执行脚本：通过 Abaqus Python 解释器调用 extract.py。
"""
import os

# Get the path to the extraction script / 获取提取脚本的路径
script_path = os.path.join(os.getcwd(), 'abaqus', 'extract.py')

# Execute the script using Abaqus Python / 使用 Abaqus Python 执行脚本
os.system(f"abaqus python {script_path}")
