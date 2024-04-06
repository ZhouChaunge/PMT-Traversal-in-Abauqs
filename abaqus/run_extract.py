# -*- coding: utf-8 -*-
import os

script_path = os.path.join(os.getcwd(), 'abaqus', 'extract.py')  # 获取当前脚本的目录
os.system(f"abaqus python {script_path}")  # 执行Python脚本文件
