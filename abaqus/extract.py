# -*- coding: utf-8 -*-
"""
Extract nodal displacement data from Abaqus ODB files and save as CSV.
从 Abaqus ODB 文件中提取节点位移数据并保存为 CSV。
"""
import os
import csv
from odbAccess import openOdb

# Build path to the models directory / 构建 models 文件夹路径
models_path = os.path.join(os.getcwd(), 'abaqus', 'models')

# Find all .odb files to determine the current job index / 查找所有 .odb 文件以确定当前任务编号
files = os.listdir(models_path)
odb_files = [f for f in files if f.endswith('.odb')]
current_odb_index = len(odb_files)

# Construct file paths / 构建文件路径
odb_filename = 'job-{}.odb'.format(current_odb_index)
csv_filename = 'job-{}.csv'.format(current_odb_index)
odb_filepath = os.path.join(models_path, odb_filename)
csv_filepath = os.path.join(models_path, csv_filename)

# Open the ODB file / 打开 ODB 文件
my_odb = openOdb(odb_filepath, readOnly=False)
step = my_odb.steps['load']

# Extract time-displacement data and write to CSV / 提取时间-位移数据并写入 CSV
with open(csv_filepath, 'wb') as file:
    writer = csv.writer(file)
    writer.writerow(['Time', 'U1'])
    for frame in step.frames:
        time = frame.frameValue
        dis_field = frame.fieldOutputs['U']
        node = my_odb.rootAssembly.instances['PART-1-1'].getNodeFromLabel(label=216)
        node_displacement = dis_field.getSubset(region=node).getScalarField(componentLabel='U1')
        for dis_val in node_displacement.values:
            u1_dis = dis_val.data
            writer.writerow([time, u1_dis])
