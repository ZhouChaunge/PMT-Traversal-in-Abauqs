# -*- coding: utf-8 -*-
import os
import csv
from odbAccess import openOdb

models_path = os.path.join(os.getcwd(), 'abaqus', 'models')  # 拼接models的文件夹下的路径
files = os.listdir(models_path)                              # 获取‘models’文件夹下的所有文件名
odb_files = [f for f in files if f.endswith('.odb')]         # 获取'models'文件夹下所有的odb文件
current_odb_index = len(odb_files)                           # 获取'models'文件夹下odb文件的数量

odb_filename = 'job-{}.odb'.format(current_odb_index)
csv_filename = 'job-{}.csv'.format(current_odb_index)
odb_filepath =os.path.join(models_path, odb_filename)
csv_filepath =os.path.join(models_path, csv_filename)
my_odb = openOdb(odb_filepath, readOnly=False)
step = my_odb.steps['load']


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