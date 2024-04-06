# -*- coding: utf-8 -*-
import os
import time
import math
import numpy as np
import pandas as pd
import subprocess
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from scipy.interpolate import griddata
from sklearn.metrics import r2_score
from itertools import product


# TODO 分别创建outputs和results文件夹
outputs_folder = 'outputs'
results_folder = 'results'
if not os.path.exists(outputs_folder):
    os.makedirs(outputs_folder)
if not os.path.exists(results_folder):
    os.makedirs(results_folder)


# TODO 覆写本构参数，生成inp文件
def read_and_modify_inp_file(m, a, i, material_line_index, inp_lines):
    line_parts = inp_lines[material_line_index].split(',')
    line_parts[5] = str(round(m, 2))  # 覆写本构参数m
    line_parts[6] = str(round(a, 3))  # 覆写本构参数a
    modified_line = ','.join(line_parts)
    inp_lines[material_line_index] = modified_line

    with open(f'abaqus\\models\\Job-{i}.inp', 'w') as f:
        f.write('\n'.join(inp_lines))
    print(f'------------------------------Step1: Job-{i}.inp文件本构参数覆写完成!------------------------------\n')


# TODO 提交inp文件,执行abaqus计算
def submit_abaqus_job(i):
    abaqus_filename, models_filename, inputs_filename = 'abaqus', 'models', 'inputs.bat'
    working_filepath = os.path.join(os.getcwd(), abaqus_filename, models_filename)
    inputs_filepath = os.path.join(os.getcwd(), abaqus_filename, inputs_filename)
    with open(inputs_filepath, 'w') as input_file:
        '''
        这里使用了gpu进行了加速计算，abaqus2019版本之后支持gpu加速，在本算例当中，gpu的计算速度相较于cpu可以提升约40%以上
        若本地未部署gpu加速相关设置，可仍采用cpu进行计算，将input_file.write的代码改为如下所示即可：
        input_file.write(f'call abaqus job=Job-{i} int user=CYCLIC.for cpus=6 ask=off')
        '''
        input_file.write(f'call abaqus job=Job-{i} int user=CYCLIC.for gpus=1 ask=off')

    run_start_time = time.time()
    subprocess.run(inputs_filepath, cwd=working_filepath)
    run_end_time = time.time()
    time_used = run_end_time - run_start_time
    print(f'------------------------------Step2: Job-{i}.inp计算完成！用时{time_used}------------------------------\n')


# TODO 提取odb文件中的旁压曲线结果，生成csv文件
def extract_results_from_odb(script_path, i):
    subprocess.call(['python', script_path])
    print(f'------------------------------Step3: job-{i}.obd节点位移数据提取完成------------------------------\n')


# TODO 计算曲线与实测数据之间的误差统计
def process_and_visualize_data(m, a, i, output_path, real_data_path, r0, h0, V0, Pv, P0):
    real_data = pd.read_excel(real_data_path, header=None, skiprows=1)
    real_P, real_V = real_data[0].values, real_data[1].values
    csv_file_path = os.path.join(output_path, f'job-{i}.csv')
    csv_file_df = pd.read_csv(csv_file_path)
    csv_file_df['P_norm'] = (csv_file_df['Time'] * Pv) / P0
    csv_file_df['V'] = (math.pi * ((r0 + ((csv_file_df['U1']) * 100)) ** 2 - r0 ** 2) * h0) + V0
    pred_P, pred_V = csv_file_df['P_norm'].values, csv_file_df['V'].values

    if not os.path.exists(f'outputs/Job-{i}_output'):
        os.mkdir(f'outputs/Job-{i}_output')

    pred_df = pd.DataFrame({'P': pred_P, 'V': pred_V})
    pred_df.to_excel(f'outputs/Job-{i}_output/Fitted_Data(m={round(m, 2)}&a={round(a, 3)}).xlsx', index=False)

    plt.clf()  # 清除当前图形窗口
    plt.plot(real_P, real_V, label='Real Curve')
    plt.plot(pred_P, pred_V, label='Fitted Curve')
    plt.xlabel('P/P0')
    plt.ylabel('V')
    plt.title(f'P-V Curve(m={round(m, 2)}&a={round(a, 3)})')
    plt.savefig(f'outputs/Job-{i}_output/PV_curve(m={round(m, 2)}&a={round(a, 3)}).png')
    plt.legend()
    plt.close()

    print(f'------------------------------Step4: 对job-{i}的拟合结果进行误差评价------------------------------\n')
    # 获取类弹性段开始位置
    dV = np.diff(real_V)
    ddV = np.diff(dV)
    positive_idx = np.argwhere(ddV > 0)
    if positive_idx.size > 0:
        start_idx = int(positive_idx[0])
        print(f'曲线的有效对比段在第{start_idx + 1}个点开始')
    else:
        idx = input("该曲线没有明显的类弹性段，请根据实际图像输入类弹性段开始点的序号索引：")
        while not idx.isdigit():
            idx = input("请输入一个数字索引：")
        print(f"曲线的有效对比段在第{idx}个点开始")
        start_idx = int(idx) - 1

    # 获取曲线的最后有效位置
    real_endpoint = real_P[-1]
    pred_endpoint = pred_P[-1]
    if pred_endpoint > real_endpoint:
        end_idx = len(real_P) - 1
        print(f"曲线的有效对比段在第{len(real_P)}个点结束\n")
    else:
        closest_elements = sorted(real_P, key=lambda x: abs(x - pred_endpoint))[:2]
        closest_elements.sort()
        end_idx = np.abs(real_P - closest_elements[0]).argmin()
        print(f"曲线的有效对比段在第{end_idx + 1}个点结束\n")

    inter_pred_f = interp1d(pred_P, pred_V, kind='linear')
    RealV_in_Fitted = inter_pred_f(real_P[start_idx:end_idx + 1])
    effective_real_V = real_V[start_idx:end_idx + 1]
    effective_pred_V = RealV_in_Fitted

    r2 = r2_score(effective_real_V, effective_pred_V)
    corr = np.corrcoef(effective_real_V, effective_pred_V)
    real_data, pred_data = np.array(effective_real_V), np.array(effective_pred_V)
    deviation_rates = abs(pred_data - real_data) / real_data
    deviation_rate = np.mean(deviation_rates)
    print(f'可决系数r2为：{r2}')
    print(f'皮尔逊相关系数Cor为：{corr[0, 1]}')
    print(f'误差率为：{deviation_rate * 100}%')

    return m, a, r2, corr[0, 1], deviation_rate


# TODO 可视化处理，绘制 m-a 参数空间上的误差曲面
def plot_error_surface(m_values, a_values, deviation_values):
    df = pd.DataFrame({'m': m_values, 'a': a_values, 'deviation': deviation_values})
    m = df['m'].values
    a = df['a'].values
    e = 100 * df['deviation'].values  # 转换为百分数

    m_min, m_max = np.min(m), np.max(m)
    a_min, a_max = np.min(a), np.max(a)
    m_grid, a_grid = np.mgrid[m_min:m_max:100j, a_min:a_max:100j]

    e_grid = griddata((m, a), e, (m_grid, a_grid), method='cubic')

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(m_grid, a_grid, e_grid, cmap='viridis', edgecolor='none')
    ax.set_xlabel('m')
    ax.set_ylabel('a')
    ax.set_zlabel('Deviation (e)')
    ax.set_title('3D Interpolated Surface Plot of Deviation (e) over m and a')
    plt.savefig(f'results/error_surface.png', dpi=1200)
    plt.show()


# TODO 主函数部分
def main():
    with open('abaqus\\basic model.inp', 'r') as f:     # 打开基础inp模型，获取本构参数的行编号
        inp_contents = f.read()

    inp_lines = inp_contents.splitlines()
    material_line_index = None
    for line_index, line in enumerate(inp_lines):
        if '*User Material, constants=9' in line:
            material_line_index = line_index + 1

    # 定义输出所需的变量
    m_values = []
    a_values = []
    r2_values = []
    Pearson_values = []
    deviation_values = []

    error_rate = 0.2  # 定义参数的容许误差率
    m_pred, a_pred = 7.5, 0.45  # 设置参数m和a的预测基准值
    m_min, m_max, m_step = m_pred*(1-error_rate), m_pred*(1+error_rate), 0.2
    a_min, a_max, a_step = a_pred*(1-error_rate), a_pred*(1+error_rate), 0.01

    m_array = np.arange(m_min, m_max + m_step, m_step)
    a_array = np.arange(a_min, a_max + a_step, a_step)

    # 设置旁压试验的参数
    output_path = os.path.join(os.getcwd(), 'abaqus', 'models')
    real_data_path = os.path.join(os.getcwd(), 'test data', 'insitu_data.xlsx')
    Height = 15  # 设置旁压试验的深度（单位m）
    r0 = 4.5     # 设置旁压器压力腔的半径（单位cm）
    h0 = 21      # 设置旁压器压力腔的长度（单位cm）
    V0 = 85      # 在旁压曲线中(P0,V0)为第一特征点，V0是初始压力P0对应的体积
    Pv = 209     # 设置旁压试验加载的压力大小（单位kPa)
    P0 = 18 * 2 + (Height - 2) * 20  # 计算旁压试验深度的原位自重应力

    # 执行主循环
    for i, (m, a) in enumerate(product(m_array, a_array), start=1):
        read_and_modify_inp_file(m, a, i, material_line_index, inp_lines.copy())
        submit_abaqus_job(i)
        extract_results_from_odb(os.path.join(os.getcwd(), "abaqus", "run_extract.py"), i)
        m, a, r2, Pearson, deviation = process_and_visualize_data(
            m, a, i, output_path, real_data_path, r0, h0, V0, Pv, P0
        )
        m_values.append(m)
        a_values.append(a)
        r2_values.append(r2)
        Pearson_values.append(Pearson)
        deviation_values.append(deviation)

    # 输出统计结果到Excel
    result_df = pd.DataFrame({
        'm': m_values,
        'a': a_values,
        'r2': r2_values,
        'Pearson': Pearson_values,
        'deviation': deviation_values
    })
    if not os.path.exists('results'):
        os.mkdir('results')
    result_df.to_excel('results/Output_Statistical.xlsx', sheet_name='Sheet1', index=False)
    plot_error_surface(m_values, a_values, deviation_values)


if __name__ == "__main__":
    main()