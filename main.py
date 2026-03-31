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


# ============================================================
# Create output directories / 创建输出文件夹
# ============================================================
OUTPUTS_FOLDER = 'outputs'
RESULTS_FOLDER = 'results'
os.makedirs(OUTPUTS_FOLDER, exist_ok=True)
os.makedirs(RESULTS_FOLDER, exist_ok=True)


def read_and_modify_inp_file(m, a, i, material_line_index, inp_lines):
    """
    Overwrite constitutive parameters in the .inp file and save as a new job file.
    覆写本构参数（m, a）到 .inp 文件中，并保存为新的 Job 文件。
    """
    line_parts = inp_lines[material_line_index].split(',')
    line_parts[5] = str(round(m, 2))   # Overwrite parameter m / 覆写本构参数 m
    line_parts[6] = str(round(a, 3))   # Overwrite parameter a / 覆写本构参数 a
    inp_lines[material_line_index] = ','.join(line_parts)

    job_path = os.path.join('abaqus', 'models', f'Job-{i}.inp')
    with open(job_path, 'w') as f:
        f.write('\n'.join(inp_lines))
    print(f'--- Step1: Job-{i}.inp constitutive parameters overwritten / 本构参数覆写完成 ---\n')


def submit_abaqus_job(i):
    """
    Generate a batch file and submit the Abaqus job for computation.
    生成批处理文件并提交 Abaqus 计算任务。
    """
    models_dir = os.path.join(os.getcwd(), 'abaqus', 'models')
    inputs_path = os.path.join(os.getcwd(), 'abaqus', 'inputs.bat')

    with open(inputs_path, 'w') as input_file:
        # GPU acceleration (default). To use CPU instead, replace gpus=1 with cpus=6.
        # GPU 加速（默认）。若使用 CPU，将 gpus=1 替换为 cpus=6。
        input_file.write(f'call abaqus job=Job-{i} int user=CYCLIC.for gpus=1 ask=off')

    start_time = time.time()
    subprocess.run(inputs_path, cwd=models_dir)
    elapsed = time.time() - start_time
    print(f'--- Step2: Job-{i}.inp completed in {elapsed:.1f}s / 计算完成，用时 {elapsed:.1f}s ---\n')


def extract_results_from_odb(script_path, i):
    """
    Extract PMT curve results from the .odb file by running the extraction script.
    通过执行提取脚本，从 .odb 文件中提取旁压曲线结果。
    """
    subprocess.call(['python', script_path])
    print(f'--- Step3: Job-{i} nodal displacement data extracted / 节点位移数据提取完成 ---\n')


def process_and_visualize_data(m, a, i, output_path, real_data_path, r0, h0, V0, Pv, P0):
    """
    Compare simulation results with in-situ data, compute error metrics, and generate plots.
    将模拟结果与实测数据对比，计算误差指标，并生成可视化图表。

    Args / 参数:
        m, a: Constitutive parameters / 本构参数
        i: Job index / 任务编号
        output_path: Path to model output CSVs / 模型输出 CSV 路径
        real_data_path: Path to in-situ measurement data / 实测数据路径
        r0: Probe cavity radius (cm) / 旁压器压力腔半径
        h0: Probe cavity length (cm) / 旁压器压力腔长度
        V0: Initial volume at P0 / 初始压力 P0 对应的体积
        Pv: Applied pressure (kPa) / 施加压力
        P0: In-situ overburden stress (kPa) / 原位自重应力

    Returns / 返回:
        Tuple of (m, a, R2, Pearson correlation, deviation rate)
    """
    # Read in-situ measurement data / 读取实测数据
    real_data = pd.read_excel(real_data_path, header=None, skiprows=1)
    real_P, real_V = real_data[0].values, real_data[1].values

    # Read simulation results / 读取模拟结果
    csv_file_path = os.path.join(output_path, f'job-{i}.csv')
    csv_df = pd.read_csv(csv_file_path)

    # Compute normalized pressure and volume / 计算归一化压力和体积
    csv_df['P_norm'] = (csv_df['Time'] * Pv) / P0
    csv_df['V'] = (math.pi * ((r0 + csv_df['U1'] * 100) ** 2 - r0 ** 2) * h0) + V0
    pred_P, pred_V = csv_df['P_norm'].values, csv_df['V'].values

    # Save fitted data / 保存拟合数据
    job_output_dir = os.path.join(OUTPUTS_FOLDER, f'Job-{i}_output')
    os.makedirs(job_output_dir, exist_ok=True)

    pred_df = pd.DataFrame({'P': pred_P, 'V': pred_V})
    pred_df.to_excel(
        os.path.join(job_output_dir, f'Fitted_Data(m={round(m, 2)}&a={round(a, 3)}).xlsx'),
        index=False
    )

    # Plot P-V curves / 绘制 P-V 曲线
    plt.clf()
    plt.plot(real_P, real_V, label='Real Curve')
    plt.plot(pred_P, pred_V, label='Fitted Curve')
    plt.xlabel('P/P0')
    plt.ylabel('V')
    plt.title(f'P-V Curve (m={round(m, 2)} & a={round(a, 3)})')
    plt.legend()
    plt.savefig(os.path.join(job_output_dir, f'PV_curve(m={round(m, 2)}&a={round(a, 3)}).png'))
    plt.close()

    print(f'--- Step4: Evaluating fit quality for Job-{i} / 对 Job-{i} 拟合结果进行误差评价 ---\n')

    # Determine the start of the quasi-elastic segment / 获取类弹性段起始位置
    dV = np.diff(real_V)
    ddV = np.diff(dV)
    positive_idx = np.argwhere(ddV > 0)
    if positive_idx.size > 0:
        start_idx = int(positive_idx[0])
        print(f'Effective comparison starts at point {start_idx + 1} / 有效对比段从第 {start_idx + 1} 个点开始')
    else:
        idx = input(
            "No clear quasi-elastic segment found. "
            "Please enter the starting index based on the curve:\n"
            "该曲线没有明显的类弹性段，请根据实际图像输入起始点序号索引："
        )
        while not idx.isdigit():
            idx = input("Please enter a numeric index / 请输入一个数字索引：")
        print(f"Effective comparison starts at point {idx} / 有效对比段从第 {idx} 个点开始")
        start_idx = int(idx) - 1

    # Determine the end of the effective comparison range / 获取有效对比段终止位置
    real_endpoint = real_P[-1]
    pred_endpoint = pred_P[-1]
    if pred_endpoint > real_endpoint:
        end_idx = len(real_P) - 1
        print(f"Effective comparison ends at point {len(real_P)} / 有效对比段在第 {len(real_P)} 个点结束\n")
    else:
        closest_elements = sorted(real_P, key=lambda x: abs(x - pred_endpoint))[:2]
        closest_elements.sort()
        end_idx = np.abs(real_P - closest_elements[0]).argmin()
        print(f"Effective comparison ends at point {end_idx + 1} / 有效对比段在第 {end_idx + 1} 个点结束\n")

    # Interpolate and compute error metrics / 插值并计算误差指标
    inter_pred_f = interp1d(pred_P, pred_V, kind='linear')
    effective_real_V = real_V[start_idx:end_idx + 1]
    effective_pred_V = inter_pred_f(real_P[start_idx:end_idx + 1])

    r2 = r2_score(effective_real_V, effective_pred_V)
    corr = np.corrcoef(effective_real_V, effective_pred_V)[0, 1]
    deviation_rate = np.mean(np.abs(effective_pred_V - effective_real_V) / effective_real_V)

    print(f'R2 / 可决系数 R2: {r2:.6f}')
    print(f'Pearson Correlation / 皮尔逊相关系数: {corr:.6f}')
    print(f'Deviation Rate / 误差率: {deviation_rate * 100:.4f}%')

    return m, a, r2, corr, deviation_rate


def plot_error_surface(m_values, a_values, deviation_values):
    """
    Visualize the error surface over the m-a parameter space.
    在 m-a 参数空间上可视化误差曲面。
    """
    df = pd.DataFrame({'m': m_values, 'a': a_values, 'deviation': deviation_values})
    m = df['m'].values
    a = df['a'].values
    e = 100 * df['deviation'].values  # Convert to percentage / 转换为百分比

    m_grid, a_grid = np.mgrid[m.min():m.max():100j, a.min():a.max():100j]
    e_grid = griddata((m, a), e, (m_grid, a_grid), method='cubic')

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(m_grid, a_grid, e_grid, cmap='viridis', edgecolor='none')
    ax.set_xlabel('m')
    ax.set_ylabel('a')
    ax.set_zlabel('Deviation (%)')
    ax.set_title('3D Error Surface over Parameter Space (m, a)')
    plt.savefig(os.path.join(RESULTS_FOLDER, 'error_surface.png'), dpi=300)
    plt.show()


def main():
    """
    Main function: traverse parameter space, run Abaqus simulations, and collect results.
    主函数：遍历参数空间，运行 Abaqus 模拟，收集结果。
    """
    # Read the base .inp model and locate the material parameter line
    # 读取基础 .inp 模型，定位本构参数所在行
    with open(os.path.join('abaqus', 'basic model.inp'), 'r') as f:
        inp_lines = f.read().splitlines()

    material_line_index = None
    for idx, line in enumerate(inp_lines):
        if '*User Material, constants=9' in line:
            material_line_index = idx + 1
            break

    if material_line_index is None:
        raise ValueError("Cannot find '*User Material, constants=9' in inp file / 无法在 inp 文件中找到材料参数行")

    # Define parameter traversal range / 定义参数遍历范围
    error_rate = 0.2                     # Allowable deviation rate / 容许误差率
    m_pred, a_pred = 7.5, 0.45           # Baseline values / 预测基准值
    m_min = m_pred * (1 - error_rate)
    m_max = m_pred * (1 + error_rate)
    m_step = 0.2
    a_min = a_pred * (1 - error_rate)
    a_max = a_pred * (1 + error_rate)
    a_step = 0.01

    m_array = np.arange(m_min, m_max + m_step, m_step)
    a_array = np.arange(a_min, a_max + a_step, a_step)

    # PMT test parameters / 旁压试验参数
    output_path = os.path.join(os.getcwd(), 'abaqus', 'models')
    real_data_path = os.path.join(os.getcwd(), 'test data', 'insitu_data.xlsx')
    Height = 15     # Test depth (m) / 试验深度
    r0 = 4.5        # Probe cavity radius (cm) / 压力腔半径
    h0 = 21          # Probe cavity length (cm) / 压力腔长度
    V0 = 85          # Initial volume at P0 / P0 对应初始体积
    Pv = 209         # Applied pressure (kPa) / 施加压力
    P0 = 18 * 2 + (Height - 2) * 20  # Overburden stress (kPa) / 原位自重应力

    # Result collectors / 结果收集
    m_values, a_values, r2_values, pearson_values, deviation_values = [], [], [], [], []

    # Main traversal loop / 主遍历循环
    extract_script = os.path.join(os.getcwd(), 'abaqus', 'run_extract.py')
    for i, (m, a) in enumerate(product(m_array, a_array), start=1):
        read_and_modify_inp_file(m, a, i, material_line_index, inp_lines.copy())
        submit_abaqus_job(i)
        extract_results_from_odb(extract_script, i)
        m_val, a_val, r2, pearson, deviation = process_and_visualize_data(
            m, a, i, output_path, real_data_path, r0, h0, V0, Pv, P0
        )
        m_values.append(m_val)
        a_values.append(a_val)
        r2_values.append(r2)
        pearson_values.append(pearson)
        deviation_values.append(deviation)

    # Export summary statistics / 导出汇总统计结果
    result_df = pd.DataFrame({
        'm': m_values,
        'a': a_values,
        'R2': r2_values,
        'Pearson': pearson_values,
        'Deviation': deviation_values
    })
    result_df.to_excel(os.path.join(RESULTS_FOLDER, 'Output_Statistical.xlsx'), sheet_name='Sheet1', index=False)
    plot_error_surface(m_values, a_values, deviation_values)


if __name__ == "__main__":
    main()
