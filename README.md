<div align="center">

# 🧪 PMT-Traversal-in-Abaqus

**基于 Abaqus UMAT 子程序的旁压试验参数遍历计算框架**

[![Python 3.7](https://img.shields.io/badge/Python-3.7-3776AB?logo=python&logoColor=white)](https://www.python.org/downloads/release/python-370/)
[![Abaqus 2021](https://img.shields.io/badge/Abaqus-2021-FF6600?logo=dassaultsystemes&logoColor=white)](https://www.3ds.com/products-services/simulia/products/abaqus/)
[![Fortran](https://img.shields.io/badge/Fortran-UMAT-734F96?logo=fortran&logoColor=white)](#)
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](LICENSE)

[**English**](README_EN.md) **|** **中文**

---

*自动化遍历本构参数 · GPU 加速有限元计算 · 误差评价与可视化*

</div>

---

## 📋 目录

- [项目简介](#-项目简介)
- [程序流程](#-程序流程)
- [环境配置](#️-环境配置)
- [项目结构](#-项目结构)
- [快速开始](#-快速开始)
- [GPU vs CPU](#-gpu-vs-cpu)
- [参考资料](#-参考资料)

---

## 🔍 项目简介

利用 Python 在 Abaqus 中基于 UMAT 子程序进行 **旁压试验（PMT，轴对称模型）的遍历计算**，自动完成从参数覆写到结果可视化的完整流程。

### 核心功能

| 步骤 | 功能 | 说明 |
|:---:|:---|:---|
| ① | **参数覆写** | 自动修改 `.inp` 文件中的本构参数 (m, a) |
| ② | **自动提交** | 通过 `.bat` 批处理提交 Abaqus 运算，支持 GPU 加速 |
| ③ | **结果提取** | 从 `.odb` 文件中自动提取节点位移数据 |
| ④ | **后处理** | 误差评价（R², Pearson, 偏差率）+ P-V 曲线可视化 |

---

## 🔄 程序流程

<div align="center">
<img src="https://raw.githubusercontent.com/ZhouChaunge/PMT-Traversal-in-Abauqs/main/images/Pseudo-code.png" alt="程序流程图" width="620"/>
</div>

---

## ⚙️ 环境配置

### 1. Python 环境

```bash
# 创建 Python 3.7 虚拟环境
conda create -n pyabaqus37 python=3.7

# 激活环境
conda activate pyabaqus37

# 安装依赖包
pip install -r requirements.txt
```

<details>
<summary>📺 <b>点击查看安装演示 GIF</b></summary>
<br/>
<div align="center">
<img src="https://raw.githubusercontent.com/ZhouChaunge/PMT-Traversal-in-Abauqs/main/images/py_env_install.gif" alt="Python环境安装演示" width="700"/>
</div>
</details>

### 2. Abaqus 及关联软件

| 软件 | 推荐版本 | 用途 |
|:---|:---|:---|
| **Abaqus** | 2021 | 有限元分析主程序 |
| **Visual Studio** | 2019 | C/C++ 编译环境 |
| **Intel Parallel Studio XE** | 2020 | Fortran 编译器 |
| **NVIDIA CUDA + cuDNN** | 与显卡匹配 | GPU 加速计算 |

> ⚠️ 必须正确关联 Fortran 编译器才能使用 UMAT 子程序。关于安装与关联，请参见 [参考资料](#-参考资料)。

---

## 📁 项目结构

```
PMT-Traversal-in-Abaqus/
├── 📄 main.py                  主程序入口
├── 📄 requirements.txt         Python 依赖清单
├── 📄 README.md                中文文档（本文件）
├── 📄 README_EN.md             English Documentation
│
├── 📂 abaqus/
│   ├── basic model.inp         Abaqus 基础模型（参数覆写的模板）
│   ├── inputs.bat              批处理提交脚本
│   ├── extract.py              ODB 结果提取（Abaqus Python 脚本）
│   ├── run_extract.py          extract.py 调用器
│   └── 📂 models/
│       └── CYCLIC.for          UMAT 子程序（上海模型，Fortran）
│
├── 📂 test data/
│   └── insitu_data.xlsx        现场实测旁压数据
│
├── 📂 outputs/                 各组参数的计算结果与拟合曲线
└── 📂 results/                 汇总统计表 + 误差曲面图
```

---

## 🚀 快速开始

> ⚠️ **路径要求**：Abaqus 不支持中文路径，请确保项目位于**全英文路径**下。

**1.** 将实测旁压数据放入 `test data/insitu_data.xlsx`

**2.** 在 `main.py` 中配置试验参数：

```python
Height = 15     # 试验深度 (m)
r0 = 4.5        # 压力腔半径 (cm)
h0 = 21         # 压力腔长度 (cm)
V0 = 85         # P0 对应初始体积
Pv = 209        # 施加压力 (kPa)

# 参数遍历范围
m_pred, a_pred = 7.5, 0.45    # 基准值
error_rate = 0.2               # 容许偏差率 ±20%
```

**3.** 运行：

```bash
conda activate pyabaqus37
python main.py
```

**4.** 查看结果：

| 输出目录 | 内容 |
|:---|:---|
| `outputs/Job-{i}_output/` | 每组参数的拟合数据 (.xlsx) 和 P-V 曲线图 (.png) |
| `results/Output_Statistical.xlsx` | 所有参数组合的 R², Pearson, 偏差率汇总 |
| `results/error_surface.png` | m-a 参数空间的三维误差曲面 |

---

## 🎮 GPU vs CPU

默认使用 GPU 加速，在 `main.py` 中可切换：

```python
# ✅ GPU 加速（默认）
input_file.write(f'call abaqus job=Job-{i} int user=CYCLIC.for gpus=1 ask=off')

# 📌 CPU 多核并行（备选）
input_file.write(f'call abaqus job=Job-{i} int user=CYCLIC.for cpus=6 ask=off')
```

> 💡 **性能参考**（笔记本实测）
>
> | 硬件 | 配置 |
> |:---|:---|
> | CPU | i7-12700H (14核 20线程) |
> | GPU | RTX 3060 Laptop (6GB) |
>
> GPU 运算相较 CPU 多核并行可节省约 **40%** 计算时间。

---

## 📚 参考资料

| 主题 | 链接 |
|:---|:---|
| Abaqus + VS + Fortran 安装关联 | [Bilibili 视频教程](https://www.bilibili.com/video/BV1Mj411r7WQ) |
| 安装关联问题解决 | [Bilibili 问题记录](https://www.bilibili.com/video/BV1Sa411g744) |
| CUDA / cuDNN 环境配置 | [Bilibili CUDA 教程](https://www.bilibili.com/video/BV1nL4y1b7oT) |
| Abaqus GPU 加速设置 | [Bilibili GPU 教程](https://www.bilibili.com/video/BV1vT4y1z74H) |

---

<div align="center">

**如果这个项目对您有帮助，请点个 ⭐ Star 支持一下！**

Made with ❤️ for Geotechnical Engineering

</div>
