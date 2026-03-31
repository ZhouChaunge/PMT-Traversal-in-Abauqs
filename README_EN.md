<div align="center">

# 🧪 PMT-Traversal-in-Abaqus

**Pressuremeter Test Parameter Traversal Framework Based on Abaqus UMAT Subroutine**

[![Python 3.7](https://img.shields.io/badge/Python-3.7-3776AB?logo=python&logoColor=white)](https://www.python.org/downloads/release/python-370/)
[![Abaqus 2021](https://img.shields.io/badge/Abaqus-2021-FF6600?logo=dassaultsystemes&logoColor=white)](https://www.3ds.com/products-services/simulia/products/abaqus/)
[![Fortran](https://img.shields.io/badge/Fortran-UMAT-734F96?logo=fortran&logoColor=white)](#)
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](LICENSE)

**English** **|** [**中文**](README.md)

---

*Automated constitutive parameter traversal · GPU-accelerated FEA · Error evaluation & visualization*

</div>

---

## 📋 Table of Contents

- [Introduction](#-introduction)
- [Workflow](#-workflow)
- [Environment Setup](#️-environment-setup)
- [Project Structure](#-project-structure)
- [Quick Start](#-quick-start)
- [GPU vs CPU](#-gpu-vs-cpu)
- [References](#-references)

---

## 🔍 Introduction

This project uses Python to perform **Pressuremeter Test (PMT, axisymmetric model) traversal calculations** in Abaqus with a UMAT subroutine, automating the entire pipeline from parameter overwriting to result visualization.

### Key Features

| Step | Feature | Description |
|:---:|:---|:---|
| ① | **Parameter Overwriting** | Automatically modify constitutive parameters (m, a) in `.inp` files |
| ② | **Auto Submission** | Submit Abaqus jobs via `.bat` batch files with GPU acceleration support |
| ③ | **Result Extraction** | Automatically extract nodal displacement data from `.odb` files |
| ④ | **Post-Processing** | Error evaluation (R², Pearson, deviation rate) + P-V curve visualization |

---

## 🔄 Workflow

<div align="center">
<img src="https://raw.githubusercontent.com/ZhouChaunge/PMT-Traversal-in-Abauqs/main/images/Pseudo-code.png" alt="Program Flowchart" width="620"/>
</div>

---

## ⚙️ Environment Setup

### 1. Python Environment

```bash
# Create a Python 3.7 virtual environment
conda create -n pyabaqus37 python=3.7

# Activate the environment
conda activate pyabaqus37

# Install dependencies
pip install -r requirements.txt
```

<details>
<summary>📺 <b>Click to view installation demo GIF</b></summary>
<br/>
<div align="center">
<img src="https://raw.githubusercontent.com/ZhouChaunge/PMT-Traversal-in-Abauqs/main/images/py_env_install.gif" alt="Python Environment Installation Demo" width="700"/>
</div>
</details>

### 2. Abaqus & Associated Software

| Software | Recommended Version | Purpose |
|:---|:---|:---|
| **Abaqus** | 2021 | Finite Element Analysis |
| **Visual Studio** | 2019 | C/C++ Build Environment |
| **Intel Parallel Studio XE** | 2020 | Fortran Compiler |
| **NVIDIA CUDA + cuDNN** | Match your GPU | GPU Acceleration |

> ⚠️ The Fortran compiler must be correctly linked to use the UMAT subroutine. See [References](#-references) for installation guides.

---

## 📁 Project Structure

```
PMT-Traversal-in-Abaqus/
├── 📄 main.py                  Main entry point
├── 📄 requirements.txt         Python dependencies
├── 📄 README.md                Chinese documentation
├── 📄 README_EN.md             English documentation (this file)
│
├── 📂 abaqus/
│   ├── basic model.inp         Abaqus base model (template for param overwriting)
│   ├── inputs.bat              Batch submission script
│   ├── extract.py              ODB result extraction (Abaqus Python script)
│   ├── run_extract.py          Runner for extract.py
│   └── 📂 models/
│       └── CYCLIC.for          UMAT subroutine (Shanghai Model, Fortran)
│
├── 📂 test data/
│   └── insitu_data.xlsx        In-situ PMT measurement data
│
├── 📂 outputs/                 Computation results & fitted curves per parameter set
└── 📂 results/                 Summary statistics + error surface plot
```

---

## 🚀 Quick Start

> ⚠️ **Path Requirement**: Abaqus does not support non-ASCII characters in file paths. Ensure the project is in an **English-only path**.

**1.** Place in-situ PMT data in `test data/insitu_data.xlsx`

**2.** Configure test parameters in `main.py`:

```python
Height = 15     # Test depth (m)
r0 = 4.5        # Probe cavity radius (cm)
h0 = 21         # Probe cavity length (cm)
V0 = 85         # Initial volume at P0
Pv = 209        # Applied pressure (kPa)

# Parameter traversal range
m_pred, a_pred = 7.5, 0.45    # Baseline values
error_rate = 0.2               # Allowable deviation rate ±20%
```

**3.** Run:

```bash
conda activate pyabaqus37
python main.py
```

**4.** Check results:

| Output Directory | Contents |
|:---|:---|
| `outputs/Job-{i}_output/` | Fitted data (.xlsx) and P-V curve plot (.png) per parameter set |
| `results/Output_Statistical.xlsx` | R², Pearson, and deviation rate summary for all combinations |
| `results/error_surface.png` | 3D error surface over the m-a parameter space |

---

## 🎮 GPU vs CPU

GPU acceleration is enabled by default. Switch in `main.py`:

```python
# ✅ GPU acceleration (default)
input_file.write(f'call abaqus job=Job-{i} int user=CYCLIC.for gpus=1 ask=off')

# 📌 CPU multi-core parallel (alternative)
input_file.write(f'call abaqus job=Job-{i} int user=CYCLIC.for cpus=6 ask=off')
```

> 💡 **Performance Reference** (laptop benchmark)
>
> | Hardware | Configuration |
> |:---|:---|
> | CPU | i7-12700H (14 cores, 20 threads) |
> | GPU | RTX 3060 Laptop (6GB VRAM) |
>
> GPU computation saves approximately **40%** time compared to CPU multi-core parallel processing.

---

## 📚 References

| Topic | Link |
|:---|:---|
| Abaqus + VS + Fortran Setup | [Bilibili Tutorial](https://www.bilibili.com/video/BV1Mj411r7WQ) |
| Installation Troubleshooting | [Bilibili Notes](https://www.bilibili.com/video/BV1Sa411g744) |
| CUDA / cuDNN Setup | [Bilibili CUDA Guide](https://www.bilibili.com/video/BV1nL4y1b7oT) |
| Abaqus GPU Acceleration | [Bilibili GPU Guide](https://www.bilibili.com/video/BV1vT4y1z74H) |

---

<div align="center">

**If this project is helpful, please give it a ⭐ Star!**

Made with ❤️ for Geotechnical Engineering

</div>
