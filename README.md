<div align="center">
<h1 align="center">PMT-Traversal-in-Abauqs</h1>
</div>
<br>

## 0.项目简介

> 利用python进行在Abaqus中基于Umat子程序进行旁压试验(轴对称模型)的遍历计算，主要流程包括<b>覆写本构参数</b>、<b>自动提交运算</b>，<b>结果导出</b>和<b>后处理</b>等，在执行求解过程中调用了gpu加速有限元计算。程序伪代码如下：

![images/Pseudo-code.png](https://raw.githubusercontent.com/ZhouChaunge/PMT-Traversal-in-Abauqs/main/images/Pseudo-code.png)
<br>

## 1.环境配置

环境配置分为两个部分，第一个是python环境配置，第二个是abaqus环境配置。

### 1.1 python环境配置

下载项目后，在项目当前文件夹目录(...\PMT-Traversal-in-Abauqs-main)下，进入cmd控制台，然后执行下面的命令：

    conda crate -n pyabaqus37 python=3.7 #创建名为pyabaqus37的python环境

    conda activate pyabaqus37            # 在当前目录下激活pyabaqus37环境
    
    pip install -r requirements.txt      # 安装依赖包

上述过程如下所示：
![images/py_env_install.gif](https://raw.githubusercontent.com/ZhouChaunge/PMT-Traversal-in-Abauqs/main/images/py_env_install.gif)

完成安装后，python环境配置完成。

### 1.2 abaqus环境配置

在本项目中，使用了gpu加速有限元计算，因此需要安装支持gpu运算的abaqus版本。由于在执行有限元计算时需要调用Umat子程序，所以也需要安装fortran编译的关联软件，可以考虑安装如下版本：

    Abaqus 2021
    Visual Stuido 2019
    Intel Parallel Studio XE 2020
关于abaqus软件安装和关联子程序方面的内容，可以参考以下资料
- [Abaqus2021+vs2019+fortran2020子程序全过程安装关联视频教程](https://www.bilibili.com/video/BV1Mj411r7WQ)
- [abaqus2021+vs2019+fortran2020安装及关联时，一些问题的解决记录](https://www.bilibili.com/video/BV1Sa411g744)

在完成abaqus及关联程序的安装后，需要配置gpu加速的相关设置，与深度学习配置gpu加速相同，需要安装nvidia的cuda和cudnn，关于nvidia cuda和cudnn的安装，可以参考以下资料：
- [【CUDA安装/多CUDA兼容】Windows深度学习环境配置](https://www.bilibili.com/video/BV1nL4y1b7oT)
- [ABAQUS调用NVIDIA显卡CUDA加速](https://www.bilibili.com/video/BV1vT4y1z74H)


### 1.3 项目文件解释

#### 1.3.1 主目录

1. 因为abaqus执行运算时不支持中文路径，需要在全英文路径下运行
   
2. 默认abaqus是使用gpu进行计算的，若尚未部署gpu运算环境，需要在main.py中进行修改，下面分别给出gpu和cpu运算的代码，请您根据自己的实际需求进行修改
   
        input_file.write(f'call abaqus job=Job-{i} int user=CYCLIC.for gpus=1 ask=off')  # 调用gpu执行计算
        input_file.write(f'call abaqus job=Job-{i} int user=CYCLIC.for cpus=6 ask=off')  # 调用cpu执行计算
    这里作者给出一点参考，以本人笔记本的基本配置，在执行gpu运算时，单个模型的运算时间相较于cpu多核并行计算可以节省约40%左右。

        CPU型号：12th Gen Intel(R)Core(TM)i7-12700H  14核20线程
        GPU型号：Nvidia GeForce RTX 3060 Laptop      6GB显存
   
3. 文件夹“test data”用于存放现场实测的数据，每次运行后依据insitu_data.xlsx中的数据进行误差计算，并将结果输出到outputs文件夹中

4. 每组参数的运算结果在文件夹“outputs”中；所有结果的统计数据在文件夹"results"中


#### 1.3.2 abaqus文件夹
1. basic model.inp是一个abaqus基础模型，每次运算均是在基础模型上修改本构参数，进而利用bat批处理文件提交计算，然后执行run_extract.py进行后处理。

2. 在子文件夹models中，包含了一个fortran编写的Umat子程序（上海模型）

3. abaqus文件夹中涉及到的文件解释如下：

    |文件名称        | 文件类型      |解释|
    |---------------|---------------|--------------------------|
    |basic model.inp|基本模型       |参数的改写均是基于该文件进行|
    |models         |子文件夹       |用于存放计算后的文件       |
    |input.bat      |批处理程序     |用于提交inp文件执行计算    |
    |extract.py     |结果提取文件   |用于从odb文件中提取结果    |
    |run_extract.py |结果提取文件   |用于执行extract.py文件     |

