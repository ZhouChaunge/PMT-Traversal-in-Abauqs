<div align="center">
<h1 align="center">PMT-Traversal-in-Abauqs</h1>
</div>

<br>
<b>项目简介</b>：利用python进行在Abaqus中基于Umat子程序进行旁压试验(轴对称模型)的遍历计算，主要流程包括<b>覆写本构参数</b>、<b>自动提交运算</b>，<b>结果导出</b>和<b>后处理</b>等，在执行求解过程中调用了gpu加速有限元计算。

[images/Pseudo-code.png](https://raw.githubusercontent.com/ZhouChaunge/PMT-Traversal-in-Abauqs/main/images/Pseudo-code.png)
<br>


# PMT-Traversal-in-Abauqs

1. 配置运行环境。在当前文件夹目录下，进入cmd控制台，然后执行下面的命令
   
    conda crate -n pyabaqus37 python=3.7
   
    conda activate pyabaqus37
   
    pip install -r requirements.txt
   

3. 运行程序 main.py


注意：

(1).需要在全英文路径下运行;

(2).默认abaqus是使用gpu进行计算的，若尚未部署gpu运算环境，需要在main.py中进行修改

(3).每组参数的运算结果在文件夹“outputs”中；所有结果的统计数据在文件夹"results"中
