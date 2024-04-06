# PMT-Traversal-in-Abauqs

1. 配置运行环境。在当前文件夹目录下，执行
   
    conda crate -n pyabaqus37 python=3.7
   
    conda activate pyabaqus37
   
    pip install -r requirements.txt
   

3. 运行程序 main.py


注意：

(1).需要在全英文路径下运行;

(2).默认abaqus是使用gpu进行计算的，若尚未部署gpu运算环境，需要在main.py中进行修改

(3).每组参数的运算结果在文件夹“outputs”中；所有结果的统计数据在文件夹"results"中
