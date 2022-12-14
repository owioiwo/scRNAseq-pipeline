## A scRNA-seq analysis pipeline script with interactive mode of R
#这是一个主要基于seurat软件包的单细胞RNA测序数据分析管道的交互式R语言脚本
#主要功能：
1、读入基因表达矩阵等测序数据；
2、检测异常基因名，并能更正部分异常基因名；
3、初始化数据，提供几个去除异常细胞的功能；
4、常规数据分析管道，数据归一化、降维、聚类、数据整合、差异基因分析、细胞注释(手动或自动)
5、Rstudio界面实时展示分析所得数据和图片，文字提示操作，交互性能好，个性化定制，0编程基础人员可轻松完成整个分析流程

主要功能及使用说明如下：
a：文件说明：
1、rawdata：文件夹，放置原始基因表达矩阵数据
2、rawdata/0_inputdata_list.txt：文件，用于存放需要处理的数据的文件名，脚本将从这里读取待处理的原始数据
3、results：文件夹，存放每个步骤的输出结果，包括处理后需要整合的数据
4、results/4_integrate_list.txt：文件，用于存放需要整合的数据的文件名，脚本将从这里读取待整合的数据
5、results/5_gene_cell_analysis.txt：文件，用于存放需要做下游分析的数据的文件名，脚本将从这里读取待分析的数据
6、main.R：主脚本文件，在R平台按要求逐行运行该脚本即可，无需操作其它脚本，需按脚本输出提示和输出图片进行交互操作
7、scripts：文件夹，内含所有主脚本所需代码，勿改动！
8、README.txt：文件，即本说明文字

b：main.R主脚本详述：
1、step A之前为设置R工作环境、加载软件包和功能包

2、step A：读取数据，分析从基因表达矩阵开始，支持表达矩阵txt文件和seurat对象的rds文件，需将数据文件放入rawdata文件夹中，并将带后缀的数据文件名存放进0_inputdata_list.txt文件中，以便脚本读取数据，使用者还可以自行在R console载入数据，赋值给main.R中该行的inputdata变量即可

3、step B：分三行脚本进行单个data的处理：
stepB第一行：检测基因名是否有异常符号、是否有科学计数法格式、是否简化成日期格式，能修正日期格式基因名，输出文件名为1output_check_gene_names.rds；
stepB第二行：输入文件为上一行的输出文件，输出文件名为2output_seurat_obj.rds，做数据初始化质控，包括：
B0、要求输入project name，作为step B最后输出文件的名称，且每个细胞会带上preject name的标签，方便整合后的识别
B1、去除少于在1%细胞中表达的基因+去除表达少于200个基因的细胞，创建seurat对象，计算mt%（线粒体基因表达百分比）
B2、使用者可在此设置两个过滤阈值、低于阈值的细胞将被删除（ncount代表细胞的基因表达总数；nfeature代表细胞的基因表达总个数，默认设置均为200）
B3、对细胞的ncount、nfeature和mt%三个值进行质控，循环计算去掉中位数三倍标准偏差外的异常数据对应的细胞
B4、对细胞的ncount、nfeature表达做线性拟合，循环计算去除拟合残差三倍标准差外的异常数据对应的细胞
stepB第三行：深度质控，输入文件为上一行输出文件，输出文件名为project_name.rds，提供三个主要功能：
B5、提供seurat包含的两种数据归一化方式：LogNormalize和SCTransform，可由使用者自行选择
B6、可由使用者设置分辨率，反复查找不同数量的cluster，进行深度质控，由使用者通过输出的图片和信息确认要去除的cluster（建议高分辨率下进行2~3轮）
B7、由使用者选择合适的分辨率，按亚群进行质控，去除亚群细胞坐标中位数三倍标准差外的异常值对应的细胞（建议进行2~3轮）
注意事项：质控轮数应按实际需要进行，B6和B7每当有删除细胞操作时都需要再次做聚类分析和find cluster操作（即脚本内选项3操作，选项3操作完才算完成1轮质控），B6和B7为两个互不影响相互独立的功能（可按实际需求或单独使用、或都使用、或都不使用）

“step A” + “step B”结合可处理多个测序数据，为数据整合做准备，无需整合的数据可直接跨过“step optional”进行“step C”

4、step optional：可选项，即数据整合，对不需要整合的数据可略过此步骤，输入数据文件均需放入results文件夹内，且对应数据带后缀的文件名需存入4_integrate_list.txt文件中（每行放一个文件名），整合功能包含了seurat提供的四种锚定策略：两种归一化策略（Log、SCT）*两种不同算法（rpca、cca）（使用者可自由选取），使用者可根据最终结果重复调整整合强度（调整k.anchor参数），输出文件名为使用者填写的project name.rds，保存在results文件夹下

5、step C：下游数据分析，输入数据文件需放入results文件夹内，且对应数据带后缀的文件名需存入5_gene_cell_analysis.txt文件中（每行放一个文件名），含以下三个主要功能：
C1、由使用者重复设置分辨率找出合适的亚群后进行差异基因分析（找markers），自动保存markers文件
C2、提供seurat所含的多种基因表达图的绘制功能（由使用者选择要画的基因和图类），需手动保存图片
C3、细胞注释功能，两种方式可选：由使用者根据基因表达情况手动注释，或者由SingleR自动注释，注释后的数据文件覆盖保存为读取时的文件
###
数据整合功能不足处：细胞量较少的数据（约少于200个细胞）在使用某种整合方式时会出错（不过seurat提供4种整合方式，总有合适的）
###
虽然脚本功能简单，也还有很多值得改进的地方和能增加的下游分析功能，暂时先这样吧,做为一个阶段的学习成果收尾 ^_^
### end
