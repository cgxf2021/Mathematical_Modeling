import numpy as np
import matplotlib.pyplot as plt
import geatpy as ea

# matplotlib字体
plt. rcParams['font.sans-serif']= ['SimHei']
plt. rcParams['axes.unicode_minus'] = False

# 求fun(x, y)的最大值   x \in [-2, 2], y \in [-2, 2]
# 定义函数fun(x, y) 
fun = lambda x, y: x * np.cos(2 * np.pi * y) + y * np.sin(2 * np.pi * x)

def reinsert(Chrom, SelCh, ObjV, ObjVSel):
    """
    重插入函数(父代精英重插入)
    Chrom: 父代种群
    SelCh: 子代种群
    ObjV: 父代目标函数值
    ObjVSel: 子代目标函数值
    return: 重插入后的种群
    """

    parent_shape = Chrom.shape                                      # 父代种群大小
    offsprint_shape = SelCh.shape                                   # 子代种群大小
    in_ObjV = np.zeros(1)                                           # 待插入的父代个体值
    in_Chrom = np.zeros(shape = (1, parent_shape[1]))               # 待插入的父代个体
    rein_count = parent_shape[0] - offsprint_shape[0]               # 重插入的个数
    rein_index = np.zeros(rein_count, np.int8)                      # 重插入索引

    # 选出父代精英
    for i in range(rein_count):
        # 获取待插入个体索引
        rein_index[i] = np.argmax(ObjV)
        # 待插入父代值插入子代值数组
        in_ObjV[0] = ObjV[rein_index[i]]
        ObjVSel = np.concatenate((ObjVSel, in_ObjV))
        # 待插入父代个体插入子代种群
        in_Chrom = Chrom[rein_index[i]].reshape(1, parent_shape[1])
        SelCh = np.concatenate((SelCh, in_Chrom), axis = 0)
        # 置最大值为-1
        ObjV[rein_index[i]] = -1

    return SelCh, ObjVSel


# 遗传算法参数设置
"""
种群大小            40
最大遗传代数        100
个体长度            20 + 20
代沟                0.95
交叉概率            0.7
变异概率            0.01
"""

# 画出fun(x, y)的图像
fig_1 = plt.figure()
ax_1 = fig_1.add_subplot(111, projection = '3d')
x = np.linspace(-2, 2, 300)
y = np.linspace(-2, 2, 300)
X, Y = np.meshgrid(x, y)
Z = fun(X, Y)
ax_1.plot_surface(X, Y, Z, cmap = 'cool')
ax_1.set_title(r"$ f(x, y) = xcos(2 \pi y) + ysin(2 \pi x) $")
ax_1.set_xlabel("$ x $")
ax_1.set_ylabel("$ y $")
ax_1.set_zlabel("$ z $")

# 定义遗传算法参数
NIND = 40                   # 种群大小
MAXGEN = 100                # 最大遗传代数
PRECI = 20                  # 个体长度
GGAP = 0.95                 # 代沟
PX = 0.7                    # 交叉概率
PM = 0.01                   # 变异概率
trace = np.zeros(shape = (3, MAXGEN))                                                               # 寻优结果的初始值
FieldD = np.array([[PRECI, PRECI], [-2, -2], [2, 2], [1, 1], [0, 0], [1, 1], [1, 1], [0, 0]])       # 区域描述器
chrom = ea.crtbp(NIND, PRECI * 2)                                                                   # 创建任意离散随机种群

# 优化
gen = 0                             # 代计数器
XY = ea.bs2real(chrom, FieldD)      # 初始种群的十进制转换
X = XY[:, 0]
Y = XY[:, 1]
ObjV = fun(X, Y)                    # 计算目标函数值

while gen < MAXGEN:
    FitnV = ea.ranking(ObjV.reshape(NIND, 1))               # 分配适应度值
    SelChIn = ea.selecting('sus', FitnV, GGAP)              # 选择
    SelCh = ea.recombin('xovsp', chrom[SelChIn], PX)        # 重组
    SelCh = ea.mutbin('BG', SelCh, None, PM)                # 变异
    XY = ea.bs2real(SelCh, FieldD)                          # 子代个体十进制转换
    X = XY[:, 0]
    Y = XY[:, 1]
    ObjVSel = fun(X, Y)                                     # 计算子代的目标函数值
    chrom, ObjV = reinsert(chrom, SelCh, ObjV, ObjVSel)     # 重插入
    XY = ea.bs2real(chrom, FieldD)
    # 获取每代的最优解以及索引
    best_V = np.amax(ObjV)
    best_I = np.argmax(ObjV)
    # 记下每代最优值
    trace[0, gen] = XY[best_I, 0]
    trace[1, gen] = XY[best_I, 1]
    trace[2, gen] = best_V
    gen = gen + 1                                           # 代数增加1

# 最优解
print("X = %f"%trace[0, -1])
print("Y = %f"%trace[1, -1])
print("Z = %f"%trace[2, -1])

# 画出每代的最优点
ax_1.scatter(xs = trace[0, -1], ys = trace[1, -1], zs = trace[2, -1], s = 30, c = 'r')
ax_1.view_init(elev = 30, azim = 60)

# 画出进化图
fig_2 = plt.figure()
ax_2 = fig_2.add_subplot(111)
ax_2.plot(trace[2, :])
ax_2.set_title("进化过程")
ax_2.set_xlabel("遗传代数")
ax_2.set_ylabel("解的变化")
plt.show()


"""
OutCome:
X = 1.762672
Y = -2.000000
Z = 3.756336
"""