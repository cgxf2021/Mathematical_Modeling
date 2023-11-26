import numpy as np
import matplotlib.pyplot as plt

# matplotlib字体
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

"""
求函数 y = x + 10 sin(5x) + 7 cos(4x) 的最大值
"""
# 定义函数
fun = lambda x: x + 10 * np.sin(5 * x) + 7 * np.cos(4 * x)

# 模拟退火判断函数
def judge(dE, T):
    if dE < 0:
        return True
    else:
        p = np.exp(-(dE / T))
        pick = np.random.uniform()
        while pick == 0:
            pick = np.random.uniform()
        if p > pick:
            return True
        else:
            return False



if __name__ == "__main__":
    # 设定参数
    section_l = 0                                       # 区间下限
    section_h = 9                                       # 区间上限
    tmp = 1e5                                           # 初始温度
    tmp_min = 1e-3                                      # 停止温度
    alpha = 0.90                                        # 降温系数

    # 画出目标函数的图像
    fig = plt.figure()
    ax = fig.add_subplot(111)
    x = np.linspace(0, 9, 200)
    ax.plot(x, fun(x), color = 'orange')
    ax.set_title(r"$ y = x + 10 sin5x + 7 cos4x $")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    
    # 生成初始随机解
    pick = np.random.uniform()
    x_old = (section_h - section_l) * pick + section_l
    s_old = fun(x_old)
    x_new, s_new = x_old, s_old

    # 初始化图片
    text = ax.text(0.3, 21, "tmp: %.4f"%tmp, fontsize = 10)

    # 计数器
    count = 0

    # 开始模拟退火
    while tmp > tmp_min:
        flag = 0
        # 随机扰动
        pick = np.random.uniform()
        delta = (pick - 0.5) * 3
        x_new = x_old + delta
        if x_new < section_l or x_new > section_h:
            x_new = x_new - 2 * delta
        s_new = fun(x_new)
        # 求差值dE
        dE = s_old - s_new
        # 判断
        if judge(dE, tmp) == True:
            s_old = s_new
            x_old = x_new
        # 只有dE < 0才降温
        if dE < 0:
            # # 作图
            ax.scatter(x_old, s_old, marker = 'o', cmap = 'rainbow')
            text.set_text("tmp: %.4f"%tmp)
            plt.pause(0.001)
            tmp = tmp * alpha
        else:
            count = count + 1
        # 当接受更差的解的概率太小时，若又长时间找不到更优的解，那么退出循环，结束算法
        if count > 10000:
            break
    
    # 打印结果
    print("maximum: ", x_old, "  ", s_old)
    # 作图
    plt.show()
