"""
数学建模(P154)  排队模型: 港口系统
"""

import numpy as np

global harTime                                 # 每艘船待在港口的平均时间
global maxHar                                  # 一艘船待在港口的最长时间
global waitTime                                # 每艘船卸货之前的平均等待时间
global maxWait                                 # 一艘船卸货之前的最长等待时间
global idleTime                                # 卸货设备空闲时间占总模拟时间的百分比

def portSimulation(n):
    """
    港口系统模拟算法
    n                   模拟中的船只总数
    returns
        harTime         每艘船待在港口的平均时间
        maxHar          一艘船待在港口的最长时间
        waitTime        每艘船卸货之前的平均等待时间
        maxWait         一艘船卸货之前的最长等待时间
        idleTime        卸货设备空闲时间占总模拟时间的百分比
    """
    between = np.zeros(n)                   # 相邻两船到达港口的时间间隔数组
    unload = np.zeros(n)                    # 每艘船的卸货时间
    arrive = np.zeros(n)                    # 每艘船到达港口的时间
    start = np.zeros(n)                     # 每艘船开始卸货的时间
    wait = np.zeros(n)                      # 每艘船到达后到开始卸货的等待时间
    finish = np.zeros(n)                    # 每艘船卸货完毕的时间
    harbor = np.zeros(n)                    # 每艘船待在港口总的时间

    BL, BH = 15, 145
    UL, UH = 45, 90
    # 初始化变量
    between[0] = np.random.randint(low = BL, high = BH)
    unload[0] = np.random.randint(low = UL, high = UH)
    harTime = unload[0]
    maxHar = unload[0]
    waitTime = 0
    maxWait = 0
    idleTime = between[0]
    # 船1
    arrive[0] = between[0]
    start[0] = arrive[0]
    finish[0] = between[0] + unload[0]
    # 船2:n
    for i in range(1, n):
        between[i] = np.random.randint(low = BL, high = BH)
        unload[i] = np.random.randint(low = UL, high = UH)
        # 计算船i到达的时间
        arrive[i] = arrive[i - 1] + between[i]
        # 计算船i与船i-1卸货完毕的时间之差
        timeDiff = arrive[i] - finish[i - 1]
        # 判断是否等待
        if timeDiff > 0:
            idleTime += timeDiff
        else:
            wait[i] = - timeDiff
        # 计算船i开始卸货的时间
        start[i] = arrive[i] + wait[i]
        # 计算船i卸货完毕的时间
        finish[i] = start[i] + unload[i]
        # 计算船i待在港口的时间
        harbor[i] = finish[i] - arrive[i]
        # 总的港口时间
        harTime += harbor[i]
        # 最大港口时间
        if maxHar < harbor[i]:
            maxHar = harbor[i]
        # 总的等待时间
        waitTime += wait[i]
        # 最大等待时间
        if maxWait < wait[i]:
            maxWait = wait[i]
    # 计算每艘船待在港口的平均时间
    harTime /= n
    # 计算每艘船卸货之前的平均等待时间
    waitTime /= n
    # 计算卸货设备空闲时间占总模拟时间的百分比
    idleTime /= finish[n - 1]

    return harTime, maxHar, waitTime, maxWait, idleTime



if __name__ == "__main__":
    i = 0
    while i < 6:
        print("%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n"%portSimulation(100))
        i += 1



"""
120.82  293.00  52.58   217.00  0.11

85.69   178.00  18.84   110.00  0.22

119.55  296.00  52.68   231.00  0.17

106.54  282.00  41.02   206.00  0.18

133.33  330.00  65.24   243.00  0.11

107.78  231.00  40.90   156.00  0.10
"""