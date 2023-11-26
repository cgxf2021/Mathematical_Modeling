from os import terminal_size
import numpy as np
import matplotlib.pyplot as plt

"""
遗传算法解决TSP问题
"""

def initPop(NIND, N):
    """
    初始化种群
    NIND: 种群大小
    N: 染色体长度
    returns: 初始化的种群
    """

    chrom = np.zeros(shape = (NIND, N), dtype = int)                         # 存储种群
    for i in range(NIND):
        chrom[i, :] = np.random.permutation(N)                  # 随机生成初始种群

    return chrom


def fitness(distance):
    """
    适应度函数
    distance: 个体TSP的距离矩阵
    returns: 个体适应度矩阵
    """

    return 1 / distance


def sus(FitnV, Nsel):
    """
    随机通用采样算法
    FitnV: 个体适应度矩阵
    Nsel: 被选择个体的数目
    returns: 被选择个体的索引号
    """

    pass











if __name__ == "__main__":
    pass
