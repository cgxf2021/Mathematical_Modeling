import numpy as np
import matplotlib.pyplot as plt

def fun(x):
    """
    目标函数值
    x: 实值变量向量
    returns: 目标函数值
    """
    return - 5 * np.sin(x[0]) * np.sin(x[1]) * np.sin(x[2]) * np.sin(x[3]) * np.sin(x[4]) - np.sin(5 * x[0]) * np.sin(5 * x[1]) * np.sin(5 * x[2]) * np.sin(5 * x[3]) * np.sin(5 * x[4]) + 8


def fitness(x):
    """
    适应度函数
    x: 实值变量向量
    returns: 适应度
    """
    return 1 / fun(x)


def select(chrom, popsize, popfitness):
    """
    选择函数(采用轮赌盘法)
    chrom: 种群信息
    popsize: 种群规模
    popfitness: 种群适应度
    returns: 经过选择后的种群和种群适应度
    """

    sumfitness = np.sum(popfitness)             # 个体适应度总和
    fitness_ratio = popfitness / sumfitness     # 个体适应度占比
    new_chrom = np.zeros(shape = chrom.shape)   # 记录选择的新个体
    # 轮赌盘算法
    for i in range(popsize):
        # 获取一个(0, 1)随机数
        pick = np.random.rand()
        while pick == 0.0:
            pick = np.random.rand()
        for j in range(popsize):
            pick = pick - fitness_ratio[j]
            if pick < 0:
                # 记录被选择到的个体
                new_chrom[i, :] = chrom[j, :]
                # 计算被选个体的适应度
                popfitness[i] = fitness(new_chrom[i, :])
                break

    return new_chrom, popfitness


def cross(chrom, chrom_len, popsize, pcross, bound):
    """
    染色体交叉函数
    chrom: 染色体
    chrom_len: 染色体长度
    popsize: 种群规模
    pcross: 交叉概率
    bound: 边界
    returns: 交叉后的染色体
    """
    
    for i in range(popsize):
        # 随机选择两个染色体进行交叉
        cross_index = np.random.randint(low = 0, high = popsize, size = 2, dtype = int)
        # 交叉概率决定是否交叉
        pick = np.random.rand()
        while pick == 0:
            pick = np.random.rand()
        if pick > pcross:
            continue
        flag = 0
        while flag == 0:
            # 随机选择交叉位置
            pos = np.random.randint(low = 0, high = chrom_len, size = 1, dtype = int)
            # 交叉开始
            pick = np.random.rand()
            v1 = chrom[cross_index[0], pos]
            v2 = chrom[cross_index[1], pos]
            aj = pick * v2 + (1 - pick) * v1
            al = pick * v1 + (1 - pick) * v2
            # 检验染色体可行性
            flag1, flag2 = 1, 1
            if aj < bound[0] or aj > bound[1]:
                flag1 = 0
            if al < bound[0] or aj > bound[1]:
                flag2 = 0
            if flag1 * flag2 == 0:
                flag = 0
            else:
                flag = 1
                chrom[cross_index[0], pos] = aj
                chrom[cross_index[1], pos] = al
    
    return chrom


def mutation(chrom, chrom_len, popsize, pmutation, pop, bound):
    """
    染色体变异函数
    chrom: 染色体
    chrom_len: 染色体长度
    popsize: 种群规模
    pmutation: 突变概率
    pop: 当前种群的进化代数和最大进化代数
    bound: 边界
    returns: 变异后的染色体
    """

    for i in range(popsize):
        # 随机选择一个染色体进行变异
        chrom_index = np.random.randint(low = 0, high = popsize, size = 1, dtype = int)
        # 变异概率决定该轮是否进行变异
        pick = np.random.rand()
        while pick == 0:
            pick = np.random.rand()
        if pick > pmutation:
            continue
        flag = 0
        while flag == 0:
            # 选取变异位置
            pos = np.random.randint(low = 0, high = chrom_len, size = 1, dtype = int)
            v = chrom[i, pos]
            v1 = v - bound[0]
            v2 = bound[1] - v
            # 开始变异
            pick = np.random.rand()
            delta = 0
            if pick >= 0.5:
                delta = v2 * (1 - pick ** ((1 - pop[0] / pop[1]) ** 2))
                if bound[0] <= v + delta <= bound[1]:
                    chrom[i, pos] = v + delta
                    flag = 1
            else:
                delta = v1 * (1 - pick ** ((1 - pop[0] / pop[1]) ** 2))
                if bound[0] <= v + delta <= bound[1]:
                    chrom[i, pos] = v + delta
                    flag = 1

    return chrom


if __name__ == "__main__":
    # 遗传算法参数
    gen = 100                                               # 进化代数
    popsize = 100                                           # 种群规模
    cross_p = 0.6                                           # 交叉概率
    mutation_p = 0.01                                       # 突变概率
    chrom_len = 5                                           # 染色体长度
    bound = np.array([0, 0.9 * np.pi])                      # 边界

    # 个体初始化
    popfitness = np.zeros(popsize)                          # 每代种群适应度
    avgfitness = np.zeros(gen + 1)                          # 每代种群平均适应度
    genfitness = np.zeros(gen + 1)                          # 每代种群最优适应度
    bestchrom = np.zeros(chrom_len)                         # 最优染色体
    bestfitness = -100                                      # 最优适应度
    # 种群初始化
    chrom = np.random.uniform(low = bound[0], high = bound[1], size = (popsize, 5))
    for i in range(popsize):
        x = chrom[i, :]
        popfitness[i] = fitness(x)
    # 找到最好的染色体
    bestindex = np.argmax(popfitness)
    genfitness[0] = popfitness[bestindex]
    # 计算平均适应度
    avgfitness[0] = np.sum(popfitness) / popsize

    # 开始进化
    for i in range(1, gen + 1):
        # 选择
        chrom, popfitness = select(chrom, popsize, popfitness)
        # 交叉
        chrom = cross(chrom, chrom_len, popsize, cross_p, bound)
        # 变异
        pop = np.array([i, gen])
        chrom = mutation(chrom, chrom_len, popsize, mutation_p, pop, bound)
        # 计算适应度
        for j in range(popsize):
            x = chrom[j, :]
            popfitness[j] = fitness(x)
        # 记录最优适应度
        bestindex = np.argmax(popfitness)
        genfitness[i] = popfitness[bestindex]
        # 记录最劣适应度
        worestindex = np.argmin(popfitness)
        # 记录最优染色体
        if bestfitness < genfitness[i]:
            bestfitness = genfitness[i]
            bestchrom = chrom[bestindex, :]
        # 重插入
        chrom[worestindex, :] = bestchrom[:]
        popfitness[worestindex] = bestfitness
        # 计算平均适应度
        avgfitness[i] = np.sum(popfitness) / popsize
        
    
    # 打印结果
    print("minimum: ", 1 / bestfitness)
    print("vector: ", bestchrom)
    # 作图
    fig = plt.figure()
    ax = fig.add_subplot(111)
    x = np.arange(0, gen + 1)
    ax.plot(x, 1 / genfitness, color = 'blue', label = "1 / best fitness")
    ax.plot(x, 1 / avgfitness, color = 'orange', label = "1 / average fitness")
    ax.legend()
    ax.grid()
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 10)
    ax.set_title("The course of evolution")
    ax.set_ylabel("value")
    ax.set_xlabel("generation")
    plt.show()