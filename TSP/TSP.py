import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


class GeneticAlgorithm(object):
    def __init__(self) -> None:
        pass

    def initPop(self, NIND, N):
        """
        description: 初始化种群
        param: NIND 种群大小
        param: N 个体染色体长度(城市个数)
        Returns: Chrom 初始种群
        """

        Chrom = np.zeros(shape=(NIND, N), dtype=int)
        for i in range(NIND):
            Chrom[i] = np.random.permutation(N)

        return Chrom

    def distance_matrix(self, position):
        """
        description: 计算距离矩阵
        param: position 每个城市的坐标
        Returns: Distance 距离矩阵
        """

        L = len(position)
        Distance = np.zeros(shape=(L, L), dtype=float)
        for i in range(L):
            for j in range(L):
                Distance[i, j] = np.sqrt(np.power(
                    position[i, 0] - position[j, 0], 2) + np.power(position[i, 1] - position[j, 1], 2))

        return Distance

    def path_distance(self, Chrom, Distance):
        """
        description: 适应度函数
        param: Chrom 种群
        param: Distance 距离矩阵
        Returns: distance 个体距离向量
        """

        NIND, N = Chrom.shape
        distance = np.zeros(shape=(NIND, 1), dtype=float)
        for i in range(NIND):
            for j in range(N - 1):
                distance[i, 0] += Distance[Chrom[i, j], Chrom[i, j + 1]]
            distance[i, 0] += Distance[Chrom[i, j + 1], Chrom[i, 0]]

        return distance

    def fitness(self, distance):
        """
        description: 适应度函数
        param: distance 个体长度向量(TSP距离)
        Returns: FitnV 个体适应度向量
        """

        return 1 / distance

    def select(self, Chrom, FitnV, GGAP):
        """
        description: 选择操作
        param: Chrom 种群
        param: FitnV 适应度
        param: GGAP 选择概率
        Returns: SelCh 被选择的个体
        """

        def sus(Nsel, FitnV):
            """
            description: 
            param: FitnV 个体适应度向量
            param: Nsel 被选择个体的数目
            Returns: sel_index 选择到的个体索引
            """

            NIND = FitnV.shape[0]
            sel_index = np.zeros(Nsel, dtype=int)           # 保存被选择的个体下标
            cumfit = np.cumsum(FitnV)                       # 累计适应度
            cumfit_one = cumfit / cumfit[-1]                # 映射到[0, 1]
            for i in range(Nsel):
                rand = np.random.rand()
                for j in range(NIND):
                    if cumfit_one[j] > rand:
                        sel_index[i] = j
                        break

            return sel_index

        NIND = Chrom.shape[0]                               # 种群个数
        # 被选择的个体数
        Nsel = int(np.max(np.array([np.ceil(NIND * GGAP), 2])))
        sel_index = sus(Nsel=Nsel, FitnV=FitnV)
        SelCh = Chrom[sel_index, :]

        return SelCh

    def recombin(self, SelCh, Pc):
        """
        description: 交叉操作
        param: SelCh 被选中的个体
        param: Pc 交叉概率
        Returns: Selch 交叉后的个体
        """

        def intercross(A, B):
            """
            description: 两条染色体交叉
            param: A 染色体A
            param: B 染色体B
            Returns: 交叉后的两条染色体A, B 
            """

            L = len(A)
            # 随机生成两个交叉点
            r1, r2 = 0, 0
            while r1 == r2:
                r1 = np.random.randint(low=1, high=L - 1)
                r2 = np.random.randint(low=1, high=L - 1)
            if r1 > r2:
                r1, r2 = r2, r1
            # 交叉
            for i in range(r1, r2):
                A[i], B[i] = B[i], A[i]
            # 解决A冲突
            for i in range(r1):
                j = r1
                while j < r2:
                    if A[i] == A[j]:
                        A[i] = B[j]
                        j = r1
                        continue
                    j += 1
            for i in range(r2, L):
                j = r1
                while j < r2:
                    if A[i] == A[j]:
                        A[i] = B[j]
                        j = r1
                        continue
                    j += 1

            # 解决B冲突
            for i in range(r1):
                j = r1
                while j < r2:
                    if B[i] == B[j]:
                        B[i] = A[j]
                        j = r1
                        continue
                    j += 1
            for i in range(r2, L):
                j = r1
                while j < r2:
                    if B[i] == B[j]:
                        B[i] = A[j]
                        j = r1
                        continue
                    j += 1

            return A, B

        Nsel = SelCh.shape[0]
        for i in range(Nsel - 1):
            rand = np.random.rand()
            if rand < Pc:
                SelCh[i, :], SelCh[i + 1,
                                   :] = intercross(SelCh[i, :], SelCh[i + 1, :])

        return SelCh

    def mutate(self, SelCh, Pm):
        """
        description: 变异操作
        param: SelCh 被选择的个体
        param: Pm 变异概率
        Returns: SelCh 变异后的个体
        """

        Nsel, L = SelCh.shape
        for i in range(Nsel):
            rand = np.random.rand()
            if rand < Pm:
                r1, r2 = 0, 0
                while r1 == r2:
                    r1 = np.random.randint(low=0, high=L)
                    r2 = np.random.randint(low=0, high=L)
                if r1 > r2:
                    r1, r2 = r2, r1
                SelCh[i, r1], SelCh[i, r2] = SelCh[i, r2], SelCh[i, r1]

        return SelCh

    def reverse(self, SelCh, Distance):
        """
        description: 进化逆转操作
        param: SelCh 被选择的个体
        param: Distance 距离矩阵
        Returns: SelCh 进化逆转后的个体
        """

        row, L = SelCh.shape
        ObjV = self.path_distance(Chrom=SelCh, Distance=Distance)
        SelCh_1 = SelCh.copy()
        for i in range(row):
            r1, r2 = 0, 0
            while r1 == r2:
                r1 = np.random.randint(low=0, high=L)
                r2 = np.random.randint(low=0, high=L)
            if r1 > r2:
                r1, r2 = r2, r1
            SelCh_1[i, r1], SelCh_1[i, r2] = SelCh_1[i, r2], SelCh_1[i, r1]
        ObjV_1 = self.path_distance(Chrom=SelCh_1, Distance=Distance)
        for i in range(row):
            if ObjV[i] > ObjV_1[i]:
                SelCh[i, :] = SelCh_1[i, :]

        return SelCh

    def reins(self, Chrom, SelCh, ObjV, Distance, opt_path):
        """
        description: 父代精英重插入
        param: Chrom 父代种群
        param: SelCh 子代种群
        param: ObjV 父代适应度
        param: Distance 距离矩阵
        param: opt_path 父代精英
        Returns: SelCh 组合后得到的新种群
        """

        NIND = Chrom.shape[0]
        Nsel = SelCh.shape[0]
        index = np.argsort(ObjV)
        for i in range(NIND - Nsel):
            SelCh = np.vstack((SelCh, Chrom[index[i]]))
        ObjV_s = self.path_distance(SelCh, Distance)
        # 用父代精英替换最bad的个体
        bad_index = np.argmax(ObjV_s)
        SelCh[bad_index, :] = opt_path[:]
        return SelCh


class TSP(object):
    def __init__(self, pop_size, city_size, position, iters, Pc, Pm, Pa) -> None:
        self.pop_size = pop_size                # 种群大小
        self.city_size = city_size              # 城市个数
        self.position = position                # 城市坐标
        self.iters = iters                      # 迭代次数
        self.Pc = Pc                            # 交叉概率
        self.Pm = Pm                            # 变异概率
        self.Pa = Pa                            # 选择概率
        self.ga = GeneticAlgorithm()            # 遗传算法实例
        self.opt_distance = np.zeros(
            shape=(self.iters, 1), dtype=float)                # 记录每代最优个体值
        self.opt_path = np.zeros(
            shape=(1, self.city_size))                         # 父代精英
        self.all_opt_path = np.zeros(
            shape=(self.iters, self.city_size), dtype=int)     # 记录每代最优个体

    def compute_optimal_path(self):
        # 初始化种群
        self.Chrom = self.ga.initPop(self.pop_size, self.city_size)
        # 计算距离矩阵
        self.Distance = self.ga.distance_matrix(self.position)
        # 初始的距离
        self.ObjV = self.ga.path_distance(self.Chrom, self.Distance)
        # 初始时最优个体值
        self.preObjV = np.min(self.ObjV)
        # 初始时最优个体
        self.opt_path = self.Chrom[np.argmin(self.ObjV)]
        # 开始进化
        for gen in range(self.iters):
            # 计算适应度向量
            self.ObjV = self.ga.path_distance(self.Chrom, self.Distance)
            # 计算当前代适应度
            self.FitnV = self.ga.fitness(self.ObjV)
            # 选择
            self.SelCh = self.ga.select(self.Chrom, self.FitnV, self.Pa)
            # 交叉
            self.SelCh = self.ga.recombin(self.SelCh, self.Pc)
            # 变异
            self.SelCh = self.ga.mutate(self.SelCh, self.Pm)
            # 逆进化
            self.SelCh = self.ga.reverse(self.SelCh, self.Distance)
            # 找出最优个体(与父代精英比较)
            if np.min(self.ObjV) < self.preObjV:
                self.preObjV = np.min(self.ObjV)
                self.opt_path = self.Chrom[np.argmin(self.ObjV)]
            # 记录每代最优路径值
            self.opt_distance[gen, 0] = self.preObjV
            # 记录每代最优路径
            self.all_opt_path[gen, :] = self.opt_path
            # 重插入
            self.Chrom = self.ga.reins(
                self.Chrom, self.SelCh, self.ObjV, self.Distance, self.opt_path)

    def draw_evolution(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        x = np.arange(self.iters)
        ax.plot(x, self.opt_distance, label="optimal", color="orange")
        ax.set_xlabel("The number of iterations")
        ax.set_ylabel("The optimal distance of each generation")
        ax.set_title("The process of evolution")
        ax.legend()
        plt.savefig("evolution.pdf", bbox_inches="tight")
        plt.show()

    def draw_distance(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        def update(j):
            path = self.all_opt_path[j, :]
            plt.cla()
            ax.scatter(self.position[:, 0],
                       self.position[:, 1], marker="o", c="b")
            ax.set_title("gen = %d" %(j + 1))
            ax.text(23, 98, "spend: %.4f"%self.opt_distance[j, 0], fontsize = 10)
            for i in range(self.city_size):
                x_1, y_1 = self.position[path[i % self.city_size]]
                x_2, y_2 = self.position[path[(i + 1) % self.city_size]]
                dx, dy = x_2 - x_1, y_2 - y_1
                ax.arrow(x_1, y_1, dx, dy, head_width=0.1,
                         head_length=0.1, fc="r", ec="r")
        gens = np.hstack((np.arange(0, iters, 5), np.array([iters - 1])))
        ani = FuncAnimation(fig, update, frames=gens, interval=150, blit=False, repeat=False)
        ani.save("tsp.gif", writer='imagemagick')
        plt.show()


if __name__ == "__main__":
    pop_size = 20
    city_size = 14
    position = np.array([[16.47, 96.10],
                         [16.47, 94.44],
                         [20.09, 92.54],
                         [22.39, 93.37],
                         [25.23, 97.24],
                         [22.00, 96.05],
                         [20.47, 97.02],
                         [17.20, 96.29],
                         [16.30, 97.38],
                         [14.05, 98.12],
                         [16.53, 97.38],
                         [21.52, 95.59],
                         [19.41, 97.13],
                         [20.09, 92.55]])
    iters = 200
    Pc = 0.7
    Pm = 0.05
    Pa = 0.8
    tsp = TSP(pop_size=pop_size, city_size=city_size,
              position=position, iters=iters, Pc=Pc, Pm=Pm, Pa=Pa)
    tsp.compute_optimal_path()
    print("最优路径: ", tsp.opt_path)
    print("花费代价: ", tsp.preObjV)
    tsp.draw_evolution()
    tsp.draw_distance()

# 最优路径:  [ 8  9  0  1 13  2  3  4  5 11  6 12  7 10]
# 花费代价:  29.340520066994223