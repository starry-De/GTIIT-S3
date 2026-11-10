# Week 7 核心知识点大总结：大型线性方程组迭代解法、收敛性理论与范数分析

---

## 模块一：为什么需要迭代法？（问题背景与复杂度分析）

在工程实际问题（如机器人关节静力学分析、桁架结构平衡受力、复杂电路节点电压求解、传热与连续介质力学的差分方程离散）中，系统最终都会归结为形如 $Ax = b$ 的线性方程组：

$$\begin{cases}
a_{11}x_1 + a_{12}x_2 + \cdots + a_{1n}x_n = b_1 \\
a_{21}x_1 + a_{22}x_2 + \cdots + a_{2n}x_n = b_2 \\
\quad \vdots \\
a_{n1}x_1 + a_{n2}x_2 + \cdots + a_{nn}x_n = b_n
\end{cases}$$

### 1. 高斯消元法（直接法）的瓶颈
* **时间开销（浮点运算量）**：对于 $n$ 阶稠密矩阵，高斯消元法（LU 分解）的运算次数约为：
  $$\sum_{k=1}^n (n - k)^2 \approx \frac{n^3}{3} = \mathcal{O}(n^3)$$
  若系统规模达到 $n = 10^6$（百万阶），即使计算机每秒能执行 $10^9$ 次浮点运算（1 GFLOPS），纯计算时间也需要：
  $$t \approx \frac{1}{3} \times \frac{(10^6)^3}{10^9} \text{ 秒} \approx 3.33 \times 10^8 \text{ 秒} \approx 10.56 \text{ 年}$$
* **内存开销（稀疏性破坏，Fill-in 现象）**：工程中的大规模矩阵绝大多数是**稀疏矩阵（Sparse Matrix）**（例如每个节点只与其相邻的少数几个节点相连，99.9% 的元素为 0）。高斯消元过程中的行变换会使原本为 0 的位置迅速被非零数填满（Fill-in），直接导致内存耗尽。

### 2. 迭代法（Iterative Methods）的破局哲学
迭代法将单变量非线性方程的“不动点迭代思想”推广至多维向量空间：
* 从一个任意初始猜测向量 $x^{(0)}$（通常取全 0 向量）出发；
* 每一轮迭代（称为一次 **Sweep**）直接利用矩阵原有的非零系数做矩阵-向量乘法；
* 每次 Sweep 的计算代价仅正比于非零元的个数（对稀疏矩阵通常只需 $\mathcal{O}(n)$ 次运算）；
* 随着迭代步数 $k \to \infty$，使得 $x^{(k)} \to x$。在实际工程中，通常迭代十几步到几十步即可达到所需的机器精度。

---

## 模块二：矩阵的加法拆分（Additive Splitting）深度推导

要建立向量的不动点迭代格式 $x^{(k)} = T x^{(k-1)} + c$，关键在于将矩阵 $A$ 拆解。

### 1. 拆分定义：$A = D - L - U$
将 $n \times n$ 阶方阵 $A$ 严格拆解为三个子矩阵的加减组合：
$$A = D - L - U$$

其中：
* **$D$（对角矩阵，Diagonal）**：
  $$D = \begin{bmatrix}
  a_{11} & 0 & \dots & 0 \\
  0 & a_{22} & \dots & 0 \\
  \vdots & \vdots & \ddots & \vdots \\
  0 & 0 & \dots & a_{nn}
  \end{bmatrix}$$
  *前提条件*：主对角元 $a_{ii} \ne 0$。其逆矩阵为对角线元素直接取倒数：
  $$D^{-1} = \text{diag}\left(\frac{1}{a_{11}}, \frac{1}{a_{22}}, \dots, \frac{1}{a_{nn}}\right)$$

* **$-L$（严格下三角部分）与 $-U$（严格上三角部分）**：
  为了使 $A = D - L - U$ 成立，$L$ 和 $U$ 在数学上定义为原矩阵非对角线元素的**相反数**：
  $$L_{ij} = \begin{cases} -a_{ij}, & i > j \\ 0, & i \le j \end{cases}, \qquad U_{ij} = \begin{cases} -a_{ij}, & i < j \\ 0, & i \ge j \end{cases}$$

### 2. 为什么定义中带有负号？（数学优雅性与代数移项）
考虑线性方程组第 $i$ 行的分量形式：
$$a_{ii}x_i + \sum_{j \ne i} a_{ij}x_j = b_i$$
为了求出当前未知数 $x_i$，必须将除 $x_i$ 之外的所有项移到等号右边：
$$a_{ii}x_i = b_i \mathbf{-} \sum_{j \ne i} a_{ij}x_j$$
注意求和符号前面的**负号**。
* 若直接拆分成全正号 $A = D + L_{\text{原}} + U_{\text{原}}$，则移项后矩阵公式为：
  $$Dx = -(L_{\text{原}} + U_{\text{原}})x + b \implies x = \mathbf{-}D^{-1}(L_{\text{原}} + U_{\text{原}})x + D^{-1}b$$
  每一步都必须拖带一个负号。
* 若采用 $A = D - L - U$，则：
  $$(D - L - U)x = b \implies Dx = \mathbf{(L + U)}x + b$$
  **负负得正**，移项产生的负号被矩阵定义提前吸收，使得迭代矩阵形式 $T = D^{-1}(L + U)$ 极其清爽整洁。

---

## 模块三：两大经典迭代算法详述

### 1. 雅可比迭代法（Jacobi Method, Algorithm 7.1）

#### 核心直觉
**“分居独立，全用旧值”**。在计算第 $k$ 轮的每一个新未知数 $x_i^{(k)}$ 时，等号右侧所有涉及其他变量的值，一律使用上一轮（第 $k-1$ 轮）计算完毕的旧值。

#### 分量迭代形式
$$x_i^{(k)} = \frac{1}{a_{ii}}\left(b_i - \sum_{j \ne i} a_{ij} x_j^{(k-1)}\right), \quad i = 1, 2, \dots, n$$

#### 矩阵迭代形式
由 $Dx = (L + U)x + b$，引入迭代步角标：
$$D x^{(k)} = (L + U)x^{(k-1)} + b$$
两边左乘 $D^{-1}$：
$$x^{(k)} = \underbrace{D^{-1}(L + U)}_{T_{jac}} x^{(k-1)} + \underbrace{D^{-1}b}_{c_{jac}}$$
* **Jacobi 迭代矩阵**：$T_{jac} = D^{-1}(L + U)$
* **常数项**：$c_{jac} = D^{-1}b$

#### 几何轨迹与工程特点
* **几何表现**：在二维平面上，由于两个分量不能同时获知对方的最新修正量，更新轨迹在真实交点两侧呈现相互震荡的 **“Z 字形”折线（Zig-zag）**。
* **工程优势**：在计算 $x_1^{(k)}, x_2^{(k)}, \dots, x_n^{(k)}$ 时，各分量彼此互不依赖，可以同时分发到多核 CPU 或 GPU 的成千上万个流处理器上进行**完全并行计算**。
* **工程代价**：计算机必须在内存中开辟两个独立数组（一个存放上一轮旧向量 $x^{(k-1)}$，一个存放正在写入的新向量 $x^{(k)}$）。

---

### 2. 高斯-赛德尔迭代法（Gauss–Seidel Method, Algorithm 7.2）

#### 核心直觉
**“贪心思想，见新就用”**。当计算机按次序依次计算 $x_1, x_2, \dots, x_n$ 时，计算到 $x_i^{(k)}$ 时，编号排在前面的分量 $x_1^{(k)}, \dots, x_{i-1}^{(k)}$ 已经更新完毕，它们显然比旧值更接近真解，因此立刻使用最新值参与计算。

#### 分量迭代形式
$$x_i^{(k)} = \frac{1}{a_{ii}}\left(b_i - \sum_{j < i} a_{ij} x_j^{(k)} - \sum_{j > i} a_{ij} x_j^{(k-1)}\right), \quad i = 1, 2, \dots, n$$

#### 矩阵迭代形式推导
将包含最新值 $x^{(k)}$ 的对角项与下三角项全部保留在等号左边：
$$a_{ii}x_i^{(k)} + \sum_{j < i} a_{ij}x_j^{(k)} = b_i - \sum_{j > i} a_{ij}x_j^{(k-1)}$$
写成矩阵块：
* 左侧对应矩阵 $D - L$ 作用于新向量 $x^{(k)}$；
* 右侧对应上三角矩阵 $U$ 作用于旧向量 $x^{(k-1)}$ 加上常数项 $b$。

因此得到矩阵隐式关系：
$$(D - L)x^{(k)} = Ux^{(k-1)} + b$$
两边同时左乘逆矩阵 $(D - L)^{-1}$：
$$x^{(k)} = \underbrace{(D - L)^{-1}U}_{T_{gs}} x^{(k-1)} + \underbrace{(D - L)^{-1}b}_{c_{gs}}$$
* **Gauss–Seidel 迭代矩阵**：$T_{gs} = (D - L)^{-1}U$
* **常数项**：$c_{gs} = (D - L)^{-1}b$

#### 几何轨迹与工程特点
* **几何表现**：由于分量更新具有序贯性，一个方向更新完后立刻引导下一个方向，轨迹在平面上呈现贴近真实解的 **“台阶形”折线（Staircase）**。
* **工程优势**：
  1. **内存极省**：只需要维护一个数组，计算出 $x_i^{(k)}$ 后直接原地覆盖（In-place overwrite）写入内存。
  2. **收敛更快**：通常情况下，高斯-赛德尔法的收敛速度大约是雅可比法的 2 倍（即做 1 步 Gauss–Seidel 的误差缩减量相当于做 2 步 Jacobi）。
* **工程劣势**：分量计算存在严格的前后依赖（数据冒险），无法像雅可比法那样进行大规模无锁并行加速。

---

## 模块四：收敛性核心理论（特征值、谱半径与收敛充要定理）

并不是所有的方程组拿来迭代都能收敛。若方程排列不当，更新过程会迅速发散（Run away）。判断收敛性的终极判据是**谱半径**。

### 1. 特征值与谱半径定义
设 $M$ 为 $n \times n$ 阶实方阵：
* **特征值（Eigenvalue）**：满足 $M v = \lambda v \; (v \ne 0)$ 的标量 $\lambda$，其代数求法为特征方程：
  $$\det(M - \lambda I) = 0$$
* **谱半径（Spectral Radius, $\rho(M)$）**：矩阵所有特征值的**模长（复数绝对值）中的最大值**：
  $$\mathbf{\rho(M) = \max_{1 \le i \le n} |\lambda_i|}$$

### 2. 线性迭代法收敛充要定理
对于形如 $x^{(k)} = T x^{(k-1)} + c$ 的线性定常迭代系统，设精确解满足 $x = T x + c$。两式相减得到误差向量的递推演化：
$$e^{(k)} = x^{(k)} - x = T(x^{(k-1)} - x) = T e^{(k-1)} = T^k e^{(0)}$$

**【定理：线性迭代法收敛充要定理】**
**序列 $\{x^{(k)}\}$ 对任意初始猜测向量 $x^{(0)}$ 都能收敛到方程组唯一解的充分必要条件是：迭代矩阵 $T$ 的谱半径严格小于 1，即：**
$$\mathbf{\rho(T) < 1}$$

* 若 $\rho(T) < 1$：$T^k \to 0$ 当 $k \to \infty$，误差必定衰减至 0。
* 若 $\rho(T) \ge 1$：存在某些初始误差方向，使得 $T^k e^{(0)}$ 保持震荡或指数级放大（发散）。
* **渐近收敛速率**：$\rho(T)$ 越接近 0，每轮迭代的误差缩减比例越大，收敛越快。

### 3. 充分条件：严格对角占优准则
在实际手算或编程前，求大矩阵的特征值非常繁琐。数值分析给出了一个无需计算特征值即可直接保证收敛的几何充分条件：

**【定理：严格对角占优矩阵收敛定理】**
**若系数矩阵 $A$ 满足严格行对角占优（Strictly Diagonally Dominant），即每一行主对角元绝对值严格大于该行其余所有元素绝对值之和：**
$$\mathbf{|a_{ii}| > \sum_{j \ne i} |a_{ij}|, \quad \forall i = 1, 2, \dots, n}$$
**则：**
1. **矩阵 $A$ 严格非奇异（可逆）；**
2. **雅可比迭代法与高斯-赛德尔迭代法对任意初值 $x^{(0)}$ 均保证绝对收敛！**

*工程启示*：在编写算法前，应优先通过交换行（Pivoting），确保主对角线上的系数占据绝对压倒优势。

---

## 模块五：停机准则与向量范数体系

由于无法预先知道真实解 $x$，算法终止必须依赖相邻两步迭代解的变化量。

### 1. 向量范数（Vector Norms）定义
用于度量多维向量 $x = [x_1, x_2, \dots, x_n]^T$ 的大小/长度：

| 范数名称 | 符号表示 | 数学定义 | 直观几何意义 |
| :--- | :--- | :--- | :--- |
| **1-范数** | $\|x\|_1$ | $\sum_{i=1}^n \|x_i\|$ | 曼哈顿距离（各分量绝对值之和） |
| **2-范数** | $\|x\|_2$ | $\sqrt{\sum_{i=1}^n x_i^2}$ | 欧几里得距离（多维空间真实直线几何长度） |
| **$\infty$-范数** | $\|x\|_\infty$ | $\max_{1 \le i \le n} \|x_i\|$ | 最大值范数（所有分量中绝对值最大的那一项） |

### 2. 相对变化停机准则（Relative Change Stopping Criterion）
课件中推荐的工业级停止判据使用无穷范数构建：
$$\frac{\|x^{(k)} - x^{(k-1)}\|_\infty}{\|x^{(k)}\|_\infty} < TOL$$
* **分子**：本轮迭代中，波动最剧烈的单一分量的绝对变化幅度；
* **分母**：本轮解向量中最大的分量绝对值（做归一化，使得判断标准不受物理量纲或数值缩放的影响）；
* 当该相对变化率首次低于预设容差 $TOL$（如 $10^{-3}$）时，判定算法收敛并输出解。

---

## 模块六：经典作业与考试题型深度解析

---

### 【例题 1】三阶系统矩阵迭代形式、谱半径与收敛判定全解

**题目描述**：
设矩阵 $A = \begin{bmatrix} 3 & 2 & -2 \\ -2 & -2 & 1 \\ 5 & -5 & 4 \end{bmatrix}$：
(a) 计算迭代矩阵 $T_{jac}$ 与 $T_{gs}$。
(b) 计算各迭代矩阵的谱半径 $\rho(T_{jac})$ 与 $\rho(T_{gs})$。
(c) 判定 Jacobi 与 Gauss–Seidel 迭代是否对任意初始向量 $x^{(0)}$ 都能收敛，并说明理论依据。

---

#### 深度解析与演算步骤

#### (a) 计算 $T_{jac}$ 与 $T_{gs}$
1. **加法拆分矩阵**：
   $$D = \begin{bmatrix} 3 & 0 & 0 \\ 0 & -2 & 0 \\ 0 & 0 & 4 \end{bmatrix}, \quad 
   D^{-1} = \begin{bmatrix} \frac{1}{3} & 0 & 0 \\ 0 & -\frac{1}{2} & 0 \\ 0 & 0 & \frac{1}{4} \end{bmatrix}$$
   严格上下三角定义：
   $$L = \begin{bmatrix} 0 & 0 & 0 \\ 2 & 0 & 0 \\ -5 & 5 & 0 \end{bmatrix}, \quad 
   U = \begin{bmatrix} 0 & -2 & 2 \\ 0 & 0 & -1 \\ 0 & 0 & 0 \end{bmatrix}$$
   合成非对角矩阵：
   $$L + U = \begin{bmatrix} 0 & -2 & 2 \\ 2 & 0 & -1 \\ -5 & 5 & 0 \end{bmatrix}, \quad 
   D - L = \begin{bmatrix} 3 & 0 & 0 \\ -2 & -2 & 0 \\ 5 & -5 & 4 \end{bmatrix}$$

2. **计算 $T_{jac} = D^{-1}(L + U)$**：
   对角阵乘法即把 $(L+U)$ 各行分别乘以对角元：
   $$T_{jac} = \begin{bmatrix} \frac{1}{3} & 0 & 0 \\ 0 & -\frac{1}{2} & 0 \\ 0 & 0 & \frac{1}{4} \end{bmatrix} \begin{bmatrix} 0 & -2 & 2 \\ 2 & 0 & -1 \\ -5 & 5 & 0 \end{bmatrix} = \begin{bmatrix} 0 & -\frac{2}{3} & \frac{2}{3} \\ -1 & 0 & \frac{1}{2} \\ -\frac{5}{4} & \frac{5}{4} & 0 \end{bmatrix}$$

3. **计算 $T_{gs} = (D - L)^{-1}U$**：
   设未知矩阵 $X = T_{gs}$，转化为解下三角方程 $(D - L)X = U$：
   $$\begin{bmatrix} 3 & 0 & 0 \\ -2 & -2 & 0 \\ 5 & -5 & 4 \end{bmatrix} \begin{bmatrix} x_{11} & x_{12} & x_{13} \\ x_{21} & x_{22} & x_{23} \\ x_{31} & x_{32} & x_{33} \end{bmatrix} = \begin{bmatrix} 0 & -2 & 2 \\ 0 & 0 & -1 \\ 0 & 0 & 0 \end{bmatrix}$$
   由左至右逐列进行前向代入求解：
   * **求解第 1 列**：右边为 $[0, 0, 0]^T$，齐次方程直接解得：
     $$x_{11} = 0, \quad x_{21} = 0, \quad x_{31} = 0$$
   * **求解第 2 列**：右边为 $[-2, 0, 0]^T$：
     * 行 1：$3x_{12} = -2 \implies x_{12} = -\frac{2}{3}$
     * 行 2：$-2(-\frac{2}{3}) - 2x_{22} = 0 \implies \frac{4}{3} - 2x_{22} = 0 \implies x_{22} = \frac{2}{3}$
     * 行 3：$5(-\frac{2}{3}) - 5(\frac{2}{3}) + 4x_{32} = 0 \implies -\frac{20}{3} + 4x_{32} = 0 \implies x_{32} = \frac{5}{3}$
   * **求解第 3 列**：右边为 $[2, -1, 0]^T$：
     * 行 1：$3x_{13} = 2 \implies x_{13} = \frac{2}{3}$
     * 行 2：$-2(\frac{2}{3}) - 2x_{23} = -1 \implies -\frac{4}{3} - 2x_{23} = -1 \implies -2x_{23} = \frac{1}{3} \implies x_{23} = -\frac{1}{6}$
     * 行 3：$5(\frac{2}{3}) - 5(-\frac{1}{6}) + 4x_{33} = 0 \implies \frac{20}{6} + \frac{5}{6} + 4x_{33} = 0 \implies 4x_{33} = -\frac{25}{6} \implies x_{33} = -\frac{25}{24}$

   拼装矩阵：
   $$T_{gs} = \begin{bmatrix} 0 & -\frac{2}{3} & \frac{2}{3} \\ 0 & \frac{2}{3} & -\frac{1}{6} \\ 0 & \frac{5}{3} & -\frac{25}{24} \end{bmatrix}$$

---

#### (b) 计算谱半径
1. **Jacobi 矩阵的谱半径**：
   计算特征方程 $\det(T_{jac} - \lambda I) = 0$：
   $$\det \begin{bmatrix} -\lambda & -\frac{2}{3} & \frac{2}{3} \\ -1 & -\lambda & \frac{1}{2} \\ -\frac{5}{4} & \frac{5}{4} & -\lambda \end{bmatrix} = -\lambda\left(\lambda^2 - \frac{5}{8}\right) + \frac{2}{3}\left(\lambda + \frac{5}{8}\right) + \frac{2}{3}\left(-\frac{5}{4} - \frac{5}{4}\lambda\right) = 0$$
   展开合并同类项：
   $$-\lambda^3 + \left(\frac{5}{8} + \frac{2}{3} - \frac{5}{6}\right)\lambda + \left(\frac{5}{12} - \frac{5}{6}\right) = -\lambda^3 + \frac{11}{24}\lambda - \frac{10}{24} = 0$$
   即：
   $$24\lambda^3 - 11\lambda + 10 = 0$$
   经数值求解，该方程有 1 个实根和 2 个共轭复根：
   $$\lambda_1 \approx -0.898, \quad \lambda_{2,3} \approx 0.449 \pm 0.523i$$
   复根模长为 $|\lambda_{2,3}| = \sqrt{0.449^2 + 0.523^2} \approx 0.689$。
   因此：
   $$\rho(T_{jac}) = \max\{|-0.898|, 0.689\} \approx \mathbf{0.898}$$

2. **Gauss–Seidel 矩阵的谱半径**：
   观察 $T_{gs}$ 第 1 列全为 0，按第 1 列降阶展开特征方程：
   $$\det(T_{gs} - \lambda I) = -\lambda \det \begin{bmatrix} \frac{2}{3} - \lambda & -\frac{1}{6} \\ \frac{5}{3} & -\frac{25}{24} - \lambda \end{bmatrix} = 0$$
   得到一个显然特征值 $\lambda_1 = 0$。
   其余两根由行列式二阶子式确定：
   $$\left(\lambda - \frac{2}{3}\right)\left(\lambda + \frac{25}{24}\right) - \left(-\frac{1}{6}\right)\left(\frac{5}{3}\right) = 0$$
   $$\lambda^2 + \left(\frac{25}{24} - \frac{16}{24}\right)\lambda - \frac{50}{72} + \frac{20}{72} = 0 \implies \lambda^2 + \frac{3}{8}\lambda - \frac{5}{12} = 0$$
   同乘 24：
   $$24\lambda^2 + 9\lambda - 10 = 0$$
   求根公式：
   $$\lambda = \frac{-9 \pm \sqrt{9^2 - 4(24)(-10)}}{2 \times 24} = \frac{-9 \pm \sqrt{81 + 960}}{48} = \frac{-9 \pm \sqrt{1041}}{48}$$
   由于 $\sqrt{1041} \approx 32.264$：
   $$\lambda_2 \approx \frac{-9 + 32.264}{48} \approx 0.485, \quad \lambda_3 \approx \frac{-9 - 32.264}{48} \approx -0.860$$
   比较模长：
   $$\rho(T_{gs}) = \max\{|0|, |0.485|, |-0.860|\} = \frac{9 + \sqrt{1041}}{48} \approx \mathbf{0.860}$$

---

#### (c) 收敛性判断与解释
* **依据**：线性迭代法对任意初始向量 $x^{(0)}$ 均收敛的充分必要条件为迭代矩阵的谱半径严格小于 1（$\rho(T) < 1$）。
* **Jacobi 结论**：由于 $\rho(T_{jac}) \approx 0.898 < 1$，**Jacobi 迭代法对任意初值必定收敛**。
* **Gauss–Seidel 结论**：由于 $\rho(T_{gs}) \approx 0.860 < 1$，**Gauss–Seidel 迭代法对任意初值必定收敛**。
* **速度比较**：因 $\rho(T_{gs}) < \rho(T_{jac}) < 1$，Gauss–Seidel 法渐近收敛速度快于 Jacobi 法。

---

### 【例题 2】松弛迭代法（Richardson Iteration）参数分析与手算

**题目描述**：
求解方程组 $Ax = b$，其中 $A = \begin{bmatrix} 3 & 2 \\ 1 & 2 \end{bmatrix}, \; b = \begin{bmatrix} 3 \\ -1 \end{bmatrix}$。采用迭代格式：
$$x^{(k+1)} = x^{(k)} + \omega(Ax^{(k)} - b)$$
(a) 求参数 $\omega$ 满足何种取值范围时，该公式对任意初始向量均收敛？
(b) 当 $\omega = -\frac{1}{4}$，从初值 $x^{(0)} = \begin{bmatrix} 0 \\ 0 \end{bmatrix}$ 出发，手算执行前两次迭代。

---

#### 深度解析与演算步骤

#### (a) 求解 $\omega$ 的收敛范围
1. **提取迭代矩阵**：
   整理公式：
   $$x^{(k+1)} = x^{(k)} + \omega A x^{(k)} - \omega b = (I + \omega A)x^{(k)} - \omega b$$
   迭代矩阵为 $T(\omega) = I + \omega A$。

2. **计算 $A$ 的特征值**：
   $$\det(A - \lambda I) = \det \begin{bmatrix} 3-\lambda & 2 \\ 1 & 2-\lambda \end{bmatrix} = (3-\lambda)(2-\lambda) - 2 = \lambda^2 - 5\lambda + 4 = 0$$
   分解因式 $(\lambda - 1)(\lambda - 4) = 0$，解得：
   $$\lambda_1 = 1, \quad \lambda_2 = 4$$

3. **求迭代矩阵 $T(\omega)$ 的特征值 $\mu$**：
   由矩阵函数谱映射定理，若 $A$ 的特征值为 $\lambda$，则 $I + \omega A$ 的特征值为 $\mu = 1 + \omega \lambda$：
   $$\mu_1 = 1 + \omega, \quad \mu_2 = 1 + 4\omega$$

4. **强制条件 $\rho(T(\omega)) < 1$**：
   所有特征值必须落在 $(-1, 1)$ 内：
   $$\begin{cases}
   -1 < 1 + \omega < 1 \implies -2 < \omega < 0 \\
   -1 < 1 + 4\omega < 1 \implies -2 < 4\omega < 0 \implies -\frac{1}{2} < \omega < 0
   \end{cases}$$
   取交集得到最终收敛范围：
   $$\mathbf{-\frac{1}{2} < \omega < 0}$$

---

#### (b) $\omega = -\frac{1}{4}$ 时的两轮手算过程
迭代更新格式为：
$$x^{(k+1)} = x^{(k)} - \frac{1}{4}(Ax^{(k)} - b)$$

* **第 1 轮迭代（$k=0$ 求 $x^{(1)}$）**：
  将 $x^{(0)} = \begin{bmatrix} 0 \\ 0 \end{bmatrix}$ 代入残差式：
  $$Ax^{(0)} - b = \begin{bmatrix} 0 \\ 0 \end{bmatrix} - \begin{bmatrix} 3 \\ -1 \end{bmatrix} = \begin{bmatrix} -3 \\ 1 \end{bmatrix}$$
  代入更新：
  $$x^{(1)} = \begin{bmatrix} 0 \\ 0 \end{bmatrix} - \frac{1}{4}\begin{bmatrix} -3 \\ 1 \end{bmatrix} = \mathbf{\begin{bmatrix} \frac{3}{4} \\ -\frac{1}{4} \end{bmatrix}} = \begin{bmatrix} 0.75 \\ -0.25 \end{bmatrix}$$

* **第 2 轮迭代（$k=1$ 求 $x^{(2)}$）**：
  先计算 $Ax^{(1)}$：
  $$Ax^{(1)} = \begin{bmatrix} 3 & 2 \\ 1 & 2 \end{bmatrix} \begin{bmatrix} \frac{3}{4} \\ -\frac{1}{4} \end{bmatrix} = \begin{bmatrix} 3(\frac{3}{4}) + 2(-\frac{1}{4}) \\ 1(\frac{3}{4}) + 2(-\frac{1}{4}) \end{bmatrix} = \begin{bmatrix} \frac{7}{4} \\ \frac{1}{4} \end{bmatrix}$$
  计算残差 $Ax^{(1)} - b$：
  $$Ax^{(1)} - b = \begin{bmatrix} \frac{7}{4} \\ \frac{1}{4} \end{bmatrix} - \begin{bmatrix} 3 \\ -1 \end{bmatrix} = \begin{bmatrix} \frac{7}{4} - \frac{12}{4} \\ \frac{1}{4} + \frac{4}{4} \end{bmatrix} = \begin{bmatrix} -\frac{5}{4} \\ \frac{5}{4} \end{bmatrix}$$
  代入更新：
  $$x^{(2)} = x^{(1)} - \frac{1}{4}\begin{bmatrix} -\frac{5}{4} \\ \frac{5}{4} \end{bmatrix} = \begin{bmatrix} \frac{3}{4} \\ -\frac{1}{4} \end{bmatrix} + \begin{bmatrix} \frac{5}{16} \\ -\frac{5}{16} \end{bmatrix} = \mathbf{\begin{bmatrix} \frac{17}{16} \\ -\frac{9}{16} \end{bmatrix}} = \begin{bmatrix} 1.0625 \\ -0.5625 \end{bmatrix}$$

---

### 【例题 3】Jacobi 分量方程构建与双步手算

**题目描述**：
写出下列系统的 Jacobi 分量迭代方程，并从全零初值 $x^{(0)} = \mathbf{0}$ 开始执行两次迭代：
$$A = \begin{bmatrix} 4 & 2 & -1 \\ 2 & 4 & 1 \\ -1 & 1 & 4 \end{bmatrix}, \quad b = \begin{bmatrix} 1 \\ -1 \\ 1 \end{bmatrix}$$

---

#### 深度解析与演算步骤

1. **分量方程改写（主元保留，其余移项）**：
   $$\begin{cases}
   x_1^{(k)} = \frac{1}{4}\left(1 - 2x_2^{(k-1)} + x_3^{(k-1)}\right) \\[6pt]
   x_2^{(k)} = \frac{1}{4}\left(-1 - 2x_1^{(k-1)} - x_3^{(k-1)}\right) \\[6pt]
   x_3^{(k)} = \frac{1}{4}\left(1 + x_1^{(k-1)} - x_2^{(k-1)}\right)
   \end{cases}$$

2. **第 1 步迭代（将 $x_1^{(0)}=0, x_2^{(0)}=0, x_3^{(0)}=0$ 代入）**：
   $$x_1^{(1)} = \frac{1}{4}(1 - 0 + 0) = \mathbf{\frac{1}{4}}$$
   $$x_2^{(1)} = \frac{1}{4}(-1 - 0 - 0) = \mathbf{-\frac{1}{4}}$$
   $$x_3^{(1)} = \frac{1}{4}(1 + 0 - 0) = \mathbf{\frac{1}{4}}$$
   $$x^{(1)} = \begin{bmatrix} \frac{1}{4} \\ -\frac{1}{4} \\ \frac{1}{4} \end{bmatrix}$$

3. **第 2 步迭代（将 $x^{(1)}$ 全体带入右端）**：
   $$x_1^{(2)} = \frac{1}{4}\left[1 - 2\left(-\frac{1}{4}\right) + \frac{1}{4}\right] = \frac{1}{4}\left(1 + \frac{1}{2} + \frac{1}{4}\right) = \frac{1}{4}\left(\frac{7}{4}\right) = \mathbf{\frac{7}{16}} = 0.4375$$
   $$x_2^{(2)} = \frac{1}{4}\left[-1 - 2\left(\frac{1}{4}\right) - \frac{1}{4}\right] = \frac{1}{4}\left(-1 - \frac{1}{2} - \frac{1}{4}\right) = \frac{1}{4}\left(-\frac{7}{4}\right) = \mathbf{-\frac{7}{16}} = -0.4375$$
   $$x_3^{(2)} = \frac{1}{4}\left[1 + \frac{1}{4} - \left(-\frac{1}{4}\right)\right] = \frac{1}{4}\left(1 + \frac{1}{4} + \frac{1}{4}\right) = \frac{1}{4}\left(\frac{3}{2}\right) = \mathbf{\frac{3}{8}} = \frac{6}{16} = 0.375$$
   $$x^{(2)} = \begin{bmatrix} \frac{7}{16} \\ -\frac{7}{16} \\ \frac{6}{16} \end{bmatrix} = \begin{bmatrix} 0.4375 \\ -0.4375 \\ 0.375 \end{bmatrix}$$

---

### 【例题 4】向量三大范数计算

**题目描述**：
计算向量 $x = \begin{bmatrix} 2 \\ 1 \\ -3 \\ 4 \end{bmatrix}$ 的 $\|x\|_1$、$\|x\|_2$ 与 $\|x\|_\infty$。

---

#### 深度解析与演算步骤

1. **1-范数（绝对值求和）**：
   $$\|x\|_1 = |2| + |1| + |-3| + |4| = 2 + 1 + 3 + 4 = \mathbf{10}$$

2. **2-范数（分量平方和开根号）**：
   $$\|x\|_2 = \sqrt{2^2 + 1^2 + (-3)^2 + 4^2} = \sqrt{4 + 1 + 9 + 16} = \sqrt{30} \approx \mathbf{5.477}$$

3. **$\infty$-范数（绝对值最大分量）**：
   $$\|x\|_\infty = \max(|2|, |1|, |-3|, |4|) = \max(2, 1, 3, 4) = \mathbf{4}$$
```

这份文档整合了第 7 周课程涉及的全部理论机制与计算技巧。你可以查看右侧编辑器查阅和保存。对于其中的任何定理证明或计算细节，若有需要深入展开讨论的部分，随时告诉我！