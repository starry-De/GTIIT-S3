# 热力学第一定律与物性评估（Thermodynamics I - Lecture 3）

---

## 模块一：状态原理与系统分类

### 1.1 什么是状态？
在热力学中，**状态（State）**是系统在某一瞬间所表现出的宏观物理状况。只要状态确定，该状态下的所有宏观热力学性质（如压力 $p$、温度 $T$、比容 $v$、内能 $u$、比焓 $h$ 等）都是唯一确定的。

### 1.2 状态原理（State Principle）
* **核心内容**：一个系统的稳定平衡状态，由它的**所有准静态功参数（Quasi-static work parameters）**和**能量（Energy, $E$）**唯一确定。
* **准静态功模式示例**：
  * **高空气球系统**：
    * 克服重力高度变化做功 $\rightarrow$ 功参数为高度 $h$
    * 膨胀做容积功 $\rightarrow$ 功参数为体积 $V$
    * 确定状态所需独立参数：$E$、$V$、$h$ 共 3 个

### 1.3 简单系统与简单可压缩系统
* **简单系统（Simple System）**：仅有 **1 种** 相关的准静态功作用模式的系统，只需 **2 个独立状态参数** 即可完全确定其状态。
* **简单可压缩系统（Simple Compressible System）**：
  * 工程热力学中最常见的工质模型（如活塞缸内的纯气体或液体）。
  * 唯一的准静态功模式为**容积胀缩功**（$p\text{d}V$ 功）。
  * **确定状态所需独立参数**：恰好 **2 个**（例如 $V$ 与 $E$，或强度性质组合 $p$ 与 $T$）。
  * 状态方程通式：
    $$f(p, T, v) = 0 \quad \text{（例如理想气体：} pv = RT \text{）}$$

> **注意：参数的独立性判定**
> 选取的两个参数必须相互独立。例如，对于定质量闭口系统，密度 $\rho$ 与比容 $v$ 存在 $\rho = 1/v$，两者互不独立，不能作为确定状态的两个独立变量。

---

## 模块二：基本热力学过程与物性参数（$u, h, c_v, c_p$）

由热力学第一定律闭口系能量方程：
$$\Delta U = Q - W$$

### 2.1 定容过程（Isochoric Process，$V = \text{常数}$）
* **物理过程**：容器刚性闭合或活塞被销钉（Pin）卡死，体积不变，容积做功 $W = 0$。
  $$\Delta U = Q$$
* **比内能（Specific Internal Energy, $u$）**：
  $u = u(T, v)$，对其全微分得：
  $$\text{d}u = \left(\frac{\partial u}{\partial T}\right)_v \text{d}T + \left(\frac{\partial u}{\partial v}\right)_T \text{d}v$$
* **定容比热容（Specific Heat at Constant Volume, $c_v$）**：
  * **定义式**：
    $$c_v = \left(\frac{\partial u}{\partial T}\right)_v \quad [\text{kJ/(kg}\cdot\text{K)}]$$
  * **物理意义**：单位质量工质在定容条件下升高 $1\text{ K}$ 所需吸收的热量。
  * **定容过程内能变化计算**：
    $$\Delta u = \int_{1}^{2} c_v \text{d}T$$

### 2.2 定压过程（Isobaric Process，$p = \text{常数}$）与焓的引入
* **物理过程**：活塞受恒定外压自由移动，系统吸热升温的同时发生膨胀对外做容积功 $W = p\Delta V$：
  $$\Delta U = Q - p\Delta V \implies Q = \Delta U + p\Delta V = \Delta(U + pV)$$
* **比焓（Specific Enthalpy, $h$）的定义**：
  为简化包含膨胀功和流动流动功的能量计算，定义新状态函数焓 $H$ 与比焓 $h$：
  $$H = U + pV \quad [\text{kJ}]$$
  $$h = u + pv \quad [\text{kJ/kg}]$$
* **定压比热容（Specific Heat at Constant Pressure, $c_p$）**：
  * 取 $h = h(T, p)$，对其全微分：
    $$\text{d}h = \left(\frac{\partial h}{\partial T}\right)_p \text{d}T + \left(\frac{\partial h}{\partial p}\right)_T \text{d}p$$
  * **定义式**：
    $$c_p = \left(\frac{\partial h}{\partial T}\right)_p \quad [\text{kJ/(kg}\cdot\text{K)}]$$
  * **定压过程焓变计算**：
    $$\Delta h = \int_{1}^{2} c_p \text{d}T$$

---

## 模块三：纯物质相变与状态区域划分

### 3.1 基本概念
* **纯物质（Pure Substance）**：化学组成在全系统内完全均匀且恒定的物质（如纯水、冰水共存体系、未发生化学反应的均相空气）。
* **相（Phase）**：物理结构与化学组成均一的聚集状态（固相、液相、气相）。

### 3.2 纯物质加热相变全过程（$T-v$ 图）
恒压下对工质持续加热，经历五个典型阶段：

```
[过冷液体]  -->  [饱和液体]  -->  [汽液两相共存区]  -->  [饱和蒸汽]  -->  [过热蒸汽]
 (T < Tsat)      (x = 0)        (T = Tsat, 沸腾)      (x = 1)       (T > Tsat)
```


1. **过冷水 / 压缩液体（Subcooled / Compressed Liquid）**：
   * 实际温度低于该压力下的饱和温度（沸点）：$T < T_{\text{sat}}(p)$。
2. **两相混合区（Liquid-Vapor Mixture）**：
   * 液体与气体处于动态平衡共存，相变过程等温等压。
   * **关键特性**：在两相区内，**压力 $p$ 与温度 $T$ 完全锁定对应，不再是两个相互独立的变量**。不能仅用 $(p, T)$ 确定状态。
3. **过热蒸汽（Superheated Vapor）**：
   * 液体完全汽化后继续吸热，实际温度高于该压力下的饱和温度：$T > T_{\text{sat}}(p)$。
4. **临界点（Critical Point）**：
   * 气相与液相能够平衡共存的最高温度 $T_c$ 与最高压力 $p_c$。高于临界温度时，物质不能通过等温压缩而液化。

### 3.3 干度（Quality, $x$）及杠杆法则
* **干度定义**：汽液两相区中，**饱和气体的质量占两相总质量的百分比**：
  $$x = \frac{m_{\text{vapor}}}{m_{\text{total}}} = \frac{m_g}{m_f + m_g} \quad (0 \le x \le 1)$$
  * $x = 0$：全部为饱和液体（点 $f$）
  * $x = 1$：全部为饱和干蒸汽（点 $g$）
* **两相区混合物性质通式（杠杆法则）**：
  若 $\Pi$ 代表比容 $v$、比内能 $u$、比焓 $h$ 或比熵 $s$：
  $$\Pi = (1 - x)\Pi_f + x\Pi_g = \Pi_f + x(\Pi_g - \Pi_f) = \Pi_f + x\Pi_{fg}$$
  * 展开式例：
    * $v = v_f + x v_{fg} = v_f + x(v_g - v_f)$
    * $u = u_f + x u_{fg} = u_f + x(u_g - u_f)$
    * $h = h_f + x h_{fg} = h_f + x(h_g - h_f)$

---

## 模块四：热力学性质查表方法与工程计算技能

### 4.1 水蒸气表结构及标识规范
* 表 A.1.1 / A.1.2：**饱和蒸汽表**（以温度 $T$ 或压力 $p$ 为入口索引）
  * 下标 $f$（Fluid）：饱和液体状态
  * 下标 $g$（Gas）：饱和蒸汽状态
  * 下标 $fg$：汽化变化量（$\Pi_{fg} = \Pi_g - \Pi_f$）
* 表 A.1.3：**过热蒸汽表**（以 $p$ 和 $T$ 共同作为入口索引）
* 表 A.1.4：**压缩液体表**

### 4.2 状态判定工作流
已知某状态的压力 $p$ 和温度 $T$ 时：
1. 查饱和表得到该压力下的饱和温度 $T_{\text{sat}}(p)$：
   * 若 $T < T_{\text{sat}}$：**压缩液体**
   * 若 $T = T_{\text{sat}}$：**饱和两相混合物**（需配合 $v$、$h$ 或 $x$ 进一步确定）
   * 若 $T > T_{\text{sat}}$：**过热蒸汽**（直接查 Table A.1.3）
2. 若已知 $p$ 与比容 $v$（或比焓 $h$）：
   * 查对应压力下的 $v_f$ 和 $v_g$：
     * $v < v_f$：压缩液体
     * $v_f \le v \le v_g$：两相混合物，反求干度 $x = \frac{v - v_f}{v_{fg}}$
     * $v > v_g$：过热蒸汽

### 4.3 线性插值法（Linear Interpolation）
当目标参数未在表格中直接列出时，利用相邻两点构成的斜率进行线性估算：
$$\frac{\Pi - \Pi_L}{\Pi_H - \Pi_L} = \frac{T - T_L}{T_H - T_L}$$
$$\Pi = \Pi_L + \left(\frac{T - T_L}{T_H - T_L}\right)(\Pi_H - \Pi_L)$$
*(注：下标 $L$ 为表格低位节点，$H$ 为高位节点)*

### 4.4 压缩液体的工程近似法
由于水等液态物质可压缩性极小，压力对其体积和分子热运动影响微弱：
* **比容与内能**：直接近似取该温度下的饱和水数值（仅取决于温度 $T$）：
  $$v(T, p) \approx v_f(T)$$
  $$u(T, p) \approx u_f(T)$$
* **比焓（包含压力做功修正）**：
  $$h(T, p) \approx u_f(T) + p v_f(T) \approx h_f(T) + v_f(T)\left[p - p_{\text{sat}}(T)\right]$$
* **常温水物性经验常数**：
  $$c_p \approx c_v \approx 4.2\text{ kJ/(kg}\cdot\text{K)}$$
  $$\Delta u \approx \Delta h \approx 4.2 \times \Delta T\quad (T\text{ in }^\circ\text{C})$$

---

## 模块五：典型例题与工程计算示范

### 例题 1：两相区干度与性质求解
**【题目】**
刚性密闭容器中盛有 $m = 5\text{ kg}$ 的水-汽混合物，温度为 $100^\circ\text{C}$。已知容器内气相体积为 $V_g = 3.3458\text{ m}^3$。
查表 Table A.1.1 可得 $100^\circ\text{C}$ 下的饱和参数：
* $p_{\text{sat}} = 101.325\text{ kPa}$
* $v_f = 0.001044\text{ m}^3/\text{kg}$
* $v_g = 1.6729\text{ m}^3/\text{kg}$
* $u_f = 418.94\text{ kJ/kg}, \quad u_{fg} = 2087.6\text{ kJ/kg}$

求：
1. 系统的干度 $x$；
2. 混合物的平均比容 $v$ 和总内能 $U$。

**【解答】**
1. **计算蒸汽质量 $m_g$ 与干度 $x$**：
   $$m_g = \frac{V_g}{v_g} = \frac{3.3458\text{ m}^3}{1.6729\text{ m}^3/\text{kg}} = 2.0\text{ kg}$$
   系统总质量为 $5\text{ kg}$，故干度为：
   $$x = \frac{m_g}{m} = \frac{2.0}{5.0} = 0.4$$

2. **计算平均比容 $v$**：
   $$v = v_f + x(v_g - v_f) = 0.001044 + 0.4 \times (1.6729 - 0.001044) \approx 0.6698\text{ m}^3/\text{kg}$$

3. **计算比内能 $u$ 与总内能 $U$**：
   $$u = u_f + x u_{fg} = 418.94 + 0.4 \times 2087.6 = 1253.98\text{ kJ/kg}$$
   $$U = m \cdot u = 5.0\text{ kg} \times 1253.98\text{ kJ/kg} = 6269.9\text{ kJ}$$

---

### 例题 2：线性插值应用
**【题目】**
已知过热蒸汽在压力 $p = 1.0\text{ MPa}$（$10\text{ bar}$）下：
* $T_L = 200^\circ\text{C}$ 时，$v_L = 0.2060\text{ m}^3/\text{kg}$
* $T_H = 240^\circ\text{C}$ 时，$v_H = 0.2275\text{ m}^3/\text{kg}$

求 $T = 215^\circ\text{C}$ 时的比容 $v$。

**【解答】**
使用线性插值公式：
$$\frac{v - v_L}{v_H - v_L} = \frac{T - T_L}{T_H - T_L}$$
代入数据：
$$\frac{v - 0.2060}{0.2275 - 0.2060} = \frac{215 - 200}{240 - 200} = \frac{15}{40} = 0.375$$
$$v = 0.2060 + 0.375 \times (0.2275 - 0.2060) = 0.2060 + 0.00806 = 0.21406\text{ m}^3/\text{kg}$$

---

### 例题 3：高压压缩液体近似计算
**【题目】**
水处于 $T = 100^\circ\text{C}$、$p = 5.0\text{ MPa}$（$5000\text{ kPa}$）的压缩液体状态。在无压缩液体表的情况下，利用饱和水参数估算其比焓 $h$，并与真实值进行比对。
已知 $100^\circ\text{C}$ 饱和参数：$p_{\text{sat}} = 101.325\text{ kPa}$，$h_f = 419.05\text{ kJ/kg}$，$v_f = 0.001044\text{ m}^3/\text{kg}$；查 Table A.1.4 真实值 $h_{\text{exact}} = 422.73\text{ kJ/kg}$。

**【解答】**
利用压缩液体比焓近似公式：
$$h(T, p) \approx h_f(T) + v_f(T)\left[p - p_{\text{sat}}(T)\right]$$
代入参数（注意压力单位必须与 $\text{kJ}$ 换算一致：$1\text{ kPa}\cdot\text{m}^3 = 1\text{ kJ}$）：
$$h \approx 419.05 + 0.001044 \times (5000 - 101.325) = 419.05 + 5.11 = 424.16\text{ kJ/kg}$$
* 与真实值 $422.73\text{ kJ/kg}$ 相比，相对误差仅约 $0.34\%$，精度完全符合工程计算要求。