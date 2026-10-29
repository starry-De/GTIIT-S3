# 工程热力学：纯物质性质与闭口系统能量分析复习笔记

---

## 一、 核心状态参数（Properties）

> **符号规范**：大写字母代表总体的广延量（Extensive），小写字母代表单位质量（$1\text{ kg}$）的强度量/比参数（Specific property）。

1. **压力 $P$（Pressure）**
   * **常用单位换算**：
     $$1\text{ bar} = 100\text{ kPa} = 0.1\text{ MPa}$$
     $$1\text{ MPa} = 1000\text{ kPa}$$
   * **计算准则**：计算功时统一将压力化为 $\text{kPa}$，体积化为 $\text{m}^3$，乘积直接为千焦（$1\text{ kPa}\cdot\text{m}^3 = 1\text{ kJ}$）。
2. **温度 $T$（Temperature）**
   * 常用单位：$^\circ\text{C}$ 或 $\text{K}$（$T(\text{K}) = T(^\circ\text{C}) + 273.15$）。
3. **比容 $v$（Specific Volume）**
   * 定义：单位质量物质所占的体积（密度 $\rho$ 的倒数）。
     $$v = \frac{V}{m} \quad (\text{m}^3/\text{kg})$$
4. **比内能 $u$（Specific Internal Energy）**
   * 定义：单位质量物质内部微观分子动能与势能之和（$\text{kJ/kg}$）。
   * 总内能：$$U = m u \quad (\text{kJ})$$
5. **比焓 $h$（Specific Enthalpy）**
   * 定义：为简化定压与流动过程计算而定义的组合状态参数。
     $$h = u + P v \quad (\text{kJ/kg})$$
   * 总焓：$$H = m h = U + P V \quad (\text{kJ})$$

---

## 二、 纯物质相变过程与 $P-v$ 相图（Vapor Dome）

真实工质（如水/水蒸气）存在相变，**严禁使用理想气体状态方程 $PV = mRT$**，一切参数需查水蒸气表（Steam Tables）。

```text
  P (压力)
  ^
  |               [临界点 Critical Point]
  |                      /\
  |                     /  \
  |   (压缩液体区)      /    \       (过热蒸汽区)
  |                   /      \
P_sat ---------------B--------C---------> D  (定压加热线)
  |                 /|   湿   |\
  |      A---------/ | 蒸汽区 | \
  |               /  | (水+汽)|  \
  +--------------+---+--------+---+---------> v (比容 m³/kg)
                 vf           vg
```

### 1. 五个关键状态与三个相区
* **压缩液体 / 过冷水（Compressed / Subcooled Liquid，点 A）**：温度低于当前压力下的沸点（$T < T_{\text{sat}}$），或压力高于当前温度下的饱和压力（$P > P_{\text{sat}}$）。水分完全为液态，极难压缩。
* **饱和液体（Saturated Liquid，点 B，下标 $f$）**：刚达到沸点、即将产生第一颗气泡的极限液态（干度 $x = 0$）。参数记为 $v_f, u_f, h_f, s_f$。
* **湿蒸汽区 / 汽液两相区（Wet Vapor / Two-phase Mixture，线段 B–C）**：沸腾相变进行中，水与蒸汽共存。**在定压条件下，相变过程温度保持不变（$T = T_{\text{sat}}$）**。
* **干饱和蒸汽（Saturated Vapor，点 C，下标 $g$）**：最后一滴水刚好完全汽化的极限气态（干度 $x = 1$）。参数记为 $v_g, u_g, h_g, s_g$。
* **过热蒸汽（Superheated Vapor，点 D）**：温度高于沸点（$T > T_{\text{sat}}$）的纯蒸汽。

### 2. 饱和钟形线（Vapor Dome）的构成
* **左边界（饱和液体线）**：不同压力下 $v_f$ 点的连线。
* **右边界（饱和蒸汽线）**：不同压力下 $v_g$ 点的连线。
* **临界点（Critical Point）**：左右边界在顶部交汇点（水：$P_c = 22.06\text{ MPa}, T_c = 373.95^\circ\text{C}$），此处 $v_f = v_g$。

---

## 三、 干度计算与状态定位法则

### 1. 干度 $x$（Quality）与杠杆公式
仅适用于**湿蒸汽两相区（拱门内部，$0 \le x \le 1$）**：
$$x = \frac{m_{\text{vapor}}}{m_{\text{total}}} = \frac{m_g}{m_f + m_g}$$

**两相区物性计算公式（加权平均）：**
$$v = v_f + x (v_g - v_f) = v_f + x v_{fg} \implies x = \frac{v - v_f}{v_g - v_f}$$
$$u = u_f + x (u_g - u_f) = u_f + x u_{fg}$$
$$h = h_f + x (h_g - h_f) = h_f + x h_{fg}$$
*(注：下标 $fg$ 代表汽化潜热/相变增量，即 $\phi_{fg} = \phi_g - \phi_f$)*。

### 2. 状态定位速查表（已知任意 2 个独立参数确定相态）

| 已知 $(P, v)$ 判据 | 已知 $(P, T)$ 判据 | 所在相区 | 查表与取值方法 |
| :--- | :--- | :--- | :--- |
| $v < v_f$ | $T < T_{\text{sat}}$ 或 $P > P_{\text{sat}}$ | **压缩液体区** | 液体性质几乎不受压力影响，近似取**同温度 $T$** 下的饱和水参数：<br>$v \approx v_f(T),\quad u \approx u_f(T)$ |
| $v = v_f$ | 题目指明“饱和液体” | **饱和液体线** | 直接查饱和表取 $v_f, u_f, h_f$（$x = 0$） |
| $v_f < v < v_g$ | $T = T_{\text{sat}}$ | **湿蒸汽两相区** | 先由 $x = \frac{v - v_f}{v_g - v_f}$ 求干度，再代入杠杆公式算 $u, h$ |
| $v = v_g$ | 题目指明“饱和蒸汽” | **饱和蒸汽线** | 直接查饱和表取 $v_g, u_g, h_g$（$x = 1$） |
| $v > v_g$ | $T > T_{\text{sat}}$ 或 $P < P_{\text{sat}}$ | **过热蒸汽区** | 查**过热蒸汽表（Superheated Tables）**双重定位 $(P, T)$ |

---

## 四、 功、热量与热力学第一定律（闭口系统）

### 1. 常见的机械约束与过程转换
* **自由浮动活塞（无卡阻）**：受力平衡 $P_{\text{内}} = P_{\text{atm}} + \frac{m_{\text{piston}}g}{A} = \text{常数} \implies$ **定压过程（Isobaric，$P = \text{const}$）**。
* **活塞触碰限位块（Stops）或刚性容器（Rigid Tank）**：边界固定无法移动 $\implies$ **定容过程（Isochoric，$V = \text{const}$）**。
* **绝热容器（Insulated）**：系统与外界无温差传热 $\implies$ **绝热过程（Adiabatic，$Q = 0$）**。

### 2. 功的计算（Work Interactions）
1. **边界容积变化功（Boundary Work, $W_b$）**：
   $$W_b = \int_{1}^{2} P \, dV$$
   * **定容过程（$V = \text{const}$）**：
     $$W_b = 0$$
   * **定压过程（$P = \text{const}$）**：
     $$W_b = P (V_2 - V_1) = m P (v_2 - v_1)$$
   * **几何意义**：$P-V$ 图上过程曲线下方的投影面积（膨胀对外做功为正，压缩对内做功为负）。
2. **轴功 / 搅拌功 / 电功（$W_{\text{shaft}}, W_{\text{mixer}}, W_{\text{elec}}$）**：
   外部通过旋转轴或导线向系统单向输入的功，不可用 $\int P dV$ 计算，需由热力学第一定律求解。

### 3. 闭口系统热力学第一定律（能量守恒方程）
**通用标准形式（规定吸热 $Q>0$，系统对外净输出功 $W_{\text{net,out}}>0$）：**
$$Q - W_{\text{net,out}} = \Delta U = m (u_2 - u_1)$$

**展开净功项（含边界输出功 $W_b$ 与搅拌输入功 $W_{\text{mixer,in}}$）：**
$$Q - (W_b - W_{\text{mixer,in}}) = m (u_2 - u_1)$$
* **定压过程特殊简化**：因 $W_b = P(V_2 - V_1)$，移项可得：
  $$Q + W_{\text{mixer,in}} = \Delta U + P\Delta V = \Delta H = m (h_2 - h_1)$$

---

## 五、 典型题型标准解题模板

### 题型 1：带限位块（Stops）的活塞-汽缸两阶段加热问题
* **已知条件**：$m = 2\text{ kg}$ 水，初态 $T_1 = 20^\circ\text{C}, P_1 = 5\text{ bar} = 500\text{ kPa}$。加热使活塞上升，在 $V_2 = 0.6\text{ m}^3$ 处碰到限位块；继续加热至成为干饱和蒸汽（状态 3）[cite: 1]。求：(a) 终态压力 $P_3$；(b) 全过程功 $W$ 与热量 $Q$[cite: 1]。
* **标准解答步骤**：
  1. **确定初态 1**：
     查 $20^\circ\text{C}$ 饱和表，$P_{\text{sat}} = 2.339\text{ kPa} < 500\text{ kPa} \implies$ **压缩液体**。
     $$v_1 \approx v_f(20^\circ\text{C}) = 0.001002\text{ m}^3/\text{kg} \implies V_1 = m v_1 = 0.002\text{ m}^3$$
     $$u_1 \approx u_f(20^\circ\text{C}) = 83.9\text{ kJ/kg}$$
  2. **拆解过程 $1 \to 2 \to 3$**：
     * $1 \to 2$（碰限位块前）：**定压膨胀**，$P_2 = P_1 = 500\text{ kPa}$，$V_2 = 0.6\text{ m}^3$[cite: 1]。
     * $2 \to 3$（碰限位块后）：**定容升温升压**，$V_3 = V_2 = 0.6\text{ m}^3$[cite: 1]。
  3. **确定终态 3 并求 $P_3$**：
     $$v_3 = \frac{V_3}{m} = \frac{0.6}{2} = 0.3\text{ m}^3/\text{kg}$$
     终态为饱和蒸汽 $\implies v_g(P_3) = v_3 = 0.3\text{ m}^3/\text{kg}$[cite: 1]。
     查压力饱和表（在 $600\text{ kPa}$ 与 $700\text{ kPa}$ 之间线性插值）：
     $$P_3 \approx 636.4\text{ kPa} = \mathbf{6.36\text{ bar}}, \quad u_3 = u_g(P_3) \approx \mathbf{2568.6\text{ kJ/kg}}$$
  4. **计算总功 $W$ 与总热量 $Q$**：
     $$W = W_{1\to2} + W_{2\to3} = P_1(V_2 - V_1) + 0 = 500 \times (0.6 - 0.002) = \mathbf{299\text{ kJ}}$$
     $$Q = \Delta U + W = m(u_3 - u_1) + W = 2 \times (2568.6 - 83.9) + 299 = \mathbf{5268.4\text{ kJ}}$$

---

### 题型 2：绝热定压搅拌（Mixer）相变问题
* **已知条件**：绝热汽缸含 $m = 1\text{ kg}$ 水蒸气，初态 $P_1 = 0.4\text{ MPa} = 400\text{ kPa}, V_1 = 20\text{ L} = 0.02\text{ m}^3$。搅拌器输入轴功，保持压力恒定，直至达到干饱和蒸汽（状态 2）[cite: 2]。求：(a) $P-v$ 图；(b) $\Delta u$ 与 $\Delta h$；(c) 所有功的交互量[cite: 2]。
* **标准解答步骤**：
  1. **确定初态 1**：
     $$v_1 = \frac{V_1}{m} = 0.02\text{ m}^3/\text{kg}$$
     查 $0.4\text{ MPa}$ 饱和表：$v_f = 0.001084,\ v_g = 0.46242,\ u_f = 604.22,\ u_g = 2553.1,\ h_f = 604.66,\ h_g = 2738.1$。
     因 $v_f < v_1 < v_g \implies$ **湿蒸汽两相区**。
     $$x_1 = \frac{v_1 - v_f}{v_g - v_f} = \frac{0.02 - 0.001084}{0.46242 - 0.001084} \approx 0.0410$$
     $$u_1 = u_f + x_1(u_g - u_f) = 684.12\text{ kJ/kg}$$
     $$h_1 = h_f + x_1(h_g - h_f) = 692.13\text{ kJ/kg}$$
  2. **确定终态 2 并求 $\Delta u, \Delta h$**：
     定压至饱和蒸汽（$P_2 = 0.4\text{ MPa}, x_2 = 1$）[cite: 2]：
     $$v_2 = v_g = 0.46242\text{ m}^3/\text{kg},\quad u_2 = u_g = 2553.1\text{ kJ/kg},\quad h_2 = h_g = 2738.1\text{ kJ/kg}$$
     $$\Delta u = u_2 - u_1 = 2553.1 - 684.12 = \mathbf{1869.0\text{ kJ/kg}}$$
     $$\Delta h = h_2 - h_1 = 2738.1 - 692.13 = \mathbf{2046.0\text{ kJ/kg}}$$
  3. **计算各项功交互（Work Interactions）**：
     * **边界膨胀功（系统对外输出）**：
       $$W_b = m P (v_2 - v_1) = 1 \times 400 \times (0.46242 - 0.02) = \mathbf{+177.0\text{ kJ}}$$
     * **搅拌轴功（外界对系统输入）**：
       由绝热条件 $Q = 0$ 及第一定律 $0 - (W_b - W_{\text{mixer,in}}) = m\Delta u$：
       $$W_{\text{mixer,in}} = m\Delta u + W_b = m\Delta h = \mathbf{2046.0\text{ kJ}} \quad (\text{或记为 } W_{\text{mixer}} = -2046.0\text{ kJ})$$
     * **净功**：
       $$W_{\text{net}} = W_b - W_{\text{mixer,in}} = 177.0 - 2046.0 = \mathbf{-1869.0\text{ kJ}}$$