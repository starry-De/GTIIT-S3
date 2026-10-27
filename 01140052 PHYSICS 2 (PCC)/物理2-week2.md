# Electrostatics 电场计算核心复习笔记

适用范围：真空中的静电场；均匀介质中可在相应模型成立时将 $\varepsilon_0$ 换为介电常数 $\varepsilon$。

本笔记依据原对话中可读取的三道题目图片重新整理，并补充必要的通用公式。原笔记正文未能从对话引用中恢复，因此此文件是重新编写的复习版。

## 1. 核心公式与方向判断

令 $\mathbf r$ 为场点位置，$\mathbf r'$ 为源点位置，从源点指向场点的位移为

$$
\mathbf R=\mathbf r-\mathbf r',\qquad R=\lVert\mathbf R\rVert.
$$

库仑定律及叠加原理：

$$
d\mathbf E=k\frac{dq}{R^3}\mathbf R,
\qquad
\mathbf E=\int d\mathbf E,
\qquad
k=\frac{1}{4\pi\varepsilon_0}\approx8.99\times10^9\ \mathrm{N\,m^2/C^2}.
$$

- 正电荷的电场从源电荷向外，负电荷的电场指向源电荷。
- $dq$ 可以带符号；用矢量公式时，不要重复添加电荷正负号。
- 先按分量相加，再求大小：$E=\sqrt{E_x^2+E_y^2+E_z^2}$。
- 电场分量可以为负；电场大小（magnitude）必须非负。
- 电场方向由源电荷决定；受力还取决于试探电荷的正负。

$$
\mathbf F=q\mathbf E.
$$

正试探电荷受力与电场同向，负试探电荷反向。计算某电荷的受力时，使用其他电荷产生的外电场，不计它自身的电场。

## 2. 连续电荷积分：固定流程

1. 画出场点与源点，明确坐标轴正方向。
2. 利用对称性判断哪些分量抵消。
3. 选积分变量，写出 $dq$。
4. 写清距离 $R$ 和投影因子，再写积分限。
5. 积分后分别报告方向、分量和大小，并检查单位。

### 2.1 电荷微元

| 分布 | 密度 | 电荷微元 | 密度单位 |
| --- | --- | --- | --- |
| 线电荷 | $\lambda$ | $dq=\lambda\,dl$ | $\mathrm{C/m}$ |
| 面电荷 | $\sigma$ | $dq=\sigma\,dA$ | $\mathrm{C/m^2}$ |
| 体电荷 | $\rho$ | $dq=\rho\,dV$ | $\mathrm{C/m^3}$ |

非均匀密度必须保留在积分内。例如沿 $x'$ 轴的细线：

$$
dq=\lambda(x')\,dx',\qquad Q=\int_a^b\lambda(x')\,dx'.
$$

总电荷 $Q$ 通常不能直接代入 $kQ/R^2$：只有点电荷、满足条件的球对称分布外部等情形才能精确这样处理；远场则可能作为近似。

### 2.2 沿直线分布的通用分量式

源点为 $(x',0)$，场点为 $(x_0,y_0)$：

$$
E_x=k\int_a^b\frac{\lambda(x')(x_0-x')}{[(x_0-x')^2+y_0^2]^{3/2}}\,dx',
$$

$$
E_y=k\int_a^b\frac{\lambda(x')y_0}{[(x_0-x')^2+y_0^2]^{3/2}}\,dx'.
$$

**记忆：分母三次方配位移分量；分母二次方配单位方向或投影因子。**

### 2.3 对称性与常用积分

关于某轴成对对称的同号等量电荷，其垂直于对称轴的分量抵消，沿轴分量相加。几何形状与电荷密度都要满足相应对称性。

常见原函数（$a>0$）：

$$
\int\frac{dx}{(x^2+a^2)^{3/2}}
=\frac{x}{a^2\sqrt{x^2+a^2}}+C,
$$

$$
\int\frac{x\,dx}{(x^2+a^2)^{3/2}}
=-\frac{1}{\sqrt{x^2+a^2}}+C,
$$

$$
\int x^n\,dx=\frac{x^{n+1}}{n+1}+C\quad(n\ne-1),
\qquad
\int\frac{dx}{x}=\ln|x|+C.
$$

若积分端点距离场点为零，应按反常积分判断收敛，不能直接假定电场发散或有限。

## 3. 坐标、投影与几何微元

### 3.1 极坐标与柱坐标

平面极坐标：

$$
x=r\cos\phi,\qquad y=r\sin\phi,\qquad dA=r\,dr\,d\phi.
$$

半径固定为 $a$ 的圆弧：$dl=a\,d\phi$。圆盘上的完整薄环：$dA=2\pi r\,dr$。

柱坐标的体积元：

$$
dV=r\,dr\,d\phi\,dz.
$$

圆柱侧面半径固定为 $a$ 时：$dA=a\,d\phi\,dz$。不要把平面圆盘、圆柱侧面和圆锥侧面的面积元混用。

### 3.2 圆锥侧面：优先用母线距离

设圆锥半顶角为 $\theta$，从顶点沿母线量出的距离为 $s$，轴向正方向取为“顶点指向底面”。则

$$
r=s\sin\theta,\qquad z=s\cos\theta.
$$

母线宽度为 $ds$ 的薄环面积是“周长乘斜宽”：

$$
dA=2\pi r\,ds=2\pi s\sin\theta\,ds.
$$

若改用轴向坐标 $z$：

$$
ds=\frac{dz}{\cos\theta},\qquad
dA=\frac{2\pi z\sin\theta}{\cos^2\theta}\,dz.
$$

在顶点处计算电场时，源点至场点距离为 $s$，轴向投影因子为 $\cos\theta$；不能将距离误写为 $z$。

## 4. 可直接调用的二级结论

以下结论均需先核对几何条件与密度条件。

### 4.1 非均匀细线端点的电场

细线位于 $0<x\le L$，$\lambda(x)=Ax^n$，其中 $A>0$，场点在原点。每个电荷元产生的电场均沿 $-x$ 方向：

$$
E_x=-kA\int_0^L x^{n-2}\,dx.
$$

当 $n>1$ 时积分收敛：

$$
E_x=-\frac{kAL^{n-1}}{n-1},\qquad
E=\frac{kAL^{n-1}}{n-1}.
$$

当 $n\le1$ 且 $A\ne0$ 时，该理想细线模型在原点的场积分发散。若密度写成 $Ax^n$，系数 $A$ 的单位应为 $\mathrm{C/m^{n+1}}$。

### 4.2 均匀圆环轴线电场

半径为 $a$、总电荷为 $Q$ 的均匀圆环，圆心在原点、轴为 $z$ 轴：

$$
E_z(z)=\frac{kQz}{(z^2+a^2)^{3/2}}.
$$

横向分量因对称性抵消；$z=0$ 时电场为零。远处 $|z|\gg a$ 时，趋近同总电荷的点电荷电场。

### 4.3 均匀圆盘轴线电场

半径为 $a$、面密度为 $\sigma$ 的圆盘，在盘面上方 $z>0$：

$$
E_z(z)=\frac{\sigma}{2\varepsilon_0}
\left(1-\frac{z}{\sqrt{z^2+a^2}}\right).
$$

可由薄环积分得到。对于 $\sigma>0$，电场背离盘面；无限大平面在任意一侧的场强大小为 $|\sigma|/(2\varepsilon_0)$。上述圆盘公式适用于 $z>0$，不能直接代入负 $z$ 使用。

### 4.4 均匀圆锥截段在原顶点的电场

仅圆锥侧面带均匀面电荷，保留母线范围 $s_1\le s\le s_2$，其中 $0<s_1<s_2$。正轴从顶点指向底面：

$$
dE_z=-k\frac{\sigma\,dA}{s^2}\cos\theta
=-2\pi k\sigma\sin\theta\cos\theta\frac{ds}{s}.
$$

所以

$$
E_z=-2\pi k\sigma\sin\theta\cos\theta
\ln\frac{s_2}{s_1}.
$$

若 $\sigma>0$，方向为沿轴背离底面。去掉上半段后，$s_1=L/2$、$s_2=L$，故

$$
E=2\pi k|\sigma|\sin\theta\cos\theta\ln2.
$$

在此几何比例与密度固定时，结果与 $L$ 无关。若均匀带电侧面延伸至顶点，则积分出现对数发散；切去顶端是保证该模型结果有限的关键条件。

## 5. 匀强电场中的偏转

设粒子以水平速度 $v_0$ 入射，初始竖直速度为零，电场仅有竖直分量 $E_y$，板长为 $\ell$。忽略重力、空气阻力及边缘场，采用非相对论模型：

$$
t=\frac{\ell}{v_0},\qquad
a_y=\frac{qE_y}{m},\qquad
y=\frac{qE_y\ell^2}{2mv_0^2}.
$$

由出口偏移量反求电荷大小：

$$
|q|=\frac{2m|y|v_0^2}{|E_y|\ell^2}.
$$

若出口后再无场飞行水平距离 $D$，到屏幕的总偏移为

$$
y_{\mathrm{screen}}=
\frac{qE_y}{mv_0^2}\left(\frac{\ell^2}{2}+\ell D\right).
$$

必须区分题目要求的是“板末端偏移”还是“屏幕偏移”。失去电子后物体带正电，失去电子数为 $N=|q|/e$，其中 $e=1.602176634\times10^{-19}\ \mathrm C$。

## 6. 原题数值核对与易错点

### 6.1 非均匀线电荷

原题 $L=0.6\ \mathrm m$，将密度理解为 $\lambda=Ax^5$，其中 $A=4\times10^{-9}\ \mathrm{C/m^6}$：

$$
E_x=-\frac{kAL^4}{4}\approx-1.165\ \mathrm{N/C},
\qquad E\approx1.165\ \mathrm{N/C}.
$$

**题目问大小，应填写正值；负号只用于表示 $x$ 分量方向。**

### 6.2 去掉上半段的圆锥侧面

使用 $\sigma=6.2\ \mathrm{\mu C/m^2}$、$\theta=30^\circ$、$q=1\ \mathrm{\mu C}$：

$$
E\approx1.051\times10^5\ \mathrm{N/C},\qquad
F=|q|E\approx0.1051\ \mathrm N.
$$

正电荷受力沿轴背离底面。这里仅计算剩余圆锥侧面的作用，不添加底面电荷。

### 6.3 墨滴偏转

先换算质量：$1.5\times10^{-8}\ \mathrm g=1.5\times10^{-11}\ \mathrm{kg}$。

使用 $v_0=70\ \mathrm{m/s}$、$\ell=0.08\ \mathrm m$、$E=5\times10^4\ \mathrm{N/C}$、$|y|=7\times10^{-4}\ \mathrm m$：

$$
|q|=\frac{2(1.5\times10^{-11})(7\times10^{-4})(70)^2}
{(5\times10^4)(0.08)^2}
\approx3.216\times10^{-13}\ \mathrm C
=0.3216\ \mathrm{pC}.
$$

以上数值保留额外位数用于核算；提交答案时按题目要求取有效数字。

### 6.4 最后检查清单

| 检查项 | 正确做法 |
| --- | --- |
| 大小与分量 | 大小非负；方向单独说明 |
| 源点至场点 | 使用 $\mathbf r-\mathbf r'$ |
| 对称性 | 同时检查几何与电荷密度 |
| 距离与投影 | 分清 $R$、轴向距离和 $\cos\theta$ |
| 圆锥面积元 | 母线变量用 $2\pi s\sin\theta\,ds$ |
| 积分端点 | 场点贴近电荷分布时检查收敛性 |
| 单位换算 | $\mathrm{g}\to\mathrm{kg}$、$\mathrm{cm}\to\mathrm m$、$\mathrm{mm}\to\mathrm m$ |
| 电荷单位 | $1\ \mathrm{\mu C}=10^{-6}\ \mathrm C$；$1\ \mathrm{pC}=10^{-12}\ \mathrm C$ |
| 量纲 | 电场为 $\mathrm{N/C}$，力为 $\mathrm N$ |
| 极限 | 圆环远场趋近点电荷；对称中心检查抵消 |

**解题主线：判方向 → 用对称性 → 选微元 → 写距离与投影 → 积分 → 查单位和极限。**
