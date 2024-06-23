---
title: 'Norm'
date: 2024-6-23
permalink: /posts/2024/06/Norm/
---

# 范数 （norm）

## 1. 向量范数

**定义 1**		对于任意 $n$ 维向量 $\pmb{x}\in\mathbb{R}^n$，若存在一个唯一的实数 $N(x) = ||\pmb{x}||\in \mathbb{R}$ 与 $x$ 对应，且满足

(1) 正定性	$||\pmb{x}||\geq 0$，且 $\forall \pmb{x}\in \mathbb{R}^n$，$||\pmb{x}||=0\iff x = 0$.

(2) 齐次性	$||\alpha \pmb{x}|| = |\alpha|\cdot ||\pmb{x}||$，$\forall \pmb{x}\in \mathbb{R}^n$，$\alpha \in \mathbb{R}$.

(3) 三角不等式	$||\pmb{x}+\pmb{y}|| \leq ||\pmb{x}||+||\pmb{y}||$，$\forall \pmb{x},\pmb{y}\in\mathbb{R}^n$.

则称 $N(\pmb{x})$ 为向量 $\pmb{x}$ 的范数。

常用的向量范数：$L_p$-范数，其一般形式为
$$
||\pmb{x}||_p = (\sum_{i = 1}^{n}x_i^p)^{\frac{1}{p}}$,\quad p>0
$$
当 $p=1,2,\infty$ 时，

(1) **$L_1$-范数**
$$
||\pmb{x}||_1 = \sum_{i = i}^{n}|x_i|
$$
$L_1$-范数表示的是向量各分量的绝对值之和，又被称为**曼哈顿距离**。

(2) **$L_2$-范数**
$$
||\pmb{x}||_2 = \sqrt{\sum_{i = 1}^{n}x_i^2}
$$
$L_1$-范数表示的是向量各分量平方和开根号，与我们平时计算两物体之间的直线距离用的勾股定理类似，该距离又被称为**欧氏距离**。

(3) **$L_\infty$-范数**
$$
||\pmb{x}||_{\infty} = \max_{1\leq i\leq n}|x_i|
$$
$L_{\infty}$-范数表示元素绝对值的最大值。

***

**一个特殊的“范数”——$L_0$-范数**
$$
||\pmb{x}||_0 = \text{card}\Big(\{x_i|x_i\neq0,\ 1\leq i\leq n\}\Big)
$$
这里的 $\text{card}(\cdot)$ 符号集合中元素的数量。$L_0$-范数表示向量中非零元素的个数，这在实际中有着广泛的应用，特别是在稀疏恢复领域。需要注意的是，虽然 $L_0$-范数叫“范数”，但它并不属于范数，因为它不满足向量范数定义的第三条，三角不等式。

***

## 2. 矩阵范数

**定义 2**		对于任意 $n$ 阶矩阵 $A\in\mathbb{R}^{n\times n}$，若存在一个唯一的实数 $N(A) = ||A||\in\mathbb{R}^{n\times n}$ 与 $A$ 对应，且满足

(1) 正定性	$||A||\geq 0$，且 $\forall A\in\mathbb{R}^{n\times n}$，$||A||=0\iff A = 0$.

(2) 齐次性	$||\alpha A|| = |\alpha|\cdot ||A||$，$\forall A\in\mathbb{R}^{n\times n}$，$\alpha \in \mathbb{R}$.

(3) 三角不等式	$||A+B|| \leq ||A||+||B||$，$\forall A,B\in\mathbb{R}^{n\times n}$.

(4) 	$AB\leq ||A||\cdot ||B||$，$\forall A,B\in\mathbb{R}^{n\times n}$.

则称 $N(A)$ 为矩阵 $A$ 的范数。

矩阵范数分为**诱导范数**和**非诱导范数**。

**定义 3**		设 $\pmb{x}\in\mathbb{R}^n$，$A\in\mathbb{R}^{n\times n}$，$||\cdot||_v$ 为一种向量范数，则 $\frac{||A\pmb{x}||_v}{||\pmb{x}||_v}$ 对所有的 $\pmb{x}\neq0$有最大值，定义
$$
||A||_v = \max_{\pmb{x}\neq0}\Big\{\frac{||A\pmb{x}||_v}{||\pmb{x}||_v}\Big\} = \max_{||\pmb{x}||_v=1}\{||A\pmb{x}||_v \}
$$
可以验证 $||A||_v$ 满足**定义 2** 的四个条件，称其为从属于给定向量范数 $||\pmb{x}||_v$ 的矩阵范数，简称从属范数或算子范数，也叫诱导范数。

常用的诱导范数有以下几种：

(1) **列和范数**
$$
||A||_1= \max_{\pmb{x}\neq0}\Big\{\frac{||A\pmb{x}||_1}{||\pmb{x}||_1}\Big\} = \max_{1<j<n}\sum_{i = 1}^n|a_{ij}|
$$
(2) **行和范数**
$$
||A||_\infty= \max_{\pmb{x}\neq0}\Big\{\frac{||A\pmb{x}||_\infty}{||\pmb{x}||_\infty}\Big\} = \max_{1<i<n}\sum_{j = 1}^n|a_{ij}|
$$
(3) **谱范数**
$$
||A||_2= \max_{\pmb{x}\neq0}\Big\{\frac{||A\pmb{x}||_2}{||\pmb{x}||_2}\Big\} = \sqrt{\lambda_{\max}(A^TA)}
$$
其中，$\lambda_{\max}(A^TA)$ 为 $A^TA$ 行列式的最大值。

常见的非诱导范数有

**Frobenius 范数**，简称 **$F$-范数**
$$
||A||_F = \sqrt{\sum_{i = 1}^n\sum_{j = 1}^{n}a_{ij}^2}
$$
容易验证，$F$-范数满足
$$
||A||_F = \sqrt{\text{tr}(A^TA)} = \sqrt{\text{tr}(AA^T)}
$$

[^1]: 李庆扬, 王能超, 易大义, 数值分析（第五版）, 清华大学出版社, 2008
[^2]: 施吉林, 刘淑珍, 陈桂芝. 计算机数值方法 (第三版)[M]. 高等教育出版社, 2009.4.
