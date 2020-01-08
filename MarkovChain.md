Markov Chain
================
Lu Yiyun
2020/1/4

  - [1. 马尔科夫过程](#马尔科夫过程)
  - [2. 马尔科夫链](#马尔科夫链)
      - [例2](#例2)
      - [例3](#例3)
      - [齐次马氏链的有限维分布](#齐次马氏链的有限维分布)
  - [3. 多步转移概率](#多步转移概率)
  - [4. 遍历性](#遍历性)
      - [例1](#例1)
  - [5. 细致平稳条件](#细致平稳条件)

# 1\. 马尔科夫过程

**马尔可夫性**，又称为无后效性，即过程在时刻\(t_0\)所处的状态为已知的条件下，过程在时刻\(t>t_0\)所处状态的条件分布和在\(t_0\)之前的状态无关。设随机过程\({X(t),t\in T}\)的状态空间为\(I\)，如果对于时间\(t\)的任意\(n\)个数值\(t_1<t_2<\cdots<t_n,n\ge3,t_i\in T\)，其数学描述为：
\[
P\{X(t_n)\le x_n|X(t_1)=x_1,X(t_2)=x_2,\cdots,X(t_{n-1})=x_{n-1}\}=P\{X(t_n)\le x_n|X(t_{n-1})=x_{n-1}\},x_n\in R
\] 也称这个过程是一个**马尔科夫过程**。

> **定理1**：设\(X(t),t\ge0\)是独立增量过程，且\(X(0)=0\)，则其也是一个马尔科夫过程。

则我们可知，**泊松过程是时间连续、状态离散的马氏过程，维纳过程是时间状态都连续的马氏过程。**

# 2\. 马尔科夫链

时间状态都离散的马氏过程称为**马尔科夫链**。此时，其条件分布律就刻画了它的性质，即对于任意的正整数\(n,r\)和\(0\le t_1<t_2<\cdots<t_r<m;t_i,m,m+n\in T_1\)，有：
\[
P\{X_{m+n}=a_j|X_{t_1}=a_{i_1},X_{t_2}=a_{i_2},\cdots,X_{t_r}=a_{i_r},X_m=a_i\}=P\{X_{m+n}=a_j|X_m=a_i\}
\]

我们将以上的概率记为\(P_{ij}(m,m+n)\)，称为马氏链在时刻\(m\)处于状态\(a_i\)条件下，在时刻\(m+n\)转移到状态\(a_j\)的**转移概率**。

并且由概率的性质我们知道： \[
\sum_{j=1}^{+\infty}{P_{ij}(m, m+n)}=1,\ i=1,2,\cdots
\]

如果这个概率只和\(i,j,n\)相关，则可以记为\(P_{ij}(n)\)，此时称此转移概率具有**平稳性**，同时称此链是**齐次的**或**时齐的**。（为什么叫时齐的，因为对于每一个时间点转移概率都是一样的）

> **后面我们只会讨论齐次马氏链。**

\(P_{ij}(n)\)称为**n步转移概率**，其中**1步转移概率**最重要，如果状态有限，则我们可以将所有的1步转移概率构成一个矩阵的形式：

![](../images/transition_matrix.png)

## 例2

![](../images/markov_exam2.png)

## 例3

![](../images/markov_exam3-1.png) ![](../images/markov_exam3-2.png)
![](../images/markov_exam3-3.png)

我对这个模拟感兴趣，就计算机模拟了一下：

``` r
library(tidyverse)

random_walk_1d_1step <- function(humans, trasition_matrix, n) {
  # 多个人进行一步
  humans_transition_prop <- trasition_matrix[humans, ]
  apply(humans_transition_prop, 1, function(x) sample.int(n, size = 1, prob = x))
}

random_walk_1d <- function(humans, steps=1000, attractors=NULL) {
  # 创建转移矩阵
  n_hum <- length(humans)
  trasition_matrix <- matrix(c(0, 1, 0, 0, 0,
                               1/3, 1/3, 1/3, 0, 0,
                               0, 1/3, 1/3, 1/3, 0,
                               0, 0, 1/3, 1/3, 1/3,
                               0, 0, 0, 1, 0),
                             nrow = 5,
                             byrow = TRUE)
  if (!is.null(attractors)) {
    trasition_matrix[attractors,] <- diag(rep(1, 5))[attractors]
  }
  # 进行游走
  res_mat <- matrix(0, nrow = steps+1, ncol = n_hum,
                    dimnames = list(NULL, str_c("human", seq_len(n_hum))))
  for (i in seq_len(steps)) {
    res_mat[i,] <- humans
    humans <- random_walk_1d_1step(humans, trasition_matrix, 5)
  }
  res_mat[i+1,] <- humans
  res_mat
}
```

``` r
humans <- c(1, 2, 3)
rw_res <- random_walk_1d(humans, 50)
as_tibble(rw_res) %>% mutate(step=1:nrow(rw_res)) %>%
  gather(key = "humans", value = "position", -step) %>% 
  ggplot() + geom_line(aes(step, position, colour = humans)) +
  theme_bw()
```

![](/images/MarkovChain/unnamed-chunk-2-1.png)<!-- -->

## 齐次马氏链的有限维分布

![](../images/markov_dist.png) ![](../images/markov_dist2.png)

# 3\. 多步转移概率

![](../images/ck_equation.png)

矩阵形式： \[P(u+v)=P(u)P(v)\]

进一步分解，得到： \[P(n)=P^n\]

**对于齐次马氏链，n步转移概率矩阵是一步转移概率矩阵的n次方，进而齐次马氏链的有限维分布由初始分布和一步转移概率完全确定。**

![](../images/2step_markov.png)

# 4\. 遍历性

![](../images/2step_markov2.png)

这说明有些马尔科夫链在经过长时间的转移后，其单时间点的分布会趋近于一个确定的分布，这种性质称为**遍历性**。

遍历性的定义： ![](../images/bianli.png)

那什么时候齐次马氏链有遍历性呢？有以下定理：

![](../images/bianli_proof.png)

注意，式(3.2)并不能保证齐次马氏链是遍历的，但可以用来求解极限分布。

> 对于一个齐次马氏链的转移矩阵，其中的元素值肯定非负。所以如果只要能够通过矩阵的幂运算将矩阵中的0填满，就可以认为这个马氏链是具有遍历性的。另外，极限分布只和转移矩阵有关，而和初始分布无关。

以下的例子介绍了一个在简单情况下判别是否具有遍历性的方法，并求解其极限分布：

## 例1

![](../images/bianli_exam1.png)

![](../images/bianli_exam1-2.png)

# 5\. 细致平稳条件

虽然上面介绍的是离散马尔科夫链，我们可以将状态空间也推广至连续空间。时间连续、状态空间离散的马尔科夫链拥有相似的性质，即也有极限分布。这时候，我们回到采样的问题里：**找到一个极限分布是目标分布的马尔科夫链，然后我们对此随机过程采样，尽量取其晚时间点的样本，此时这些样本的分布接近目标分布。**

这其实就是马尔科夫链蒙特卡洛方法的基本思想。为了能够完成这样的任务，我们需要能够对指定的分布来构造一个马尔科夫链（对于齐次马氏链，就是确定转移概率），幸运地是，我们拥有这样的一个充分条件：

> 细致平稳条件：如果具有遍历性的马尔科夫链的转移概率\(P\{X_{i+1}=x|X_i=y\}\)和一个分布\(\pi(x)\)满足以下条件（称为细致平稳条件），则该分布\(\pi(x)\)就是该马尔科夫链的极限分布：
> \[
> \pi(y)P\{X_{i+1}=x|X_i=y\}=\pi(x)P\{X_{i+1}=y|X_i=x\}
> \]

可以简单地证明其符合极限分布的定义。

所以我们只需要构建符合上述细致平稳分布的转移概率就可以了。

关于马尔科夫链，还可以进行探索的方向：

  - 马尔科夫链的收敛速度；
  - 矩阵论；
  - 连续马尔科夫过程等。
