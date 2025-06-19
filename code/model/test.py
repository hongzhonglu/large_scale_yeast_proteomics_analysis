import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

# 分数列表
scores = [76, 80, 78, 91, 96, 85, 92, 84, 71, 77, 77, 74, 79, 77, 69, 68, 85, 80, 71, 81, 84, 74]

# 计算均值和标准差
mu, std = np.mean(scores), np.std(scores)

# 创建一个从最小值到最大值的等差数列作为x轴
x = np.linspace(min(scores) - 5, max(scores) + 5, 100)

# 绘制直方图
plt.hist(scores, bins=10, density=True, alpha=0.6, color='g')

# 绘制正态分布曲线
p = norm.pdf(x, mu, std)
plt.plot(x, p, 'k', linewidth=2)

# 标题与标签
plt.title('Normal Distribution of Scores')
plt.xlabel('Scores')
plt.ylabel('Probability Density')

# 显示图形
plt.show()
