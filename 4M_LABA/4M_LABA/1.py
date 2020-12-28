import pandas as pd
import numpy as np

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt


df = pd.read_csv('excelent.txt')
str_n = df.columns[3]
str_m = df.columns[4]
df = df.drop([str_n, str_m], axis=1)
df.columns=['x', 't', 'v']
df.info()
print(df)

n = int(str_n)
m = int(str_m)
print(n,m)

X = []
T = []
V = []
for i in range(n + 1):
    x = []
    t = []
    v = []
    for j in range(m + 1):
        x.append(df.loc[i * (m + 1) + j, 'x'])
        t.append(df.loc[i * (m + 1) + j, 't'])
        v.append(df.loc[i * (m + 1) + j, 'v'])
    X.append(x)
    T.append(t)
    V.append(v)
X = np.array(X)
T = np.array(T)
V = np.array(V)


fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')
# Grab some test data.
#X, Y, Z = axes3d.get_test_data(0.05)
#print(X)
# Plot a basic wireframe.
ax.plot_surface(T, X, V)
ax.set_xlabel('T(время)')
ax.set_ylabel('X(координата сечения)')
ax.set_zlabel('V(температура)')
ax.set_title('График полученной функции')
plt.show()