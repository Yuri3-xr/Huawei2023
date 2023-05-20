import sys
import random
import math

SETN = 5000
SETM = 5000
SETT = 10000
SETR = 35000
SETP = 80
SETD = 1000

N = 10
M = 20
T = 100
R = 0
P = 20
D = random.randint(2, SETD)

fa = [i for i in range(N)]
edge_list = []
task_list = []

def get_fa(x):
    if x == fa[x]:
        return x
    else:
        return get_fa(fa[x])

def gen_span_tree():
    cnt = N - 1
    while (cnt > 0):
        while (True):
            u = random.randint(0, N - 1)
            v = random.randint(0, N - 1)
            fu = get_fa(u)
            fv = get_fa(v)
            if (fu == fv):
                continue
            else:
                fa[fu] = fv
                edge_list.append((u, v, random.randint(1, D)))
                break
        cnt -= 1

def random_add_edge():
    cnt = M - N + 1
    while (cnt > 0):
        while (True):
            u = random.randint(0, N - 1)
            v = random.randint(0, N - 1)
            if (u == v):
                continue
            else:
                edge_list.append((u, v, random.randint(1, D)))
                break
        cnt -= 1

def random_task():
    cnt = T
    while (cnt > 0):
        while (True):
            u = random.randint(0, N - 1)
            v = random.randint(0, N - 1)
            if (u == v):
                continue
            else:
                task_list.append((u, v))
                break
        cnt -= 1

gen_span_tree()
random_add_edge()
random_task()

r_list = []

for i in range(T):
    r_list.append(random.randint(1, 10))
    R += r_list[i]
    i += 1

print(str(N) + ' ' + str(M) + ' ' + str(T) + ' ' + str(R) + ' ' + str(P) + ' ' + str(D))

for edge in edge_list:
    print(str(edge[0]) + ' ' + str(edge[1]) + ' ' + str(edge[2]))

cnt = 0
for task in task_list:
    print(str(task[0]) + ' ' + str(task[1]) + ' ' + str(r_list[cnt]))
    cnt += 1


