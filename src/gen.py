import sys
import random

SETN = 500
SETM = 5000
SETT = 1000
SETP = 10
SETD = 1000

# N = random.randint(5, SETN)
# M = random.randint(N - 1, N)
# T = random.randint(500, SETT)
# P = random.randint(2, SETP)
# D = random.randint(2, SETD)

N = 4500
M = 5000
T = 10000
P = 80
D = 1000

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

print(str(N) + ' ' + str(M) + ' ' + str(T) + ' ' + str(P) + ' ' + str(D))

for edge in edge_list:
    print(str(edge[0]) + ' ' + str(edge[1]) + ' ' + str(edge[2]))

for task in task_list:
    print(str(task[0]) + ' ' + str(task[1]))
