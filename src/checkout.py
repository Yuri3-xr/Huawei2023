import sys
import os

if len(sys.argv) < 3:
    exit()

exe = str(sys.argv[1])
input_file = str(sys.argv[2])
output_file = "./output.txt"

os.system('./' + exe + ' < ' + input_file + ' > ' + output_file)

edge_list = []
mat = []
edge_index = {}
channel_list = []
task_list = []
INF = 1e13
N = 0
M = 0 
T = 0
P = 0
D = 0
Y = 0
ans = 0

def add_edge(u, v, id, w):
    edge_list[u][v] = {'id': id, 'w': w}
    edge_index[id] = (u, v)
    mat[u][v] = min(mat[u][v], w)

def check(id, channel, path_list, station_list):
    vis = [0 for i in range(N)]
    for station in station_list:
        vis[station] = 1
    task = task_list[id]
    lst = task[0]
    to = task[1]
    now_dis = 0
    
    for edge_id in path_list:
        u = edge_index[edge_id][0]
        v = edge_index[edge_id][1]
        if u != lst:
            t = u
            u = v
            v = t
        if channel_list[edge_id][channel] != 0:
            return -1
        channel_list[edge_id][channel] = 1
        if u != lst:
            return -1
        now_dis += edge_list[u][v]['w']
        if now_dis > D:
            return -1
        if vis[v] == True:
            now_dis = 0
            pass
        lst = v

    if lst != to:
        return -1
    
    return 0


input_object = open(input_file)
try:
    line_list = input_object.readlines()
    first_line = line_list[0]
    par_list = first_line.split(' ')
    if len(par_list) < 5:
        print('error')
        exit()

    N = int(par_list[0])
    M = int(par_list[1])
    T = int(par_list[2])
    P = int(par_list[3])
    D = int(par_list[4])

    edge_list = [{} for i in range(N)]

    for i in range(M):
        channel_list.append([])
        channel_list[i] = [0 for i in range(P)]

    for i in range(N):
        mat.append([])
        mat[i] = [INF for i in range(N)]

    for i in range(1, M + 1):
        edge_info_list = line_list[i].split(' ')
        u = int(edge_info_list[0])
        v = int(edge_info_list[1])
        w = int(edge_info_list[2])
        add_edge(u, v, i - 1, w)
        add_edge(v, u, i - 1, w)
    
    for i in range(M + 1, 1 + M + T):
        task_info_list = line_list[i].split(' ')
        u = int(task_info_list[0])
        v = int(task_info_list[1])
        task_list.append((u, v))

finally:
    input_object.close()

output_object = open(output_file)
try:
    line_list = output_object.readlines()
    first_line = line_list[0]
    par_list = first_line.split(' ')
    Y = int(par_list[0])
    ans += Y * 1e6

    for i in range(Y):
        channel_list.append([])
        channel_list[i + M] = [0 for i in range(P)]
    
    for i in range(1, Y + 1):
        edge_info_list = line_list[i].split(' ')
        u = int(edge_info_list[0])
        v = int(edge_info_list[1])
        add_edge(u, v, i + M - 1, mat[u][v])
        add_edge(v, u, i + M - 1, mat[v][u])
    
    for i in range(Y + 1, 1 + Y + T):
        info_list = line_list[i].split(' ')
        channel = int(info_list[0])
        path_len = int(info_list[1])
        station_len = int(info_list[2])
        path_list = [int(info_list[i]) for i in range(3, 3 + path_len)]
        station_list = [int(info_list[i]) for i in range(3 + path_len, 3 + path_len + station_len)]
        ans += station_len * 100
        ans += path_len * 1
        res = check(i - Y - 1, channel, path_list, station_list)
        if res == -1:
            print(i - Y - 1)
            print('CHECK: ILLEGAL!')
            exit()

finally:
    output_object.close()


print('CHECK: OK!')

