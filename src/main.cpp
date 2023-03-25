#include <bits/stdc++.h>

using i64 = std::int64_t;

constexpr int INF = std::numeric_limits<int>::max()/2.0;
constexpr uint32_t SEED = 229;

std::mt19937 seed(SEED);

inline int randomInt(int l,int r){
    // return Int [l,r);
    std::uniform_int_distribution<int> RNG(l, r + 1);
    return RNG(seed);
}

struct Edge {
    int id;  // 边的标号
    int from;
    int to;
    int distance;                  // 距离
    std::vector<int> markChannel;  //  表示i信道被markChannel[i]任务是用了
    int cntChannel;                // 记录当前剩余的信道数量
    Edge() = delete;
    Edge(int P) : markChannel(P, -1), cntChannel(P){};
    Edge(int from, int to, int distance, int P, int id)
        : from(from),
          to(to),
          distance(distance),
          markChannel(P, -1),
          cntChannel(P),
          id(id){};
};

struct Graph {
    std::vector<std::vector<Edge>> adj;
    std::vector<std::vector<int>> mat;
    int P;
    int cnt = 0;
    Graph(int n, int P) : adj(n), P(P), mat(n,std::vector<int>(n, INF)) {};
    inline void addEdge(int from, int to, int d) {
        adj[from].emplace_back(from, to, d, P, cnt);
        adj[to].emplace_back(to, from, d, P, cnt);
        mat[from][to] = mat[to][from] = std::min(mat[from][to], d);
        ++cnt;
        return;
    };
};

struct Task {
    int id;
    int from, to;
    int channel = -1;
    std::vector<int> pathEdge = {};
    std::vector<int> pathNode = {};
    std::vector<int> station = {};
};

/*
N: 点数
M: 边数
T: 任务数
P: 信道上限
D: 衰减上限
*/
int N, M, T, P, D;
std::vector<Task> task;

int findChannel(Graph &G, int from) {
    std::vector<int> possibleChannel;
    for (int i = 0; i < P; ++i) {
        for (auto e : G.adj[from]) {
            if (e.markChannel[i] == -1) {
                possibleChannel.emplace_back(i);
                break;
            }
        }
    }
    if (!possibleChannel.size()) return -1;
    std::sort(possibleChannel.begin(), possibleChannel.end());
    possibleChannel.resize(std::unique(possibleChannel.begin(), possibleChannel.end()) - possibleChannel.begin());
    int randNum = randomInt(0, (int)possibleChannel.size());
    return possibleChannel[randNum];
}

bool randomAddEdge(Graph &G, int from) {
    if (!G.adj[from].size()) return false;
    int randNum = randomInt(0, (int)G.adj[from].size());
    for (auto e : G.adj[from]) {
        --randNum;
        if (randNum <= 0) {
            G.addEdge(from, e.to, G.mat[from][e.to]);
            return true;
        }
    }
    return false;
}

bool randomAddEdge(Graph &G, int from, int &ret) {
    if (!G.adj[from].size()) return false;
    int randNum = randomInt(0, (int)G.adj[from].size());
    for (auto e : G.adj[from]) {
        --randNum;
        if (randNum <= 0) {
            // G.addEdge(from, e.to, G.mat[from][e.to]);
            ret = e.to;
            return true;
        }
    }
    return false;
}

std::vector<int> bfs(Graph &G, Task &t) {
    std::vector<int > lst(N, -1);
    std::vector<int> vis(N, 0);
    std::queue<int> Q;
    while (t.channel == -1) {
        t.channel = findChannel(G, t.from); // 随机分配信道
        randomAddEdge(G, t.from); // 随机加边
    }
    Q.push(t.from);
    while (!Q.empty()) {
        int from = Q.front(); Q.pop();
        vis[from] = 1;
        // std::cerr << from << std::endl;
        if (from == t.to) break;
        for (auto e : G.adj[from]) {
            int to = e.to;
            if (vis[to]) continue;
            if (!vis[to]) Q.push(to);
            vis[to] = 1;
            lst[to] = from;
        }
    }
    return lst;
}

void bruteForce(Graph &G) {
    for (int i = 0; i < T; ++i) {
        auto lst = bfs(G, task[i]);
        int nowNode = task[i].to;
        std::cerr << task[i].from << " " << task[i].to << std::endl;
        std::vector<int> dis;
        while (nowNode != -1) {
            std::cerr << nowNode << std::endl;
            task[i].pathNode.push_back(nowNode);
            int lstNode = lst[nowNode];
            if (lstNode == -1) {
                break;
            }
            if (edgeId != -1) {
                task[i].pathEdge.push_back(edgeId);
            } else {
                G.addEdge(lstNode, nowNode, G.mat[lstNode][nowNode]);
                edgeId = G.cnt - 1;
            }
            for (auto &e : G.adj[nowNode]) {
                if (e.id == edgeId) {
                    e.markChannel[task[i].channel] = task[i].id;
                    dis.push_back(e.distance);
                }
            }
            for (auto &e : G.adj[lstNode]) {
                if (e.id == edgeId) {
                    e.markChannel[task[i].channel] = task[i].id;
                }
            }
            nowNode = lstNode;
        }
        std::reverse(task[i].pathEdge.begin(), task[i].pathEdge.end());
        std::reverse(task[i].pathNode.begin(), task[i].pathNode.end());
        std::reverse(dis.begin(), dis.end());
        int nowDis = 0;
        for (int j = 0; j < (int)dis.size() - 1; ++j) {
            if (nowDis + dis[j] > D) {
                task[i].station.push_back(task[i].pathNode[j]);
            }
        }
    }
}

void output(Graph &G) {
    std::cout << G.cnt - M << "\n";
    for (int from = 0; from < N; ++from) {
        for (auto e : G.adj[from]) {
            int to = e.to;
            if (e.id < M) continue;
            if (from > to) continue;
            std::cout << from << " " << to << "\n";
        }
    }
    for (int i = 0; i < T; ++i) {
        std::cout << task[i].channel << " " << task[i].pathEdge.size() << " " << task[i].station.size();
        for (auto edgeId : task[i].pathEdge) {
            std::cout << " " << edgeId;
        }
        for (auto stationId : task[i].station) {
            std::cout << " " << stationId;
        }
        std::cout << "\n";
    }
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    std::cin >> N >> M >> T >> P >> D;
    Graph G(N, P);
    for (int i = 0; i < M; i++) {
        int _u, _v, _d;
        std::cin >> _u >> _v >> _d;
        G.addEdge(_u, _v, _d);
    }

    for (int i = 0; i < T; i++) {
        int _from, _to;
        std::cin >> _from >> _to;
        task.push_back({i, _from, _to});
    }
    bruteForce(G);
    output(G);
    return 0;
}