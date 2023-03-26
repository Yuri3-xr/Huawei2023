#include <bits/stdc++.h>

using i64 = std::int64_t;

constexpr int INF = std::numeric_limits<int>::max() / 2.0;
constexpr uint32_t SEED = 229;

std::mt19937 seed(SEED);

inline int randomInt(int l, int r) {
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
    std::vector<std::pair<Edge, std::pair<int, int>>>
        edgeSet;  // 第二个pair表示u,v(u<v) {u的位置，v的位置}
    int P;
    int cnt = 0;
    Graph(int n, int P) : adj(n), P(P), mat(n, std::vector<int>(n, INF)){};
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
    int totalSumChannel = 0;
};

/*
N: 点数
M: 边数
T: 任务数
P: 信道上限
D: 衰减上限
*/
int N, M, T, P, D;
std::vector<Task> taskList;

std::vector<int> newBfs(const Graph& G, int from, int to) {
    // 若该两点不连通，则需要bfs出一条最短路，然后把这条最短路上的边重新加入一遍
    // 启发式去做bfs，因为这些边需要被加，所以每次拓展时，按照当前的剩余信道数排序拓展
    // 返回新加的边的点序列
    std::vector<int> vis(N, 0);
    std::vector<int> last(N, -1);

    std::queue<int> Q;
    Q.push(from);

    while (!Q.empty()) {
        auto curNode = Q.front();
        vis[curNode] = 1;
        Q.pop();
        if (curNode == to) break;
        std::vector<std::pair<int, int>> nextNodeList;
        for (const auto& edge : G.adj[curNode]) {
            if (!vis[edge.to]) {
                nextNodeList.emplace_back(edge.cntChannel, edge.to);
                last[edge.to] = curNode;
            }
        }
        // 启发式，剩余信道少的边先拓展
        std::sort(begin(nextNodeList), end(nextNodeList));
        for (const auto& [cnt, node] : nextNodeList) {
            Q.push(node);
        }
    }

    int pNode = to;
    std::vector<int> ret;
    while (pNode != -1) {
        ret.push_back(pNode);
        pNode = last[pNode];
    }
    std::reverse(begin(ret), end(ret));

    return ret;
}

std::vector<std::pair<int, int>> singleChannelBfs(const Graph& G, int from,
                                                  int to, int P) {
    // 对于P信道单独bfs
    // 同样是启发式bfs
    // 注意这里返回的是前驱数组，如果不可道，则返回空
    std::vector<int> vis(N, 0);
    std::vector<std::pair<int, int>> last(N, {-1, -1});

    std::queue<int> Q;
    Q.push(from);

    while (!Q.empty()) {
        auto curNode = Q.front();
        vis[curNode] = 1;
        Q.empty();
        if (curNode == to) break;
        std::vector<std::pair<int, int>> nextNodeList;
        for (const auto& edge : G.adj[curNode]) {
            if (edge.markChannel[P] != -1) continue;
            if (!vis[edge.to]) {
                nextNodeList.emplace_back(-edge.cntChannel, edge.to);
                last[edge.to] = {curNode, edge.id};
            }
        }
        // 启发式，剩余信道多的边先拓展
        std::sort(begin(nextNodeList), end(nextNodeList));
        for (const auto& [cnt, node] : nextNodeList) {
            Q.push(node);
        }
    }

    if (vis[to] == 0) return {};
    return last;
}

void solveSingleTask(Graph& G, Task& task) {
    for (int p = 0; p < P; p++) {
        auto curLast = singleChannelBfs(G, task.from, task.to, p);
        // if (curLast.empty()) {
        //     auto addPath = newBfs(G, task.from, task.to);
        //     for (int i = 0; i < (int)addPath.size() - 1; i++) {
        //         G.addEdge(addPath[i], addPath[i + 1],
        //                   G.mat[addPath[i]][addPath[i + 1]]);
        //     }
        // }

        int curNode = task.to;
        std::vector<int> resPathEdge, resNodeEdge;
        int sumChannel;
        while (curNode != -1) {
            resNodeEdge.push_back(curNode);
        }
    }
}

void outputAnswer(const Graph& G) {
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
        std::cout << taskList[i].channel << " "
                  << (int)taskList[i].pathEdge.size() << " "
                  << (int)taskList[i].station.size();
        for (auto edgeId : taskList[i].pathEdge) {
            std::cout << " " << edgeId;
        }
        for (auto stationId : taskList[i].station) {
            std::cout << " " << stationId;
        }
        std::cout << "\n";
    }

    return;
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
        taskList.push_back({i, _from, _to});
    }

    return 0;
}