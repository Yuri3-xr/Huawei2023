#include <bits/stdc++.h>

using i64 = std::int64_t;

struct Edge {
    int id;  // 边的标号
    int to;
    int distance;                  // 距离
    std::vector<int> markChannel;  //  表示i信道被markChannel[i]任务是用了
    int cntChannel;                // 记录当前剩余的信道数量
    Edge() = delete;
    Edge(int P) : markChannel(P, -1), cntChannel(P){};
    Edge(int to, int distance, int P, int id)
        : to(to),
          distance(distance),
          markChannel(P, -1),
          cntChannel(P),
          id(id){};
};

struct Graph {
    std::vector<std::vector<Edge>> adj;
    int P;
    int cnt = 0;
    Graph(int n, int P) : adj(n), P(P){};
    inline void addEdge(int from, int to, int d) {
        adj[from].emplace_back(to, d, P, cnt++);
        adj[to].emplace_back(from, d, P, cnt++);
        return;
    };
};

int N, M, T, P, D;
std::vector<std::pair<int, int>> task;

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
        task.emplace_back(_from, _to);
    }
    return 0;
}