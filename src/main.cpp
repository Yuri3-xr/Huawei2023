#include <bits/stdc++.h>

using i64 = std::int64_t;

constexpr int INF = std::numeric_limits<int>::max() / 2.0;
constexpr uint32_t SEED = 120;

std::mt19937 seed(SEED);

inline int randomInt(int l, int r) {
    // return Int [l,r);
    std::uniform_int_distribution<int> RNG(l, r - 1);
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
    std::vector<int> deg;
    int P;
    int cnt = 0;
    Graph(int n, int P)
        : adj(n), deg(n), P(P), mat(n, std::vector<int>(n, INF)){};
    inline void addEdge(int from, int to, int d) {
        adj[from].emplace_back(from, to, d, P, cnt);
        adj[to].emplace_back(to, from, d, P, cnt);

        mat[from][to] = mat[to][from] = std::min(mat[from][to], d);

        if (from > to) std::swap(from, to);
        edgeSet.push_back(
            {Edge(from, to, d, P, cnt),
             {(int)adj[from].size() - 1, (int)adj[to].size() - 1}});

        deg[from]++, deg[to]++;

        ++cnt;
        return;
    };
};

struct Task {
    int id;
    int from, to;
    int shortestPathLen;
    int channel = -1;
    std::vector<int> pathEdge = {};
    std::vector<int> pathNode = {};
    std::vector<int> station = {};
    std::vector<int> dis = {};
    int totalMinChannel = 0;
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

int solveTaskDistance(const Graph& G, int from, int to) {
    std::vector<int> dis(N, -1);

    std::queue<int> Q;
    Q.push(from);
    dis[from] = 0;

    while (!Q.empty()) {
        auto curNode = Q.front();
        Q.pop();
        for (const auto& edge : G.adj[curNode]) {
            if (dis[edge.to] == -1) {
                dis[edge.to] = dis[curNode] + 1;
                Q.push(edge.to);
            }
        }
    }

    return dis[to];
}

struct BFSNode {
    int cntChannel, nextDeg, nextNode;
    friend bool operator<(const BFSNode& opA, const BFSNode& opB) {
        return std::make_pair(opA.cntChannel, -opA.nextDeg) <
               std::make_pair(opB.cntChannel, -opB.nextDeg);
    }
};
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
        std::vector<BFSNode> nextNodeList;
        for (const auto& edge : G.adj[curNode]) {
            if (!vis[edge.to]) {
                nextNodeList.push_back(
                    {edge.cntChannel, G.deg[edge, to], edge.to});
                last[edge.to] = curNode;
            }
        }
        // 启发式，剩余信道少的边先拓展
        std::sort(begin(nextNodeList), end(nextNodeList));
        for (const auto& [cnt, nextDeg, node] : nextNodeList) {
            Q.push(node);
            vis[node] = 1;
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

struct UnionFind {
    std::vector<int> data;
    UnionFind(int N) : data(N, -1) {}

    int find(int k) { return data[k] < 0 ? k : data[k] = find(data[k]); }

    int unite(int x, int y) {
        if ((x = find(x)) == (y = find(y))) return false;
        if (data[x] > data[y]) std::swap(x, y);
        data[x] += data[y];
        data[y] = x;
        return true;
    }

    int size(int k) { return -data[find(k)]; }

    int same(int x, int y) { return find(x) == find(y); }
};

std::vector<std::pair<int, int>> singleChannelBfs(
    const Graph& G,
    const std::vector<std::pair<Edge, std::pair<int, int>>>& sortedEdgeSet,
    int from, int to, int p) {
    // 注意这里返回的是前驱数组，如果不可道，则返回空
    std::vector<int> mstEdgeSeg(G.cnt);  // 桶记录当前边在不在MST中
    UnionFind uf(N);

    for (const auto& [edge, _] : sortedEdgeSet) {
        if (edge.markChannel[p] != -1) continue;
        if (!uf.same(edge.from, edge.to)) {
            uf.unite(edge.from, edge.to);
            mstEdgeSeg[edge.id] = 1;
        }
        if (uf.same(from, to)) break;
    }
    if (!uf.same(from, to)) return {};

    std::vector<int> vis(N, 0);
    std::vector<std::pair<int, int>> last(N, {-1, -1});

    std::queue<int> Q;
    Q.push(from);

    while (!Q.empty()) {
        auto curNode = Q.front();
        vis[curNode] = 1;
        Q.pop();
        if (curNode == to) break;
        for (const auto& edge : G.adj[curNode]) {
            if (mstEdgeSeg[edge.id] == 0) continue;
            if (!vis[edge.to]) {
                last[edge.to] = {curNode, edge.id};
                vis[edge.to] = 1;
                Q.push(edge.to);
            }
        }
    }

    if (vis[to] == 0) return {};
    return last;
}

void solveSingleTask(Graph& G, Task& task) {
    auto curEdgeSet = G.edgeSet;
    std::sort(begin(curEdgeSet), end(curEdgeSet),
              [&](const auto& A, const auto& B) {
                  return A.first.cntChannel > B.first.cntChannel;
              });

    for (int p = 0; p < P; p++) {
        auto curLast = singleChannelBfs(G, curEdgeSet, task.from, task.to, p);
        if (curLast.empty()) continue;

        int curNode = task.to;
        std::vector<int> resPathEdge, resPathNode, resDis;
        int minChannel = INF;

        while (curNode != -1) {
            resPathNode.push_back(curNode);
            auto [prevNode, prevEdge] = curLast[curNode];
            curNode = prevNode;
            if (prevEdge != -1) {
                resPathEdge.push_back(prevEdge);
                resDis.push_back(G.edgeSet[prevEdge].first.distance);
                minChannel =
                    std::min(minChannel, G.edgeSet[prevEdge].first.cntChannel);
            }
        }
        std::reverse(begin(resPathEdge), end(resPathEdge));
        std::reverse(begin(resPathNode), end(resPathNode));
        std::reverse(begin(resDis), end(resDis));

        if (minChannel > task.totalMinChannel) {
            task.totalMinChannel = minChannel;
            task.pathEdge = std::move(resPathEdge);
            task.pathNode = std::move(resPathNode);
            task.dis = std::move(resDis);
            task.channel = p;
        }
    }

    if (task.pathNode.empty()) {
        auto addPath = newBfs(G, task.from, task.to);
        for (int i = 0; i < (int)addPath.size() - 1; i++) {
            G.addEdge(addPath[i], addPath[i + 1],
                      G.mat[addPath[i]][addPath[i + 1]]);
        }
        auto curEdgeSet = G.edgeSet;
        std::sort(begin(curEdgeSet), end(curEdgeSet),
                  [&](const auto& A, const auto& B) {
                      return A.first.cntChannel > B.first.cntChannel;
                  });
        int randChannel = randomInt(0, P);
        auto curLast =
            singleChannelBfs(G, curEdgeSet, task.from, task.to, randChannel);

        int curNode = task.to;
        std::vector<int> resPathEdge, resPathNode, resDis;
        int minChannel = 0;
        while (curNode != -1) {
            resPathNode.push_back(curNode);
            auto [prevNode, prevEdge] = curLast[curNode];
            curNode = prevNode;
            if (prevEdge != -1) {
                resPathEdge.push_back(prevEdge);
                resDis.push_back(G.edgeSet[prevEdge].first.distance);
                minChannel =
                    std::min(minChannel, G.edgeSet[prevEdge].first.cntChannel);
            }
        }
        std::reverse(begin(resPathEdge), end(resPathEdge));
        std::reverse(begin(resPathNode), end(resPathNode));
        std::reverse(begin(resDis), end(resDis));

        task.totalMinChannel = minChannel;
        task.pathEdge = std::move(resPathEdge);
        task.pathNode = std::move(resPathNode);
        task.dis = std::move(resDis);
        task.channel = randChannel;
    }

    int finalChannel = task.channel;
    for (auto edgeId : task.pathEdge) {
        G.edgeSet[edgeId].first.markChannel[finalChannel] = task.id;
        G.edgeSet[edgeId].first.cntChannel--;

        int _u = G.edgeSet[edgeId].first.from;
        int _v = G.edgeSet[edgeId].first.to;
        if (_u > _v) std::swap(_u, _v);

        int _uId = G.edgeSet[edgeId].second.first;
        int _vId = G.edgeSet[edgeId].second.second;

        G.adj[_u][_uId].markChannel[finalChannel] = task.id;
        G.adj[_u][_uId].cntChannel--;
        G.adj[_v][_vId].markChannel[finalChannel] = task.id;
        G.adj[_v][_vId].cntChannel--;
    }

    int nowDis = 0;
    for (int i = 0; i < (int)task.pathNode.size() - 1; i++) {
        if (nowDis + task.dis[i] > D) {
            task.station.push_back(task.pathNode[i]);
            nowDis = 0;
        }
        nowDis += task.dis[i];
    }

    return;
}

void solveAllTask(Graph& G) {
    for (int i = 0; i < T; i++) {
        solveSingleTask(G, taskList[i]);
    }
    return;
}

void outputAnswer(const Graph& G) {
    std::cout << G.cnt - M << "\n";

    std::vector<std::tuple<int, int, int>> edgeList;
    for (int from = 0; from < N; ++from) {
        for (auto e : G.adj[from]) {
            int to = e.to;
            if (e.id < M) continue;
            if (from > to) continue;
            edgeList.emplace_back(e.id, from, to);
        }
    }

    std::sort(begin(edgeList), end(edgeList));
    for (auto [id, u, v] : edgeList) {
        std::cout << u << " " << v << '\n';
    }

    std::vector<int> taskIndex(T);
    for (int i = 0; i < T; ++i) {
        taskIndex[taskList[i].id] = i;
    }
    for (int pt = 0; pt < T; ++pt) {
        int i = taskIndex[pt];
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
        auto _dis = solveTaskDistance(G, _from, _to);
        taskList.push_back({i, _from, _to, _dis});
    }

    std::sort(begin(taskList), end(taskList), [&](auto cmpA, auto cmpB) {
        return cmpA.shortestPathLen > cmpB.shortestPathLen;
    });

    solveAllTask(G);
    outputAnswer(G);

    return 0;
}