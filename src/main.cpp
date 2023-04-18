#include <bits/stdc++.h>

using i64 = std::int64_t;

constexpr int INF = std::numeric_limits<int>::max() / 2.0;
constexpr uint32_t SEED = 120;

constexpr int EDGE_WEIGHT = 1;
constexpr int ADDED_EDGE_WEIGHT = 1e6;

std::mt19937 seed(SEED);
struct Timer {
    std::chrono::high_resolution_clock::time_point st;

    Timer() { reset(); }

    void reset() { st = std::chrono::high_resolution_clock::now(); }

    std::chrono::milliseconds::rep elapsed() {
        auto ed = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::milliseconds>(ed - st)
            .count();
    }
};

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
    std::vector<std::unordered_map<int, int>> mat;
    std::vector<std::pair<Edge, std::pair<int, int>>>
        edgeSet;  // 第二个pair表示u,v(u<v) {u的位置，v的位置}
    int P;
    int cnt = 0;
    Graph(int n, int P) : adj(n), P(P), mat(n){};
    inline void addEdge(int from, int to, int d) {
        adj[from].emplace_back(from, to, d, P, cnt);
        adj[to].emplace_back(to, from, d, P, cnt);

        if (mat[from].find(to) == mat[from].end()) mat[from].insert({to, INF});
        if (mat[to].find(from) == mat[to].end()) mat[to].insert({from, INF});

        mat[from][to] = mat[to][from] = std::min(mat[from][to], d);

        if (from > to) std::swap(from, to);
        edgeSet.push_back(
            {Edge(from, to, d, P, cnt),
             {(int)adj[from].size() - 1, (int)adj[to].size() - 1}});

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
    int totalSumDistance = 0;
    bool addEdgeFlag = false;
};

/*
N: 点数
M: 边数
T: 任务数
P: 信道上限
D: 衰减上限
*/
int N, M, T, P, D;

//====================================================

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
    int cntChannel, distance, nextNode;
    friend bool operator<(const BFSNode& opA, const BFSNode& opB) {
        return std::make_pair(opA.cntChannel, opA.distance) <
               std::make_pair(opB.cntChannel, opB.distance);
    }
};
std::pair<std::vector<int>, std::vector<int>> newBfs(const Graph& G, int from,
                                                     int to) {
    // 若该两点不连通，则需要bfs出一条最短路，然后把这条最短路上的边重新加入一遍
    // 启发式去做bfs，因为这些边需要被加，所以每次拓展时，按照当前的剩余信道数排序拓展
    // 返回新加的边的点序列
    std::vector<int> vis(N, 0);
    std::vector<std::pair<int, int>> last(N, {-1, -1});

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
                    {edge.cntChannel, edge.distance, edge.to});
                last[edge.to].first = curNode;
                last[edge.to].second = edge.id;
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
    std::vector<int> pathNode;
    std::vector<int> pathEdge;
    while (pNode != -1) {
        pathNode.push_back(pNode);
        if (last[pNode].second != -1) {
            pathEdge.push_back(last[pNode].second);
        }
        pNode = last[pNode].first;
    }
    std::reverse(begin(pathNode), end(pathNode));
    std::reverse(begin(pathEdge), end(pathEdge));

    return std::make_pair(pathNode, pathEdge);
}

//==================================================== dijkstra

inline int cost(int x) {
    if (M > N * 30) {
        return exp(1.0 * x / 20);
    } else {
        return x;
    }
}

std::vector<std::pair<int, int>> singleChannelDijkstra(const Graph& G, int from,
                                                       int to, int p,
                                                       int& toDis) {
    std::vector<int> vis(N, 0);
    std::vector<std::pair<int, int>> last(N, {-1, -1});
    std::vector<int> dis(N, INF);

    dis[from] = 0;
    std::priority_queue<std::pair<int, int>> Q;
    Q.push(std::make_pair(-dis[from], from));

    while (!Q.empty()) {
        auto [negDis, u] = Q.top();
        Q.pop();
        dis[u] = -negDis;

        if (vis[u]) continue;
        vis[u] = true;

        if (u == to) break;

        for (const auto& edge : G.adj[u]) {
            int v = edge.to;
            if (vis[v]) continue;
            if (edge.cntChannel >= P && edge.id >= M) continue;  // deleted
            int nxtDis = dis[u], lastEdgeId = -1;
            if (edge.markChannel[p] == -1) {
                nxtDis += EDGE_WEIGHT + cost(P - edge.cntChannel);
                lastEdgeId = edge.id;
            } else {
                nxtDis += ADDED_EDGE_WEIGHT;
            }
            if (nxtDis < dis[v]) {
                dis[v] = nxtDis;
                last[v] = {u, lastEdgeId};
                Q.push(std::make_pair(-dis[v], v));
            }
        }
    }

    toDis = dis[to];
    if (vis[to] == 0) return {};
    return last;
}

void solveSingleTask(Graph& G, Task& task) {
    int toDis = INF, bestChannel = -1;
    for (int p = 0; p < P; ++p) {
        int nowToDis = INF;
        auto curLast =
            singleChannelDijkstra(G, task.from, task.to, p, nowToDis);
        if (nowToDis < toDis) {
            toDis = nowToDis;
            bestChannel = p;
        }
    }

    assert(toDis < INF);
    auto curLast =
        singleChannelDijkstra(G, task.from, task.to, bestChannel, toDis);
    int curNode = task.to;
    std::vector<int> resPathEdge, resPathNode, resDis;
    int minChannel = 0;
    int sumDistance = 0;

    while (curNode != -1) {
        resPathNode.push_back(curNode);
        auto [prevNode, prevEdge] = curLast[curNode];
        if (prevNode == -1) break;
        if (prevEdge == -1) {
            prevEdge = G.cnt;
            task.addEdgeFlag = true;
            int u = curNode, v = prevNode;
            G.addEdge(u, v, G.mat[u][v]);
        }
        curNode = prevNode;
        resPathEdge.push_back(prevEdge);
        resDis.push_back(G.edgeSet[prevEdge].first.distance);
        minChannel = std::min(minChannel, G.edgeSet[prevEdge].first.cntChannel);
        // sumDistance += G.edgeSet[prevEdge].first.cntChannel;
    }

    std::reverse(begin(resPathEdge), end(resPathEdge));
    std::reverse(begin(resPathNode), end(resPathNode));
    std::reverse(begin(resDis), end(resDis));

    task.totalMinChannel = minChannel;
    task.totalSumDistance = toDis;
    task.pathEdge = std::move(resPathEdge);
    task.pathNode = std::move(resPathNode);
    task.dis = std::move(resDis);
    task.channel = bestChannel;

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
}

void fix(Graph& G, Task& task) {
    if (!task.addEdgeFlag) return;

    int toDis = INF, bestChannel = -1;
    for (int p = 0; p < P; ++p) {
        int nowToDis = INF;
        auto curLast =
            singleChannelDijkstra(G, task.from, task.to, p, nowToDis);
        if (nowToDis < toDis) {
            toDis = nowToDis;
            bestChannel = p;
        }
    }

    if (toDis > task.totalSumDistance - ADDED_EDGE_WEIGHT * 0.8) {
        return;
    }

    int lastChannel = task.channel;
    for (auto edgeId : task.pathEdge) {
        G.edgeSet[edgeId].first.markChannel[lastChannel] = -1;
        G.edgeSet[edgeId].first.cntChannel++;

        int _u = G.edgeSet[edgeId].first.from;
        int _v = G.edgeSet[edgeId].first.to;
        if (_u > _v) std::swap(_u, _v);

        int _uId = G.edgeSet[edgeId].second.first;
        int _vId = G.edgeSet[edgeId].second.second;

        G.adj[_u][_uId].markChannel[lastChannel] = -1;
        G.adj[_u][_uId].cntChannel++;
        G.adj[_v][_vId].markChannel[lastChannel] = -1;
        G.adj[_v][_vId].cntChannel++;
    }

    task.pathEdge.clear();
    task.pathNode.clear();
    task.station.clear();
    task.dis.clear();

    solveSingleTask(G, task);
}

void solveAllTask(Graph& G, std::vector<Task>& taskList) {
    for (int i = 0; i < T; i++) {
        solveSingleTask(G, taskList[i]);
    }

    // #ifdef dijkstra
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < T; ++i) {
            fix(G, taskList[i]);
        }
    }
    // #endif
    return;
}

i64 solveScoreTask(const Graph& G, const std::vector<Task>& taskList) {
    i64 ret = 0;
    std::vector<std::tuple<int, int, int>> edgeList;
    for (int from = 0; from < N; ++from) {
        for (auto e : G.adj[from]) {
            int to = e.to;
            if (e.id < M) continue;
            if (from > to) continue;
            if (e.cntChannel >= P && e.id >= M) continue;
            edgeList.emplace_back(e.id, from, to);
        }
    }
    ret += 1LL * edgeList.size() * 1E6;
    for (const auto& p : taskList) {
        ret += 1LL * p.pathEdge.size() * 1;
        ret += 100LL * p.station.size();
    }
    return ret;
}

int cnt[100000];
void outputAnswer(const Graph& G, std::vector<Task>& taskList) {
    // std::cout << G.cnt - M << "\n";

    std::vector<std::tuple<int, int, int>> edgeList;
    for (int from = 0; from < N; ++from) {
        for (auto e : G.adj[from]) {
            int to = e.to;
            if (e.id < M) continue;
            if (from > to) continue;
            if (e.cntChannel >= P && e.id >= M) {
                ++cnt[e.id];
                continue;
            }
            edgeList.emplace_back(e.id, from, to);
        }
    }
    for (int i = 0; i < 100000; ++i) {
        cnt[i] += cnt[i - 1];
    }

    std::cout << edgeList.size() << "\n";
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
            std::cout << " " << edgeId - cnt[edgeId];
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
    Timer timer;
    timer.reset();
    std::cin >> N >> M >> T >> P >> D;
    Graph G(N, P);
    for (int i = 0; i < M; i++) {
        int _u, _v, _d;
        std::cin >> _u >> _v >> _d;
        G.addEdge(_u, _v, _d);
    }

    std::vector<Task> taskList;
    for (int i = 0; i < T; i++) {
        int _from, _to;
        std::cin >> _from >> _to;
        auto _dis = solveTaskDistance(G, _from, _to);
        taskList.push_back({i, _from, _to, _dis});
    }

    // constexpr int timeIWantToTry = 30;
    // std::vector<Graph> GList(timeIWantToTry, G);
    // std::vector<std::vector<Task>> listTaskList(timeIWantToTry, taskList);
    // std::vector<unsigned int> seedList;
    // for (int i = 1; i <= timeIWantToTry; i++) {
    //     seedList.push_back(i);
    // }

    // int bestSeed = 0, bestScore = INF;

    // for (int i = 0; i < timeIWantToTry; i++) {
    //     std::srand(seedList[i]);
    //     // std::random_shuffle(begin(listTaskList[i]), end(listTaskList[i]));
    //     std::random_shuffle(begin(listTaskList[i]), end(listTaskList[i]));
    //     for (auto& p : listTaskList[i]) {
    //         p.shortestPathLen += rand() % 4;
    //     }

    //     std::sort(begin(listTaskList[i]), end(listTaskList[i]),
    //               [&](auto cmpA, auto cmpB) {
    //                   return cmpA.shortestPathLen > cmpB.shortestPathLen;
    //               });

    //     solveAllTask(GList[i], listTaskList[i]);

    //     auto curScore = solveScoreTask(GList[i], listTaskList[i]);
    //     // std::cerr << i << " " << curScore << std::endl;
    //     if (curScore < bestScore) {
    //         bestScore = curScore;
    //         bestSeed = i;
    //     }
    //     unsigned int singleRunTime = (timer.elapsed() + i) / (i + 1);

    //     if (timer.elapsed() + singleRunTime >= 1000 * 110) break;
    // }

    // // std::cerr << bestSeed << std::endl;

    // // std::cerr << bestScore << std::endl;
    // outputAnswer(GList[bestSeed], listTaskList[bestSeed]);

    // std::cerr << timer.elapsed() << "ms" << std::endl;

    return 0;
}