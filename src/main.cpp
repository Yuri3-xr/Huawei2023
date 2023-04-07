#include <bits/stdc++.h>

using i64 = std::int64_t;

constexpr int INF = std::numeric_limits<int>::max() / 2.0;
constexpr uint32_t SEED = 120;

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
};

/*
N: 点数
M: 边数
T: 任务数
P: 信道上限
D: 衰减上限
*/
int N, M, T, P, D;

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

std::vector<std::pair<int, int>> singleChannelBfs(const Graph& G, int from,
                                                  int to, int p) {
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
        Q.pop();
        if (curNode == to) break;
        std::vector<BFSNode> nextNodeList;

        for (const auto& edge : G.adj[curNode]) {
            if (edge.markChannel[p] != -1) continue;
            if (!vis[edge.to]) {
                nextNodeList.push_back(
                    {-edge.cntChannel, edge.distance, edge.to});
                last[edge.to] = {curNode, edge.id};
            }
        }
        // 启发式，剩余信道多的边先拓展
        std::sort(begin(nextNodeList), end(nextNodeList));
        for (const auto& [cnt, nextDeg, node] : nextNodeList) {
            Q.push(node);
            vis[node] = 1;
        }
    }

    if (vis[to] == 0) return {};
    return last;
}

int bestChannel(Graph& G, Task& task, std::vector<int>& pathNode,
                std::vector<int>& pathEdge) {
    int bestChannel = -1, bestCnt = -1;
    for (int p = 0; p < P; ++p) {
        int cnt = 0;
        for (int i = 0; i < (int)pathEdge.size(); ++i) {
            auto edge = G.edgeSet[pathEdge[i]].first;
            if (edge.markChannel[p] != -1) {
                ++cnt;
            }
        }
        if (bestChannel == -1) {
            bestChannel = p;
            bestCnt = cnt;
        } else if (bestCnt > cnt) {
            bestChannel = p;
            bestCnt = cnt;
        }
    }
    for (int i = 0; i < (int)pathEdge.size(); ++i) {
        int u = pathNode[i], v = pathNode[i + 1];
        auto edge = G.edgeSet[pathEdge[i]].first;
        if (edge.markChannel[bestChannel] != -1) {
            G.addEdge(u, v, G.mat[u][v]);
        }
    }
    return bestChannel;
}

void solveSingleTask(Graph& G, Task& task) {
    for (int p = 0; p < P; p++) {
        auto curLast = singleChannelBfs(G, task.from, task.to, p);
        if (curLast.empty()) continue;

        int curNode = task.to;
        std::vector<int> resPathEdge, resPathNode, resDis;
        int minChannel = INF;  // ？？
        int sumDistance = 0;
        int cntEdge = 0;
        while (curNode != -1) {
            resPathNode.push_back(curNode);
            auto [prevNode, prevEdge] = curLast[curNode];
            curNode = prevNode;
            if (prevEdge != -1) {
                resPathEdge.push_back(prevEdge);
                resDis.push_back(G.edgeSet[prevEdge].first.distance);
                minChannel +=
                    (minChannel, G.edgeSet[prevEdge].first.cntChannel);
                sumDistance += G.edgeSet[prevEdge].first.distance;
                cntEdge++;
            }
        }
        minChannel = (minChannel + cntEdge + 1) / cntEdge;
        std::reverse(begin(resPathEdge), end(resPathEdge));
        std::reverse(begin(resPathNode), end(resPathNode));
        std::reverse(begin(resDis), end(resDis));

        if (std::make_pair(minChannel, -sumDistance) >
            std::make_pair(task.totalMinChannel, -task.totalSumDistance)) {
            task.totalMinChannel = minChannel;
            task.totalSumDistance = sumDistance;
            task.pathEdge = std::move(resPathEdge);
            task.pathNode = std::move(resPathNode);
            task.dis = std::move(resDis);
            task.channel = p;
        }
    }

    if (task.pathNode.empty()) {
        auto [addPathNode, addPathEdge] = newBfs(G, task.from, task.to);

        int trueChannel = bestChannel(
            G, task, addPathNode,
            addPathEdge);  // this variable's name should be modified.

        auto curLast = singleChannelBfs(G, task.from, task.to, trueChannel);

        int curNode = task.to;
        std::vector<int> resPathEdge, resPathNode, resDis;
        int minChannel = 0;
        int sumDistance = 0;

        while (curNode != -1) {
            resPathNode.push_back(curNode);
            auto [prevNode, prevEdge] = curLast[curNode];
            curNode = prevNode;
            if (prevEdge != -1) {
                resPathEdge.push_back(prevEdge);
                resDis.push_back(G.edgeSet[prevEdge].first.distance);
                minChannel =
                    std::min(minChannel, G.edgeSet[prevEdge].first.cntChannel);
                sumDistance += G.edgeSet[prevEdge].first.cntChannel;
            }
        }

        std::reverse(begin(resPathEdge), end(resPathEdge));
        std::reverse(begin(resPathNode), end(resPathNode));
        std::reverse(begin(resDis), end(resDis));

        task.totalMinChannel = minChannel;
        task.totalSumDistance = sumDistance;
        task.pathEdge = std::move(resPathEdge);
        task.pathNode = std::move(resPathNode);
        task.dis = std::move(resDis);
        task.channel = trueChannel;
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

void solveAllTask(Graph& G, std::vector<Task>& taskList) {
    for (int i = 0; i < T; i++) {
        solveSingleTask(G, taskList[i]);
    }
    return;
}

i64 solveScoreTask(const Graph& G, const std::vector<Task>& taskList) {
    i64 ret = 0;
    ret += 1LL * (G.cnt - M) * 1E6;
    for (const auto& p : taskList) {
        ret += 1LL * p.pathEdge.size() * 1;
        ret += 100LL * p.station.size();
    }
    return ret;
}
void outputAnswer(const Graph& G, std::vector<Task>& taskList) {
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

    constexpr int timeIWantToTry = 30;
    std::vector<Graph> GList(timeIWantToTry, G);
    std::vector<std::vector<Task>> listTaskList(timeIWantToTry, taskList);
    std::vector<unsigned int> seedList;
    for (int i = 1; i <= timeIWantToTry; i++) {
        seedList.push_back(i);
    }

    int bestSeed = 0, bestScore = INF;

    for (int i = 0; i < timeIWantToTry; i++) {
        std::srand(seedList[i]);
        // std::random_shuffle(begin(listTaskList[i]), end(listTaskList[i]));
        std::random_shuffle(begin(listTaskList[i]), end(listTaskList[i]));
        for (auto& p : listTaskList[i]) {
            p.shortestPathLen += rand() % 4;
        }

        std::sort(begin(listTaskList[i]), end(listTaskList[i]),
                  [&](auto cmpA, auto cmpB) {
                      return cmpA.shortestPathLen > cmpB.shortestPathLen;
                  });

        solveAllTask(GList[i], listTaskList[i]);

        auto curScore = solveScoreTask(GList[i], listTaskList[i]);
        std::cerr << i << " " << curScore << std::endl;
        if (curScore < bestScore) {
            bestScore = curScore;
            bestSeed = i;
        }
        unsigned int singleRunTime = (timer.elapsed() + i) / (i + 1);

        if (timer.elapsed() + singleRunTime >= 1000 * 110) break;
    }

    // std::cerr << bestScore << std::endl;
    outputAnswer(GList[bestSeed], listTaskList[bestSeed]);

    // std::cerr << timer.elapsed() << "ms" << std::endl;

    return 0;
}