#include <bits/stdc++.h>

using i64 = std::int64_t;

constexpr i64 INF = 1e18;
constexpr int INT_INF = std::numeric_limits<int>::max() / 2.0;
constexpr uint32_t SEED = 229;

constexpr int EDGE_WEIGHT = 10;
constexpr int ADDED_EDGE_WEIGHT = 1e6;
constexpr int PRE_ADDED_EDGE_WEIGHT = 9;

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
    int distance;  // 距离
    int hop;
    std::vector<int> markChannel;  //  表示i信道被markChannel[i]任务是用了
    int cntChannel;                // 记录当前剩余的信道数量
    bool deleted = false;
    bool preAdd = false;
    bool added = false;
    Edge() = delete;
    Edge(int P) : markChannel(P, -1), cntChannel(P){};
    Edge(int from, int to, int distance, int hop, int P, int id)
        : from(from),
          to(to),
          distance(distance),
          hop(hop),
          markChannel(P, -1),
          cntChannel(P),
          id(id){};
};

struct Graph {
    std::vector<std::vector<Edge>> adj;
    std::vector<std::unordered_map<int, std::tuple<int, int, int>>> mat;
    std::vector<std::pair<Edge, std::pair<int, int>>>
        edgeSet;  // 第二个pair表示u,v(u<v) {u的位置，v的位置}
    int P;
    int cnt = 0;
    Graph(int n, int P) : adj(n), P(P), mat(n){};
    inline void addEdge(int from, int to, int d, int h) {
        adj[from].emplace_back(from, to, d, h, P, cnt);
        adj[to].emplace_back(to, from, d, h, P, cnt);

        if (mat[from].find(to) == mat[from].end()) {
            mat[from].insert({to, {-1, INT_INF, INT_INF}});
        }
        if (mat[to].find(from) == mat[to].end()) {
            mat[to].insert({from, {-1, INT_INF, INT_INF}});
        }

        auto [preId, preD, preH] = mat[from][to];
        if (preD > d) {
            mat[from][to] = mat[to][from] = std::make_tuple(cnt, d, h);
        }

        if (from > to) std::swap(from, to);
        edgeSet.push_back(
            {Edge(from, to, d, h, P, cnt),
             {(int)adj[from].size() - 1, (int)adj[to].size() - 1}});

        ++cnt;
        return;
    };
};

struct Task {
    int id;
    int from, to;
    int shortestPathLen;
    int r;
    int channel = -1;
    std::vector<int> pathEdge = {};
    std::vector<int> pathNode = {};
    std::vector<int> station = {};
    std::vector<int> dis = {};
    std::vector<int> hop = {};
    int totalMinChannel = 0;
    i64 totalSumDistance = 0;
    int addEdge = 0;
};

/*
N: 点数
M: 边数
T: 任务数
P: 信道上限
D: 衰减上限
*/
int N, M, T, R, P, D, H;
int cnt[100000];

//====================================================

int solveTaskDistance(const Graph &G, int from, int to) {
    std::vector<int> dis(N, -1);

    std::queue<int> Q;
    Q.push(from);
    dis[from] = 0;

    while (!Q.empty()) {
        auto curNode = Q.front();
        Q.pop();
        for (const auto &edge : G.adj[curNode]) {
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
    friend bool operator<(const BFSNode &opA, const BFSNode &opB) {
        return std::make_pair(opA.cntChannel, opA.distance) <
               std::make_pair(opB.cntChannel, opB.distance);
    }
};
std::pair<std::vector<int>, std::vector<int>> newBfs(const Graph &G, int from,
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
        for (const auto &edge : G.adj[curNode]) {
            if (!vis[edge.to]) {
                nextNodeList.push_back(
                    {edge.cntChannel, edge.distance, edge.to});
                last[edge.to].first = curNode;
                last[edge.to].second = edge.id;
            }
        }
        // 启发式，剩余信道少的边先拓展
        std::sort(begin(nextNodeList), end(nextNodeList));
        for (const auto &[cnt, nextDeg, node] : nextNodeList) {
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

inline int cost(int x) { return x; }

namespace RadixHeap {
template <typename Key, typename Val>
struct RadixHeap {
    using uint = typename std::make_unsigned<Key>::type;
    static constexpr int bit = sizeof(Key) * 8;
    std::array<std::vector<std::pair<uint, Val>>, bit + 1> vs;
    std::array<uint, bit + 1> ms;

    int s;
    uint last;

    RadixHeap() : s(0), last(0) { std::fill(begin(ms), end(ms), uint(-1)); }

    bool empty() const { return s == 0; }

    int size() const { return s; }

    inline uint64_t getbit(uint a) const { return 64 - __builtin_clzll(a); }

    void push(const uint &key, const Val &val) {
        s++;
        uint64_t b = getbit(key ^ last);
        vs[b].emplace_back(key, val);
        ms[b] = std::min(key, ms[b]);
    }

    std::pair<uint, Val> pop() {
        if (ms[0] == uint(-1)) {
            int idx = 1;
            while (ms[idx] == uint(-1)) idx++;
            last = ms[idx];
            for (auto &p : vs[idx]) {
                uint64_t b = getbit(p.first ^ last);
                vs[b].emplace_back(p);
                ms[b] = std::min(p.first, ms[b]);
            }
            vs[idx].clear();
            ms[idx] = uint(-1);
        }
        --s;
        auto res = vs[0].back();
        vs[0].pop_back();
        if (vs[0].empty()) ms[0] = uint(-1);
        return res;
    }
};
}  // namespace RadixHeap
std::vector<std::pair<int, int>> singleChannelDijkstra(
    const Graph &G, int from, int to, int p, int taskId, i64 &toDis,
    std::vector<std::unordered_map<int, int>> &edgeMapList,
    bool fixFlag = false) {
    std::vector<int> vis(N, 0);
    std::vector<std::pair<int, int>> last(N, {-1, -1});
    std::vector<i64> dis(N, INF);

    dis[from] = 0;
    // RadixHeap::RadixHeap<int, int> Q;
    std::priority_queue<std::pair<i64, int>> Q;
    Q.push({-dis[from], from});

    while (!Q.empty()) {
        auto [negDis, u] = Q.top();
        Q.pop();
        dis[u] = -negDis;

        if (vis[u]) continue;
        vis[u] = true;

        if (u == to) break;

        for (const auto &edge : G.adj[u]) {
            int v = edge.to;
            if (vis[v]) continue;
            int nxtDis = dis[u], lastEdgeId = -1;
            if (edge.deleted) continue;  // has been removed.
            if (fixFlag) {
                if (edge.markChannel[p] != -1)
                    continue;  // this channel has been used.
            }
            if (edge.markChannel[p] == -1 &&
                edgeMapList[taskId].find(edge.id) ==
                    edgeMapList[taskId].end()) {
                nxtDis += cost(P - edge.cntChannel);
                if (edge.preAdd) {
                    nxtDis += PRE_ADDED_EDGE_WEIGHT;
                } else {
                    nxtDis += EDGE_WEIGHT;
                }
                lastEdgeId = edge.id;
            } else {
                nxtDis += ADDED_EDGE_WEIGHT;
            }
            if (nxtDis < dis[v]) {
                dis[v] = nxtDis;
                last[v] = {u, lastEdgeId};
                Q.push({-dis[v], v});
            }
        }
    }

    toDis = dis[to];
    if (vis[to] == 0) return {};
    return last;
}

void completeTask(Graph &G, Task &task,
                  std::vector<std::pair<int, int>> &curLast, i64 toDis,
                  int bestChannel,
                  std::vector<std::unordered_map<int, int>> &edgeMapList) {
    int curNode = task.to;
    std::vector<int> resPathEdge, resPathNode, resDis, resHop;
    int minChannel = 0;
    int sumDistance = 0;

    task.addEdge = 0;
    while (curNode != -1) {
        resPathNode.push_back(curNode);
        auto [prevNode, prevEdge] = curLast[curNode];
        if (prevNode == -1) break;
        if (prevEdge == -1) {
            prevEdge = G.cnt;
            ++task.addEdge;
            int u = curNode, v = prevNode;
            auto [id, dis, hop] = G.mat[u][v];
            G.addEdge(u, v, dis, hop);
        }
        curNode = prevNode;
        resPathEdge.push_back(prevEdge);

        edgeMapList[task.id].insert({prevEdge, 1});

        resDis.push_back(G.edgeSet[prevEdge].first.distance);
        resHop.push_back(G.edgeSet[prevEdge].first.hop);
        minChannel = std::min(minChannel, G.edgeSet[prevEdge].first.cntChannel);
    }

    std::reverse(begin(resPathEdge), end(resPathEdge));
    std::reverse(begin(resPathNode), end(resPathNode));
    std::reverse(begin(resDis), end(resDis));
    std::reverse(begin(resHop), end(resHop));

    task.totalMinChannel = minChannel;
    task.totalSumDistance = toDis;
    task.pathEdge = std::move(resPathEdge);
    task.pathNode = std::move(resPathNode);
    task.dis = std::move(resDis);
    task.hop = std::move(resHop);
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

    int nowDis = 0, nowHop = 0;

    for (int i = 0; i < (int)task.pathNode.size() - 1; i++) {
        if (nowDis + task.dis[i] > D || nowHop + task.hop[i] > H) {
            task.station.push_back(task.pathNode[i]);
            nowDis = 0;
            nowHop = 0;
        }
        nowDis += task.dis[i];
        nowHop += task.hop[i];
    }
}

void solveSingleTask(Graph &G, Task &task, std::vector<int> &channelList,
                     std::vector<std::unordered_map<int, int>> &edgeMapList) {
    i64 toDis = INF;
    int bestChannel = -1;

    if (task.r == 2 && channelList[task.id] != -1) {
        bestChannel = channelList[task.id];
    } else {
        std::vector<int> pList;
        for (int i = 0; i < P; i++) {
            pList.push_back(i);
        }

        std::random_shuffle(begin(pList), end(pList));
        pList.resize(std::min(15, P));
        for (auto p : pList) {
            i64 nowToDis = INF;
            auto curLast = singleChannelDijkstra(
                G, task.from, task.to, p, task.id, nowToDis, edgeMapList);
            if (nowToDis < toDis) {
                toDis = nowToDis;
                bestChannel = p;
            }
        }
    }

    auto curLast = singleChannelDijkstra(G, task.from, task.to, bestChannel,
                                         task.id, toDis, edgeMapList);

    completeTask(G, task, curLast, toDis, bestChannel, edgeMapList);
    channelList[task.id] = bestChannel;
}

void solveAllTask(Graph &G, std::vector<Task> &taskList,
                  std::vector<std::unordered_map<int, int>> &edgeMapList) {
    std::vector<int> channelList(T, -1);
    for (int i = 0; i < (int)taskList.size(); i++) {
        solveSingleTask(G, taskList[i], channelList, edgeMapList);
    }
}

i64 solveScoreTask(const Graph &G, const std::vector<Task> &taskList) {
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
    for (const auto &p : taskList) {
        ret += 1LL * p.pathEdge.size() * 1;
        ret += 100LL * p.station.size();
    }
    return ret;
}

void outputAnswer(Graph &G, std::vector<Task> &taskList) {
    std::vector<std::tuple<int, int, int>> edgeList;
    memset(cnt, 0, sizeof cnt);

    for (int from = 0; from < N; ++from) {
        for (auto e : G.adj[from]) {
            int to = e.to;
            if (e.id < M) continue;
            if (from > to) continue;
            if (e.deleted) {
                ++cnt[e.id];
                continue;
            }
            if (e.id >= M && e.cntChannel >= P) {
                ++cnt[e.id];
                continue;
            }
            edgeList.emplace_back(e.id, from, to);
        }
    }
    for (int i = 1; i < 100000; ++i) {
        cnt[i] += cnt[i - 1];
    }

    std::cout << edgeList.size() << "\n";
    std::sort(begin(edgeList), end(edgeList));
    for (auto [id, u, v] : edgeList) {
        auto [copyId, copyD, copyH] = G.mat[u][v];
        std::cout << copyId << "\n";
    }

    std::vector<std::vector<int>> taskIndex(T);
    for (int i = 0; i < (int)taskList.size(); ++i) {
        taskIndex[taskList[i].id].push_back(i);
    }
    for (int pt = 0; pt < T; ++pt) {
        for (auto i : taskIndex[pt]) {
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
    }

    return;
}

void init(Graph &G, std::vector<Task> &taskList, Graph &lastG,
          std::vector<std::unordered_map<int, int>> &edgeMapList, int t = 0) {
    if (!t) {
        memset(cnt, 0, sizeof cnt);
        for (auto task : taskList) {
            i64 toDis = INF;
            auto curLast = singleChannelDijkstra(G, task.from, task.to, task.id,
                                                 0, toDis, edgeMapList, true);
            int curNode = task.to;

            while (curNode != -1) {
                auto [prevNode, prevEdge] = curLast[curNode];
                if (prevEdge != -1) {
                    ++cnt[prevEdge];
                }
                curNode = prevNode;
            }
        }

        for (int edgeId = 0; edgeId < M; ++edgeId) {
            if (cnt[edgeId] >= P) {
                int u = G.edgeSet[edgeId].first.from;
                int v = G.edgeSet[edgeId].first.to;
                if (u > v) std::swap(u, v);

                auto [id, dis, hop] = G.mat[u][v];
                G.addEdge(u, v, dis, hop);
                int nowEdgeId = G.cnt - 1;
                G.edgeSet[nowEdgeId].first.preAdd = true;
                int uId = G.edgeSet[nowEdgeId].second.first;
                int vId = G.edgeSet[nowEdgeId].second.second;

                G.adj[u][uId].preAdd = true;
                G.adj[v][vId].preAdd = true;

                G.edgeSet[edgeId].first.added = true;
            }
        }
        return;
    }

    std::vector<Edge> edgeList;

    for (int i = M; i < lastG.cnt; ++i) {
        if (lastG.edgeSet[i].first.deleted) continue;
        if (lastG.edgeSet[i].first.cntChannel >= P) continue;
        edgeList.emplace_back(lastG.edgeSet[i].first);
    }

    std::sort(edgeList.begin(), edgeList.end(),
              [&](auto A, auto B) { return A.cntChannel < B.cntChannel; });

    double rate = 0.72;
    for (int i = 0; i < (int)((edgeList.size()) * rate); ++i) {
        int u = edgeList[i].from;
        int v = edgeList[i].to;
        if (u > v) std::swap(u, v);

        auto [id, dis, hop] = G.mat[u][v];
        G.addEdge(u, v, dis, hop);
        int nowEdgeId = G.cnt - 1;
        G.edgeSet[nowEdgeId].first.preAdd = true;
        int uId = G.edgeSet[nowEdgeId].second.first;
        int vId = G.edgeSet[nowEdgeId].second.second;

        G.adj[u][uId].preAdd = true;
        G.adj[v][vId].preAdd = true;
    }
}

//=================== sort by weight

namespace WeightSort {

std::vector<std::pair<int, int>> dijkstra(const Graph &G,
                                          const std::vector<int> &edgeWeight,
                                          int from, int to, i64 &toDis) {
    std::vector<int> vis(N, 0);
    std::vector<std::pair<int, int>> last(N, {-1, -1});
    std::vector<i64> dis(N, INF);

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

        for (const auto &edge : G.adj[u]) {
            int v = edge.to;
            if (vis[v]) continue;
            int nxtDis = dis[u] + edgeWeight[edge.id];
            if (nxtDis < dis[v]) {
                dis[v] = nxtDis;
                last[v] = {u, edge.id};
                Q.push(std::make_pair(-dis[v], v));
            }
        }
    }

    toDis = dis[to];
    if (vis[to] == 0) return {};
    return last;
}

void sortByWeight(const Graph &G, std::vector<Task> &taskList) {
    std::vector<int> edgeWeight(G.cnt, 1);
    std::vector<int> edgeCnt(G.cnt, 0);

    for (auto &task : taskList) {
        i64 toDis;
        auto curLast = dijkstra(G, edgeWeight, task.from, task.to, toDis);
        if (!curLast.size()) {
            continue;
        }
        int curNode = task.to;

        task.addEdge = 0;
        while (curNode != -1) {
            auto [prevNode, prevEdge] = curLast[curNode];
            if (prevEdge == -1) break;
            edgeCnt[prevNode] = edgeCnt[prevNode] + 1;
            curNode = prevNode;
        }
    }

    for (int i = 0; i < G.cnt; ++i) {
        edgeWeight[i] = (edgeCnt[i]) * 2;
    }

    for (auto &task : taskList) {
        i64 toDis = 0;
        auto curLast = dijkstra(G, edgeWeight, task.from, task.to, toDis);
        task.shortestPathLen = toDis;
    }
}

}  // namespace WeightSort

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    Timer timer;
    timer.reset();
    std::cin >> N >> M >> T >> R >> P >> D >> H;

    Graph G(N, P);
    for (int i = 0; i < M; i++) {
        int _u, _v, _d, _h;
        std::cin >> _u >> _v >> _d >> _h;
        G.addEdge(_u, _v, _d, _h);
    }

    std::vector<Task> taskList;
    std::vector<Task> realTaskList;
    for (int i = 0; i < T; i++) {
        int _from, _to, _r;
        std::cin >> _from >> _to >> _r;
        auto _dis = solveTaskDistance(G, _from, _to);
        realTaskList.push_back({i, _from, _to, _dis, _r});
    }

    for (auto task : realTaskList) {
        for (int i = 0; i < task.r; ++i) {
            taskList.push_back(task);
        }
    }

    if (M * T * P >= 1E9) {
        WeightSort::sortByWeight(G, taskList);
    }

    int timeIWantToTry = 3;
    if (N <= 500)
        timeIWantToTry = 50;
    else if (N <= 1000)
        timeIWantToTry = 10;
    else if (N <= 2000)
        timeIWantToTry = 10;
    else
        timeIWantToTry = 3;

    if (N <= 100) timeIWantToTry = 20;
    std::vector<Graph> GList(timeIWantToTry, G);
    std::vector<std::vector<Task>> listTaskList(timeIWantToTry, taskList);
    std::vector<unsigned int> seedList;
    for (int i = 1; i <= timeIWantToTry; i++) {
        seedList.push_back(i);
    }

    int bestSeed = 0;
    i64 bestScore = -1;

    for (int i = 0; i < timeIWantToTry; i++) {
        std::vector<std::unordered_map<int, int>> edgeMapList(T);
        std::srand(seedList[i]);
        std::random_shuffle(begin(listTaskList[i]), end(listTaskList[i]));
        std::random_shuffle(begin(listTaskList[i]), end(listTaskList[i]));
        std::sort(begin(listTaskList[i]), end(listTaskList[i]),
                  [&](auto cmpA, auto cmpB) {
                      int weightA = 0, weightB = 0;
                      if (cmpA.r == 2) weightA += 10;
                      if (cmpB.r == 2) weightB += 10;
                      return cmpA.shortestPathLen + weightA >
                             cmpB.shortestPathLen + weightB;
                  });

        if (i) {
            std::sort(begin(listTaskList[i - 1]), end(listTaskList[i - 1]),
                      [&](auto A, auto B) {
                          if (A.addEdge > B.addEdge) return true;
                          if (A.addEdge < B.addEdge) return false;
                          return A.shortestPathLen > B.shortestPathLen;
                      });

            for (int j = 0; j < (int)listTaskList[i].size(); ++j) {
                listTaskList[i][j].dis = listTaskList[i - 1][j].dis;
                listTaskList[i][j].from = listTaskList[i - 1][j].from;
                listTaskList[i][j].to = listTaskList[i - 1][j].to;
                listTaskList[i][j].id = listTaskList[i - 1][j].id;
                listTaskList[i][j].r = listTaskList[i - 1][j].r;
            }
        }

        if (N <= 1000) {
            if ((i % 5) == 1) {
                init(GList[i], listTaskList[i], GList[bestSeed], edgeMapList,
                     1);
            }
        }
        if (M * T * P >= 1E9) {
            std::sort(begin(listTaskList[i]), end(listTaskList[i]),
                      [&](auto cmpA, auto cmpB) {
                          int weightA = 0, weightB = 0;
                          if (cmpA.r == 2) weightA += 8;
                          if (cmpB.r == 2) weightB += 8;
                          return cmpA.shortestPathLen + weightA >
                                 cmpB.shortestPathLen + weightB;
                      });
        }

        solveAllTask(GList[i], listTaskList[i], edgeMapList);

        auto curScore = solveScoreTask(GList[i], listTaskList[i]);
        if (curScore < bestScore || bestScore == -1) {
            bestScore = curScore;
            bestSeed = i;
        }
        unsigned int singleRunTime = (timer.elapsed() + i) / (i + 1);

        if (timer.elapsed() + singleRunTime >= 1000 * 115) break;
    }

    outputAnswer(GList[bestSeed], listTaskList[bestSeed]);

    return 0;
}