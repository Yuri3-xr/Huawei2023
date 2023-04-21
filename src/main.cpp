#include <bits/stdc++.h>

using i64 = std::int64_t;

constexpr int INF = std::numeric_limits<int>::max() / 2.0;
constexpr uint32_t SEED = 229;

constexpr int EDGE_WEIGHT = 10;
constexpr int ADDED_EDGE_WEIGHT = 1e6;
constexpr int PRE_ADDED_EDGE_WEIGHT = 10;

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
    bool deleted = false;
    bool preAdd = false;
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
    int addEdge = 0;
};

/*
N: 点数
M: 边数
T: 任务数
P: 信道上限
D: 衰减上限
*/
int N, M, T, P, D;
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

inline int cost(int x) {
    if (M > N * 30) {
        return exp(1.0 * x / 20);
    } else {
        return x;
    }
}

std::vector<std::pair<int, int>> singleChannelDijkstra(const Graph &G, int from,
                                                       int to, int p,
                                                       int &toDis,
                                                       bool fixFlag = false) {
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

        for (const auto &edge : G.adj[u]) {
            int v = edge.to;
            if (vis[v]) continue;
            int nxtDis = dis[u], lastEdgeId = -1;
            if (edge.deleted) continue;  // has been removed.
            if (fixFlag) {
                if (edge.markChannel[p] != -1)
                    continue;  // this channel has been used.
            }
            if (edge.markChannel[p] == -1) {
                nxtDis += cost(P - edge.cntChannel);
                if (edge.preAdd) {
                    nxtDis += EDGE_WEIGHT;
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
                Q.push(std::make_pair(-dis[v], v));
            }
        }
    }

    toDis = dis[to];
    if (vis[to] == 0) return {};
    return last;
}

void completeTask(Graph &G, Task &task,
                  std::vector<std::pair<int, int>> &curLast, int toDis,
                  int bestChannel) {
    int curNode = task.to;
    std::vector<int> resPathEdge, resPathNode, resDis;
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

void solveSingleTask(Graph &G, Task &task) {
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

    auto curLast =
        singleChannelDijkstra(G, task.from, task.to, bestChannel, toDis);

    completeTask(G, task, curLast, toDis, bestChannel);
}

//=========== fix begin

void deleteEdgeToFix(Graph &G, std::vector<Task> &taskList,
                     std::vector<int> &taskIndex, int edgeId) {
    std::vector<int> taskToFix;

    // for (auto edgeId : edgeList) {
    G.edgeSet[edgeId].first.deleted = true;
    int u = G.edgeSet[edgeId].first.from;
    int v = G.edgeSet[edgeId].first.to;
    if (u > v) std::swap(u, v);

    int uId = G.edgeSet[edgeId].second.first;
    int vId = G.edgeSet[edgeId].second.second;

    G.adj[u][uId].deleted = true;
    G.adj[v][vId].deleted = true;

    auto &edge = G.edgeSet[edgeId].first;
    for (int p = 0; p < P; ++p) {
        if (edge.markChannel[p] != -1) {
            taskToFix.emplace_back(taskIndex[edge.markChannel[p]]);
            // std::cerr << edge.markChannel[p] << std::endl;
        }
    }
    // }

    taskToFix.resize(std::unique(taskToFix.begin(), taskToFix.end()) -
                     taskToFix.begin());
    // for (auto taskId : taskToFix) {
    //     std::cerr << taskId << " ";
    // }
    // std::cerr << std::endl;

    bool flag = true;
    std::vector<int> channelUsed(P, 0);
    for (auto taskId : taskToFix) {
        auto task = taskList[taskId];
        int toDis = INF, bestChannel = -1;
        for (int p = 0; p < P; ++p) {
            if (channelUsed[p]) continue;
            int nowToDis = INF;
            auto curLast =
                singleChannelDijkstra(G, task.from, task.to, p, nowToDis, true);
            // std::cerr << nowToDis << std::endl;
            if ((curLast.size()) && (nowToDis < toDis)) {
                toDis = nowToDis;
                bestChannel = p;
            }
        }

        if (bestChannel == -1) {
            // fix failed.
            flag = false;
            break;
        }

        channelUsed[bestChannel] = taskId;
    }

    memset(cnt, 0, sizeof cnt);
    if (flag) {
        // std::cerr << "fix: ok" << std::endl;
        for (int p = 0; p < P; ++p) {
            if (channelUsed[p]) {
                auto &task = taskList[channelUsed[p]];
                int toDis = INF;
                auto curLast = singleChannelDijkstra(G, task.from, task.to, p,
                                                     toDis, true);

                int curNode = task.to;

                task.station.clear();
                completeTask(G, task, curLast, toDis, p);
            }
        }
        // int _cnt = 0;
        // for (int i = 0; i < G.cnt; ++i) _cnt += (cnt[i] > 0);
        // std::cerr << _cnt << std::endl;
        // if (_cnt < (int)(edgeList.size())) {
        //     for (int p = 0; p < P; ++p) {
        //         if (channelUsed[p]) {
        //             auto &task = taskList[channelUsed[p]];
        //             int toDis = INF;
        //             auto curLast = singleChannelDijkstra(G, task.from,
        //             task.to, p, toDis, true); int curNode = task.to;
        //             completeTask(G, task, curLast, toDis, p);
        //         }
        //     }
        // }
        return;
    }

    // for (auto edgeId : edgeList) {
    G.edgeSet[edgeId].first.deleted = false;

    // int u = G.edgeSet[edgeId].first.from;
    // int v = G.edgeSet[edgeId].first.to;
    // if (u > v) std::swap(u, v);

    // int uId = G.edgeSet[edgeId].second.first;
    // int vId = G.edgeSet[edgeId].second.second;
    G.adj[u][uId].deleted = false;
    G.adj[v][vId].deleted = false;
    // }
}

//=========== fix end

void solveAllTask(Graph &G, std::vector<Task> &taskList) {
    for (int i = 0; i < T; i++) {
        solveSingleTask(G, taskList[i]);
    }

    // fix init.
    std::vector<int> taskIndex(T);
    for (int i = 0; i < T; ++i) {
        taskIndex[taskList[i].id] = i;
    }

    std::vector<int> edgeList;

    int H = 1;
    for (int i = M; i < G.cnt; ++i) {
        if (P - G.edgeSet[i].first.cntChannel <= H) {
            // edgeList.emplace_back(i);
            deleteEdgeToFix(G, taskList, taskIndex, i);
        }
    }
    // std::cerr << edgeList.size() << std::endl;
    // if (edgeList.size()) {
    //     std::cerr << "/??" << std::endl;
    //     deleteEdgeToFix(G, taskList, taskIndex, edgeList);
    //     std::cerr << "/##" << std::endl;
    // }
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

void outputAnswer(const Graph &G, std::vector<Task> &taskList) {
    // std::cout << G.cnt - M << "\n";

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

void init(Graph &G, std::vector<Task> &taskList, Graph &lastG, int t = 0) {
    if (!t) {
        memset(cnt, 0, sizeof cnt);
        for (auto task : taskList) {
            int toDis = INF;
            auto curLast =
                singleChannelDijkstra(G, task.from, task.to, 0, toDis, true);
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
                G.edgeSet[edgeId].first.deleted = true;
                int u = G.edgeSet[edgeId].first.from;
                int v = G.edgeSet[edgeId].first.to;
                if (u > v) std::swap(u, v);

                G.addEdge(u, v, G.mat[u][v]);
                int nowEdgeId = G.cnt - 1;
                G.edgeSet[nowEdgeId].first.preAdd = true;
                int uId = G.edgeSet[nowEdgeId].second.first;
                int vId = G.edgeSet[nowEdgeId].second.second;

                G.adj[u][uId].preAdd = true;
                G.adj[v][vId].preAdd = true;
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

    // std::cerr << edgeList.size() << std::endl;
    double rate = 0.5;
    for (int i = 0; i < (int)((edgeList.size()) * rate); ++i) {
        int u = edgeList[i].from;
        int v = edgeList[i].to;
        if (u > v) std::swap(u, v);

        G.addEdge(u, v, G.mat[u][v]);
        int nowEdgeId = G.cnt - 1;
        G.edgeSet[nowEdgeId].first.preAdd = true;
        int uId = G.edgeSet[nowEdgeId].second.first;
        int vId = G.edgeSet[nowEdgeId].second.second;

        G.adj[u][uId].preAdd = true;
        G.adj[v][vId].preAdd = true;
    }
}

// //=================== tarjan

constexpr int MAX_N = 5e3 + 7;
constexpr int MAX_M = 5e3 + 7;
int dfn[MAX_N], low[MAX_N], bridge[MAX_M], belong[MAX_N];
std::stack<int> st;
int order, num;

void tarjan(const Graph &G, int u, int fa) {
    dfn[u] = low[u] = ++order;
    bool flag = true;
    st.push(u);
    // for (int i = head[u]; i; i = e[i].next) {
    for (auto edge : G.adj[u]) {
        int v = edge.to;
        if (v == fa && flag) {
            flag = false;
            continue;
        }
        if (!dfn[v]) {
            tarjan(G, v, u);
            low[u] = std::min(low[u], low[v]);
            if (low[v] > dfn[u]) {
                bridge[edge.id] = 1;
            }
        } else if (dfn[v] < dfn[u]) {
            low[u] = std::min(low[u], dfn[v]);
        }
    }
    if (dfn[u] == low[u]) {
        num++;
        int tmp;
        do {
            tmp = st.top();
            belong[tmp] = num;
            st.pop();
        } while (tmp != u);
    }
}

std::vector<int> E[MAX_N];
int siz[MAX_N], dep[MAX_N], son[MAX_N], top[MAX_N], fa[MAX_N];
int fat[MAX_N];

int get_fa(int x) { return x == fat[x] ? x : get_fa(fat[x]); }

void dfs1(int u) {
    siz[u] = 1;
    // for (int i = head[u]; i; i = e[i].next) {
    for (auto v : E[u]) {
        if (v == fa[u]) continue;
        fa[v] = u;
        dep[v] = dep[u] + 1;
        dfs1(v);
        siz[u] += siz[v];
        if (siz[v] > siz[son[u]]) son[u] = v;
    }
}

void dfs2(int u, int tp) {
    top[u] = tp;
    if (son[u]) dfs2(son[u], tp);
    for (auto v : E[u]) {
        if (v == fa[u] || v == son[u]) continue;
        dfs2(v, v);
    }
}

int lca(int u, int v) {
    while (top[u] != top[v]) {
        if (dep[top[u]] > dep[top[v]]) std::swap(u, v);
        v = fa[top[v]];
    }
    return dep[u] < dep[v] ? u : v;
}

std::vector<int> bccWeight(const Graph &G, std::vector<Task> &taskList) {
    for (int i = 0; i < N; ++i) {
        if (!dfn[i]) {
            tarjan(G, i, -1);
        }
    }

    for (int i = 1; i <= num; ++i) fat[i] = i;
    for (int u = 0; u < N; ++u) {
        for (auto edge : G.adj[u]) {
            int v = edge.to;
            int uu = belong[u], vv = belong[v];
            int fu = get_fa(uu), fv = get_fa(vv);
            if (fu == fv) continue;
            E[uu].emplace_back(vv);
            E[vv].emplace_back(uu);
            fat[uu] = vv;
        }
    }

    dfs1(1);
    dfs2(1, 1);

    std::vector<int> dis;
    for (auto task : taskList) {
        int u = task.from, v = task.to;
        int fu = belong[u], fv = belong[v];
        int _lca = lca(fu, fv);
        dis.emplace_back(dep[fu] + dep[fv] - dep[_lca] * 2);
    }

    return dis;
}

//===================

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

    int timeIWantToTry = 30;
    if (N <= 500)
        timeIWantToTry = 100;
    else if (N <= 1000)
        timeIWantToTry = 100;
    else if (N <= 2000)
        timeIWantToTry = 100;
    else
        timeIWantToTry = 30;

    if (N <= 100) timeIWantToTry = 200;
    std::vector<Graph> GList(timeIWantToTry, G);
    std::vector<std::vector<Task>> listTaskList(timeIWantToTry, taskList);
    std::vector<unsigned int> seedList;
    for (int i = 1; i <= timeIWantToTry; i++) {
        seedList.push_back(i);
    }

    // std::vector<int> bccDis = bccWeight(G, taskList);

    int bestSeed = 0, bestScore = INF;

    for (int i = 0; i < timeIWantToTry; i++) {
        std::srand(seedList[i]);
        std::random_shuffle(begin(listTaskList[i]), end(listTaskList[i]));
        std::random_shuffle(begin(listTaskList[i]), end(listTaskList[i]));
        for (auto &p : listTaskList[i]) {
            p.shortestPathLen += rand() % 4;
        }

        std::sort(begin(listTaskList[i]), end(listTaskList[i]),
                  [&](auto cmpA, auto cmpB) {
                      // if (bccDis[cmpA.id] > bccDis[cmpB.id]) return true;
                      // if (bccDis[cmpA.id] < bccDis[cmpB.id]) return false;
                      return cmpA.shortestPathLen > cmpB.shortestPathLen;
                  });

        if (i) {
            std::sort(begin(listTaskList[i - 1]), end(listTaskList[i - 1]),
                      [&](auto A, auto B) {
                          if (A.addEdge > B.addEdge) return true;
                          if (A.addEdge < B.addEdge) return false;
                          return A.shortestPathLen > B.shortestPathLen;
                      });

            for (int j = 0; j < T; ++j) {
                listTaskList[i][j].dis = listTaskList[i - 1][j].dis;
                listTaskList[i][j].from = listTaskList[i - 1][j].from;
                listTaskList[i][j].to = listTaskList[i - 1][j].to;
                listTaskList[i][j].id = listTaskList[i - 1][j].id;
            }
        }

        if (N <= 1000) {
            if ((i % 5) == 1) {
                init(GList[i], listTaskList[i], GList[bestSeed]);
            }
        } else if (1.0 * (M - N) / N <= 0.25) {
            init(GList[i], listTaskList[i], GList[bestSeed]);
        }
        // init(GList[i], listTaskList[i], GList[bestSeed], 0);
        solveAllTask(GList[i], listTaskList[i]);

        auto curScore = solveScoreTask(GList[i], listTaskList[i]);
        // std::cerr << i << " " << curScore << std::endl;
        if (curScore < bestScore) {
            bestScore = curScore;
            bestSeed = i;
        }
        unsigned int singleRunTime = (timer.elapsed() + i) / (i + 1);

        if (timer.elapsed() + singleRunTime >= 1000 * 110) break;
    }

    // std::cerr << bestSeed << std::endl;

    // std::cerr << bestScore << std::endl;
    outputAnswer(GList[bestSeed], listTaskList[bestSeed]);

    // std::cerr << timer.elapsed() << "ms" << std::endl;

    return 0;
}