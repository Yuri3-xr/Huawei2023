#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

constexpr double EPS = 1E-10;
constexpr int numRobots = 4;
constexpr std::array<std::string, 5> opName = {"forward", "rotate", "but",
                                               "sell", "destroy"};

template <class T>
int sgn(T x) {
    return (x > 0) - (x < 0);
}

template <class T>
bool equalZero(T x) {
    return std::fabs(x) < EPS;
}

/*----------------------------------头文件，常数及宏定义-------------------------------*/

struct Vec {
    double x, y;
    Vec() : x(0), y(0){};
    Vec(double x, double y) : x(x), y(y){};
    bool operator<(Vec p) const { return std::tie(x, y) < std::tie(p.x, p.y); }
    bool operator==(Vec p) const {
        return std::tie(x, y) == std::tie(p.x, p.y);
    }
    Vec operator+(Vec p) const { return Vec(x + p.x, y + p.y); }
    Vec operator-(Vec p) const { return Vec(x - p.x, y - p.y); }
    Vec operator*(double d) const { return Vec(x * d, y * d); }
    Vec operator/(double d) const { return Vec(x / d, y / d); }
    double dot(Vec p) const { return x * p.x + y * p.y; }
    double cross(Vec p) const { return x * p.y - y * p.x; }
    double cross(Vec a, Vec b) const { return (a - *this).cross(b - *this); }
    double dist2() const { return x * x + y * y; }
    double dist() const { return sqrt((double)dist2()); }
    // angle to x-axis in interval [-pi, pi]
    double angle() const { return atan2(y, x); }
    Vec unit() const { return *this / dist(); }  // makes dist()=1
    Vec perp() const { return Vec(-y, x); }      // rotates +90 degrees
    Vec perpInv() const { return Vec(y, -x); }   // rotates -90 degrees

    Vec normal() const { return perp().unit(); }
    // returns point rotated 'a' radians ccw around the origin
    Vec rotate(double a) const {
        return Vec(x * cos(a) - y * sin(a), x * sin(a) + y * cos(a));
    }
    friend std::ostream& operator<<(std::ostream& os, Vec p) {
        return os << "(" << p.x << "," << p.y << ")";
    }
};

double dist(Vec A, Vec B) {
    auto rx = A.x - B.x;
    auto ry = A.y - B.y;
    return sqrt(rx * rx + ry * ry);
};

double dist2(Vec A, Vec B) {
    auto rx = A.x - B.x;
    auto ry = A.y - B.y;
    return (rx * rx + ry * ry);
}

/*----------------------------------坐标类-------------------------------*/

struct Robots {
    Vec cd;            // 坐标
    double angV;       // 角速度（弧度制），
    Vec lineV;         // 线速度（向量）
    double direction;  // 方向（弧度制）
    int goods;  // 当前的携带的货物，0表示当前没有携带货物
    int currentWork;      // 当前前往的工作台，-1表示当前闲置
    int currentWorkshop;  // 当前所处的工作台
    double timeCoeff;     // 时间价值系数
    double crashCoeff;    // 碰撞价值系数
    Robots(){};
    Robots(int currentWorkshop, double x, double y, double angV, double lineVx,
           double lineVy, double direction, int goods, int currentWork,
           double timeCoeff, double crashCoeff)
        : currentWorkshop(currentWorkshop),
          cd(x, y),
          angV(angV),
          lineV(lineVx, lineVy),
          direction(direction),
          goods(goods),
          timeCoeff(timeCoeff),
          crashCoeff(crashCoeff){};

    std::pair<int, std::pair<Vec, double>> getTrack() {
        /*
            给定当前坐标以及朝向和角速度，线速度，判断机器人轨迹
            返回一个pair
            first：表示直线（1）/圆弧（2）
            second：
                （1）情况下，返回当前坐标和线速度的模长
                （2）情况下，返回圆心坐标和半径
        */
        if (equalZero(angV)) return {1, {cd, lineV.dist()}};
        double r = std::abs(lineV.dist() / (2 * angV));

        Vec ret = (sgn(angV) > 0) ? lineV.perpInv() : lineV.perp();
        ret = ret + cd;
        return {2, {ret, r}};
    };

    friend std::ostream& operator<<(std::ostream& os, Robots x) {
        os << std::fixed << std::setprecision(6);
        os << "当前坐标为：" << x.cd << std::endl;
        os << "当前角速度为：" << x.angV << std::endl;
        os << "当前线速度为：" << x.lineV << std::endl;
        os << "当前朝向为：" << x.direction << std::endl;
        os << "当前携带的货物为：" << x.goods << std::endl;
        os << "当前前往的工作台为：" << x.currentWork << std::endl;
        return os;
    }
};

/*----------------------------------机器人类-------------------------------*/

struct Workshops {
    Vec cd;
    int remTime;   // 剩余生产时间
    int matState;  // 原材料状态
    int repState;  // 产品格状态
    int type;      // 工作台种类
    Workshops(){};
    Workshops(double x, double y, int remTime, int matState, int repStata,
              int type)
        : cd(x, y),
          remTime(remTime),
          matState(matState),
          repState(repState),
          type(type){};
};

/*----------------------------------工作台类-------------------------------*/

struct Operate {
    int type;
    std::vector<std::variant<int, double>> parameter;
    Operate() = delete;
    Operate(int type, std::vector<std::variant<int, double>> parameter)
        : type(type), parameter(parameter) {
        assert(0 <= type && type <= 4);
        if (0 <= type && type <= 1) assert(parameter.size() == 1);
        if (2 <= type && type <= 4) assert(parameter.size() == 2);
    };
    friend std::ostream& operator<<(std::ostream& os, Operate p) {
        os << opName[p.type] << " ";
        if (p.parameter.size() == 1)
            os << std::get<int>(p.parameter[0]);
        else
            os << std::get<int>(p.parameter[0]) << " "
               << std::get<double>(p.parameter[1]);
        os << std::endl;
        return os;
    }
};

/*----------------------------------操作类-------------------------------*/

int frameID = 0;
int64_t money = 200'000;
int numWorkshops;

std::array<Robots, numRobots> robot;  // 4个机器人
std::vector<Workshops> workshop;      // 工作台

/*----------------------------------全局数据结构-------------------------------*/

bool readUntilOK() {
    std::string line;
    while (getline(std::cin, line)) {
        if (line == "ok") return true;
    }

    return false;
}

void updateWorkshopState() {
    std::cin >> money;
    std::cin >> numWorkshops;

    workshop.resize(numWorkshops);

    for (int i = 0; i < numWorkshops; i++) {
        double x, y;
        int remTime, matState, repStata, type;
        std::cin >> type >> x >> y >> remTime >> matState >> repStata;
        workshop[i] = Workshops(x, y, remTime, matState, repStata, type);
    }
}

void updateRobotState() {
    for (int i = 0; i < numRobots; i++) {
        int currentWorkshop, goods;
        double timeCoeff, crashCoeff;
        double angV;
        Vec lineV;
        double direction;
        double x, y;

        std::cin >> currentWorkshop >> goods >> timeCoeff >> crashCoeff >>
            angV >> lineV.x >> lineV.y >> direction >> x >> y;

        robot[i] = Robots(currentWorkshop, x, y, angV, lineV.x, lineV.y,
                          direction, goods, -1, timeCoeff, crashCoeff);
    }
}

std::vector<Operate> distributePlan() {
    /*
        初步想法是贪心地搞，具体实现见naive_greedy分支
    */
}
/*----------------------------------操作函数-------------------------------*/

int main() {
    readUntilOK();
    std::cout << "OK" << std::endl;

    while (std::cin >> frameID) {
        updateWorkshopState();
        updateRobotState();

        readUntilOK();
        std::cout << frameID << std::endl;

        auto plan = distributePlan();
        for (const auto& out : plan) std::cout << out;

        std::cout << "OK" << std::endl;
        std::cout.flush();
    }
    return 0;
}
