#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <tuple>
#include <utility>
#include <vector>
using namespace std;

constexpr double EPS = 1E-10;

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
    int goods;  // 当前的携带的货物，-1表示当前没有携带货物
    int currentWork;  // 当前前往的工作台，-1表示当前闲置
    Robots() = delete;
    Robots(double x, double y, double angV, double lineVx, double lineVy,
           double direction, int goods = -1, int currentWork = -1)
        : cd(x, y),
          angV(angV),
          lineV(lineVx, lineVy),
          direction(direction),
          goods(goods){};

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

struct Workshop {
    Vec cd;
    int remTime;   // 剩余生产时间
    int matState;  // 原材料状态
    int repState;  // 产品格状态
    Workshop() = delete;
    Workshop(double x, double y, int remTime, int matState, int repStata)
        : cd(x, y), remTime(remTime), matState(matState), repState(repState){};
};

/*----------------------------------工作台类-------------------------------*/

bool readUntilOK() {
    char line[1024];
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] == 'O' && line[1] == 'K') {
            return true;
        }
        // do something
    }
    return false;
}

int main() {
    readUntilOK();
    puts("OK");
    fflush(stdout);
    int frameID;
    while (scanf("%d", &frameID) != EOF) {
        readUntilOK();
        printf("%d\n", frameID);
        int lineSpeed = 3;
        double angleSpeed = 1.5;
        for (int robotId = 0; robotId < 4; robotId++) {
            printf("forward %d %d\n", robotId, lineSpeed);
            printf("rotate %d %f\n", robotId, angleSpeed);
        }
        printf("OK\n", frameID);
        fflush(stdout);
    }
    return 0;
}
