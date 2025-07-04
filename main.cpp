#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <iomanip>

// ==== Декларації ====
double fun(double x, double y, double z);
double Grs(double u, double v, double w);
double Rrz(double p, double q, double r);
double Alg2_Rrz(double p, double q, double r);
double Alg3_Rrz(double p, double q, double r);
double Alg4(double x, double y, double z);
double T_interp(double x);
double U_interp(double x);

// ==== main ====
int main() {
    double x, y, z;
    std::cout << "Enter x, y, z: ";
    if (!(std::cin >> x >> y >> z)) {
        std::cerr << "Невірний ввід\n";
        return 1;
    }
    double res = fun(x, y, z);
    std::cout << "fun(x,y,z) = "
        << std::fixed << std::setprecision(5)
        << res << "\n";
    return 0;
}

// ==== Алг.1 ====
double fun(double x, double y, double z) {
    return x * Grs(x, y, z) + y * Grs(x, z, y);
}

double Grs(double u, double v, double w) {
    return 0.1389 * Rrz(u, v, w)
        + 1.8389 * Rrz(u - v, w, v);
}

// ==== Алг.1 п.3–6 з двома різними умовами ====
double Rrz(double p, double q, double r) {
    double c1 = r * r + p * q;
    if (c1 <= 0) {
        // 5.1 z^2+xy ≤0
        return Alg2_Rrz(p, q, r);
    }
    double c2 = p * p + r * q;
    if (c2 <= 0) {
        // 5.2 x^2+zy ≤0
        return Alg3_Rrz(p, q, r);
    }
    // Інакше — читаємо T та U; якщо файл недоступний, Alg4
    double Tu = T_interp(p);
    double Uu = U_interp(p);
    double Tv = T_interp(q);
    double Uv = U_interp(q);
    return (Tu + Uv) - (Tv + Uu);
}

// ==== Алг.2 (без рекурсії!) ====
double Alg2_Rrz(double p, double q, double r) {
    // 1) Qrz1:
    double Q = (r > p ? p * r * q : q * r * p);
    // 2) Srs1:
    double S1 = (std::fabs(p) < 1 ? Q * r * p : Q * q * r);
    // 3) Rrz = x*y*Qrz1 чи x*z*Qrz1 (п.1 Алг.2)
    if (p > q) {
        return p * q * Q;       // приклад: x*y*Qrz1(y,z)
    }
    else {
        return p * r * Q;       // приклад: x*z*Qrz1(x,y)
    }
}

// ==== Алг.3 (без рекурсії!) ====
double Alg3_Rrz(double p, double q, double r) {
    // 1) Qrz2:
    double Q = (r > p ? p * r * q : q * r * p);
    // 2) Srs2:
    double S2 = (std::fabs(p) < 1 ? Q * r * p : Q * q * r);
    // 3) Rrz = x*y*Qrz2 чи y*z*Qrz2 (п.1 Алг.3)
    if (p > q) {
        return p * q * Q;       // приклад
    }
    else {
        return q * r * Q;
    }
}

// ==== Алг.4 ====
double Alg4(double x, double y, double z) {
    return 1.3498 * x + 2.2362 * y * z - 2.348 * x * y;
}

// ==== Інтерполяція ====
struct Table {
    std::vector<double> xs, ts, us;
};

bool loadTable(const char* fname, Table& T) {
    std::ifstream f(fname);
    if (!f.is_open()) return false;
    T.xs.clear(); T.ts.clear(); T.us.clear();
    std::string line;
    while (std::getline(f, line)) {
        std::istringstream ss(line);
        double x, v1, v2;
        if (ss >> x >> v1 >> v2) {
            T.xs.push_back(x);
            T.ts.push_back(v1);
            T.us.push_back(v2);
        }
    }
    return true;
}

double interp(const Table& T, double x, const std::vector<double>& vals) {
    int n = (int)T.xs.size();
    if (n == 0) return NAN;
    if (x <= T.xs.front()) return vals.front();
    if (x >= T.xs.back())  return vals.back();
    for (int i = 0; i < n - 1; ++i) {
        if (x >= T.xs[i] && x <= T.xs[i + 1]) {
            double x0 = T.xs[i], y0 = vals[i];
            double x1 = T.xs[i + 1], y1 = vals[i + 1];
            return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
        }
    }
    return NAN;
}

double T_interp(double x) {
    Table tbl;
    const char* fn = (std::fabs(x) <= 1 ? "dat_X_1_1.dat"
        : x < -1 ? "dat_X00_1.dat"
        : "dat_X1_00.dat");
    if (!loadTable(fn, tbl)) return Alg4(x, 0, 0);
    return interp(tbl, x, tbl.ts);
}

double U_interp(double x) {
    Table tbl;
    const char* fn = (std::fabs(x) <= 1 ? "dat_X_1_1.dat"
        : x < -1 ? "dat_X00_1.dat"
        : "dat_X1_00.dat");
    if (!loadTable(fn, tbl)) return Alg4(x, 0, 0);
    return interp(tbl, x, tbl.us);
}
