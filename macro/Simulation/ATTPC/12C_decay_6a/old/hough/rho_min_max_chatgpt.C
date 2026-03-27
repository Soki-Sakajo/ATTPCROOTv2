#include "TMinuit.h"
#include "TMath.h"
#include <iostream>

// グローバル変数（TMinuit用）
double f_to_minimize(const double *vars) {
    double a = vars[0]; // psi
    double b = vars[1]; // theta
    double c = vars[2]; // phi
    double x = vars[3];
    double y = vars[4];
    double z = vars[5];

    double A = TMath::Cos(a) * TMath::Cos(c) - TMath::Sin(a) * TMath::Cos(b) * TMath::Sin(c);
    double B = TMath::Sin(a) * TMath::Cos(c) + TMath::Cos(a) * TMath::Cos(b) * TMath::Sin(c);
    double C = TMath::Sin(b) * TMath::Sin(c);

    return x * A + y * B + z * C;
}

// Minuit用の関数
void minuit_function(int &npar, double *gin, double &fval, double *vars, int iflag) {
    fval = f_to_minimize(vars);
}

// 最適化を行う関数
void rho_min_max_chatgpt() {
  //void find_extremum() {
    TMinuit minuit(6); // 6変数
    minuit.SetFCN(minuit_function);

    // パラメータの初期値と範囲を設定
    double initial[6] = {TMath::Pi(), TMath::Pi() / 4, TMath::Pi(), 0, 0, 500};
    double step[6] = {0.1, 0.1, 0.1, 1.0, 1.0, 10.0};
    double min[6] = {0, 0, 0, -280, -280, 0};
    double max[6] = {2 * TMath::Pi(), TMath::Pi() / 2, 2 * TMath::Pi(), 280, 280, 1000};

    for (int i = 0; i < 6; ++i) {
        minuit.DefineParameter(i, Form("var%d", i + 1), initial[i], step[i], min[i], max[i]);
    }

    // 最大値探索
    minuit.Migrad();
    double max_val, max_err;
    double max_coords[6];
    for (int i = 0; i < 6; ++i) {
        minuit.GetParameter(i, max_coords[i], max_err);
    }
    max_val = f_to_minimize(max_coords);

    // 出力
    std::cout << "最大値: " << max_val << "\n位置: ";
    for (int i = 0; i < 6; ++i) {
        std::cout << max_coords[i] << " ";
    }
    std::cout << std::endl;

    // 最小値探索
    minuit.SetFCN([](int &npar, double *gin, double &fval, double *vars, int iflag) {
        fval = -f_to_minimize(vars); // 符号を反転
    });
    minuit.Migrad();
    double min_val;
    double min_coords[6];
    for (int i = 0; i < 6; ++i) {
        minuit.GetParameter(i, min_coords[i], max_err);
    }
    min_val = -f_to_minimize(min_coords); // 符号を元に戻す

    // 出力
    std::cout << "最小値: " << min_val << "\n位置: ";
    for (int i = 0; i < 6; ++i) {
        std::cout << min_coords[i] << " ";
    }
    std::cout << std::endl;
}
