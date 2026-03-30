//example 1
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH3D.h>
#include <TVector3.h>
#include <TMath.h>
#include <TRandom3.h>
#include <vector>
#include <iostream>

class PointCloud {
public:
    std::vector<TVector3> points; // 点群データ

    void addPoint(double x, double y, double z) {
        points.emplace_back(x, y, z);
    }

    TVector3 getMeanValue() const {
        TVector3 mean(0, 0, 0);
        for (const auto& p : points) {
            mean += p;
        }
        return mean * (1.0 / points.size());
    }

    void shiftToOrigin(const TVector3& shift) {
        for (auto& p : points) {
            p -= shift;
        }
    }
};

class HoughSpace {
public:
    TH3D* votingSpace; // 3Dヒストグラムとしての投票空間

    HoughSpace(double dx, double maxX, double maxY, double maxZ) {
        votingSpace = new TH3D("votingSpace", "Hough Voting Space",
                               maxX / dx, -maxX, maxX,
                               maxY / dx, -maxY, maxY,
                               maxZ / dx, -maxZ, maxZ);
    }

    ~HoughSpace() {
        delete votingSpace;
    }

    void addVote(const TVector3& point) {
        votingSpace->Fill(point.X(), point.Y(), point.Z());
    }

    void subtractVote(const TVector3& point) {
        int binX = votingSpace->GetXaxis()->FindBin(point.X());
        int binY = votingSpace->GetYaxis()->FindBin(point.Y());
        int binZ = votingSpace->GetZaxis()->FindBin(point.Z());
        double currentValue = votingSpace->GetBinContent(binX, binY, binZ);
        votingSpace->SetBinContent(binX, binY, binZ, currentValue - 1);
    }

    TVector3 getHighestVote() {
        int binX, binY, binZ;
        votingSpace->GetMaximumBin(binX, binY, binZ);
        double x = votingSpace->GetXaxis()->GetBinCenter(binX);
        double y = votingSpace->GetYaxis()->GetBinCenter(binY);
        double z = votingSpace->GetZaxis()->GetBinCenter(binZ);
        return TVector3(x, y, z);
    }
};

void detectLines(PointCloud& pointCloud, double dx) {
    // 投票空間を初期化
    TVector3 minP(-10, -10, -10); // 点群の範囲（仮定）
    TVector3 maxP(10, 10, 10);
    HoughSpace hough(dx, maxP.X(), maxP.Y(), maxP.Z());

    // 点群を投票空間に追加
    for (const auto& point : pointCloud.points) {
        hough.addVote(point);
    }

    // 最も得票の多い直線を検出
    TVector3 highestVote = hough.getHighestVote();
    std::cout << "Highest voted line passes through: "
              << highestVote.X() << ", " << highestVote.Y() << ", " << highestVote.Z() << std::endl;

    // 投票結果の可視化
    TCanvas* canvas = new TCanvas("HoughSpace", "Hough Voting Space", 800, 600);
    hough.votingSpace->Draw("LEGO");
    canvas->SaveAs("hough_space.png");
}

int main() {
    // 点群データを生成
    PointCloud pointCloud;
    TRandom3 randGen;
    for (int i = 0; i < 100; i++) {
        pointCloud.addPoint(randGen.Uniform(-10, 10),
                            randGen.Uniform(-10, 10),
                            randGen.Uniform(-10, 10));
    }

    // ホフ変換を適用
    double dx = 0.5; // グリッド間隔
    detectLines(pointCloud, dx);

    return 0;
}

//-------------------------------------------------------------------------------------------------------------

//example 2
#include "TH2F.h"
#include "TCanvas.h"
#include <cmath>
#include <vector>
#include <tuple>
#include <iostream>

// Hough変換のしきい値
const double THRESHOLD = 100;

// エッジデータを持つヒストグラムを生成する関数（ダミーデータ）
TH2F* GenerateEdgeData(int sizeX, int sizeY) {
    TH2F* edgeData = new TH2F("edgeData", "Edge Data", sizeX, 0, sizeX, sizeY, 0, sizeY);
    // サンプルデータとして直線を描く
    for (int x = 20; x < 80; ++x) {
        edgeData->Fill(x, 2 * x + 10);  // 直線 y = 2x + 10
    }
    return edgeData;
}

// Hough空間を生成する関数
TH2F* PerformHoughTransform(TH2F* edgeData, int thetaBins, int rhoBins) {
    // Hough空間を初期化
    double maxRho = std::sqrt(std::pow(edgeData->GetXaxis()->GetXmax(), 2) +
                              std::pow(edgeData->GetYaxis()->GetXmax(), 2));
    TH2F* houghSpace = new TH2F("houghSpace", "Hough Space", thetaBins, -90, 90, rhoBins, -maxRho, maxRho);

    // エッジデータをHough空間に投票
    for (int xBin = 1; xBin <= edgeData->GetNbinsX(); ++xBin) {
        for (int yBin = 1; yBin <= edgeData->GetNbinsY(); ++yBin) {
            if (edgeData->GetBinContent(xBin, yBin) > 0) {  // エッジ画素のみ処理
                double x = edgeData->GetXaxis()->GetBinCenter(xBin);
                double y = edgeData->GetYaxis()->GetBinCenter(yBin);
                for (int thetaBin = 1; thetaBin <= thetaBins; ++thetaBin) {
                    double theta = houghSpace->GetXaxis()->GetBinCenter(thetaBin) * M_PI / 180.0;  // 度→ラジアン
                    double rho = x * std::cos(theta) + y * std::sin(theta);
                    houghSpace->Fill(theta * 180.0 / M_PI, rho);  // 度単位で投票
                }
            }
        }
    }
    return houghSpace;
}

// Hough空間からピークを検出する関数
std::vector<std::tuple<double, double>> FindHoughPeaks(TH2F* houghSpace, double threshold) {
    std::vector<std::tuple<double, double>> peaks;
    for (int thetaBin = 1; thetaBin <= houghSpace->GetNbinsX(); ++thetaBin) {
        for (int rhoBin = 1; rhoBin <= houghSpace->GetNbinsY(); ++rhoBin) {
            double value = houghSpace->GetBinContent(thetaBin, rhoBin);
            if (value > threshold) {  // しきい値を超えた場所をピークとして検出
                double theta = houghSpace->GetXaxis()->GetBinCenter(thetaBin);
                double rho = houghSpace->GetYaxis()->GetBinCenter(rhoBin);
                peaks.emplace_back(theta, rho);
            }
        }
    }
    return peaks;
}

// メイン関数
void HoughTransformExample() {
    // 入力データ（エッジデータ）を生成
    TH2F* edgeData = GenerateEdgeData(100, 100);

    // Hough変換を実行
    TH2F* houghSpace = PerformHoughTransform(edgeData, 180, 100);

    // Hough空間からピークを検出
    auto peaks = FindHoughPeaks(houghSpace, THRESHOLD);

    // 結果を表示
    TCanvas* c1 = new TCanvas("c1", "Hough Transform", 800, 600);
    c1->Divide(2, 1);
    c1->cd(1);
    edgeData->Draw("COLZ");  // エッジデータの可視化
    c1->cd(2);
    houghSpace->Draw("COLZ");  // Hough空間の可視化

    // 検出されたピークを表示
    for (const auto& peak : peaks) {
        double theta = std::get<0>(peak);
        double rho = std::get<1>(peak);
        std::cout << "Detected line: θ = " << theta << " degrees, ρ = " << rho << std::endl;
    }
}
