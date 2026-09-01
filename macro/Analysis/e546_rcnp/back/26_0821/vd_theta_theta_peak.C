#include <TF1.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TGraphErrors.h>

void vd_theta_theta_peak(){

    //setdata
    /*
    vector<pair<Double_t, Double_t>> peak1_data={
        {3.06,95.0},{3.50,88.0},{3.89,81.0}, {4.00,80.0},{4.10,78.0},{4.20,77.0},{4.50,73.0}
    };
    vector<pair<Double_t, Double_t>> peak2_data={
        {3.06,101.0},{3.50,93.0},{3.89,87.0}, {4.00,86.0},{4.10,84.0},{4.20,83.0},{4.50,79.0}
    };
    vector<pair<Double_t, Double_t>> peak3_data={
        {3.06,106.0},{3.50,98.0},{3.89,92.5}, {4.00,91.0},{4.10,89.0},{4.20,88.0},{4.50,85.0}
    };
    */

    vector<pair<Double_t, Double_t>> peak1_data={
        {3.89,81.0}, {4.00,80.0},{4.10,78.0},{4.20,77.0}
    };
    vector<pair<Double_t, Double_t>> peak2_data={
        {3.89,87.0}, {4.00,86.0},{4.10,84.0},{4.20,83.0}
    };
    vector<pair<Double_t, Double_t>> peak3_data={
        {3.89,92.5}, {4.00,91.0},{4.10,89.0},{4.20,88.0}
    };

    Int_t n1 = peak1_data.size();
    Int_t n2 = peak2_data.size();
    Int_t n3 = peak3_data.size();

    vector<vector<Double_t>> peak1(2,vector<Double_t>(n1,0));
    vector<vector<Double_t>> peak2(2,vector<Double_t>(n2,0));
    vector<vector<Double_t>> peak3(2,vector<Double_t>(n3,0));

    //like fill for TGraph
    if (n1 != n2 || n2 != n3 || n3 != n1){
        std::cout << std::endl;
        std::cout << "Something is strange. need to check." << std::endl;
        gSystem->Exit(0);
    }
    for(Int_t i = 0; i < n1; i++){
        peak1.at(0).at(i) = peak1_data[i].first;
        peak1.at(1).at(i) = peak1_data[i].second;
        peak2.at(0).at(i) = peak2_data[i].first;
        peak2.at(1).at(i) = peak2_data[i].second;
        peak3.at(0).at(i) = peak3_data[i].first;
        peak3.at(1).at(i) = peak3_data[i].second;
    }
    std::cout << "Filled data" << std::endl;

    //set TGraph
    TGraph *h_peak1 = new TGraph(n1, peak1.at(0).data(), peak1.at(1).data());
    TGraph *h_peak2 = new TGraph(n2, peak2.at(0).data(), peak2.at(1).data());
    TGraph *h_peak3 = new TGraph(n3, peak3.at(0).data(), peak3.at(1).data());

    //set line
    TF1 *f1 = new TF1("f1", "[0]+[1]*x", 3.0, 5.0);
    TF1 *f2 = new TF1("f2", "[0]+[1]*x", 3.0, 5.0);
    TF1 *f3 = new TF1("f3", "[0]+[1]*x", 3.0, 5.0);
    TF1 *y90 = new TF1("y90", "90", 3.0, 5.0);

    //fit
    h_peak1->Fit(f1, "R");
    h_peak2->Fit(f2, "R");
    h_peak3->Fit(f3, "R");

    //set style
    h_peak1->SetMarkerStyle(22);
    h_peak1->SetMarkerSize(1.2);
    h_peak1->SetMarkerColor(kGreen+2);
    h_peak1->SetLineColor(kGreen+2);
    f1->SetNpx(1000);
    f1->SetLineColor(kGreen+2);
    h_peak2->SetMarkerStyle(21);
    h_peak2->SetMarkerSize(1.2);
    h_peak2->SetMarkerColor(kBlue);
    h_peak2->SetLineColor(kBlue);
    f2->SetNpx(1000);
    f2->SetLineColor(kBlue);
    h_peak3->SetMarkerStyle(20);
    h_peak3->SetMarkerSize(1.2);
    h_peak3->SetMarkerColor(kRed);
    h_peak3->SetLineColor(kRed);
    f3->SetNpx(1000);
    f3->SetLineColor(kRed);
    y90->SetLineColor(kBlack);
    y90->SetLineWidth(3);

    //draw
    TCanvas *c = new TCanvas("c","c");
    //graph
    h_peak1->GetXaxis()->SetNdivisions(510);
    h_peak1->GetYaxis()->SetNdivisions(510);
    h_peak1->Draw("AP");
    h_peak2->Draw("P SAME");
    h_peak3->Draw("P SAME");
    //line
    f1->Draw("SAME");
    f2->Draw("SAME");
    f3->Draw("SAME");
    y90->Draw("same");
    //legend
    TLegend *le1 = new TLegend(0.7,0.7,0.9,0.9);
    /*
    le1 -> AddEntry(h_peak1, "peak1", "p");
    le1 -> AddEntry(h_peak2, "peak2", "p");
    le1 -> AddEntry(h_peak3, "peak3", "p");
    */
    le1 -> AddEntry(h_peak3, "gs-gs", "p");
    le1 -> AddEntry(h_peak2, "gs-ex", "p");
    le1 -> AddEntry(h_peak1, "ex-ex", "p");
    le1 -> Draw();
    //axis range
    //    h_peak1->GetXaxis()->SetLimits(3.0,5.0);
    //    h_peak1->GetYaxis()->SetRangeUser(70.0,110.0);
    h_peak1->GetXaxis()->SetLimits(3.8,4.3);
    h_peak1->GetYaxis()->SetRangeUser(75.0,95.0);

    c -> Update();
    std::cout << std::endl;

    //calculate y = 90
    Double_t y0 = 90;
    Double_t a1 = f1 -> GetParameter(0);
    Double_t b1 = f1 -> GetParameter(1);
    Double_t a2 = f2 -> GetParameter(0);
    Double_t b2 = f2 -> GetParameter(1);
    Double_t a3 = f3 -> GetParameter(0);
    Double_t b3 = f3 -> GetParameter(1);

    if(b1 == 0 || b2 == 0 || b3 == 0){
        std::cout << "Error: slope is 0" << std::endl;
        gSystem -> Exit(0);
    }

    Double_t x1 = (y0 - a1) / b1;
    Double_t x2 = (y0 - a2) / b2;
    Double_t x3 = (y0 - a3) / b3;

    std::cout << std::setprecision(4);
    std::cout << "y = 90, x1: " << x1 << ", x2: " << x2 << ", x3: " << x3 << std::endl;
    std::cout << std::endl;

}