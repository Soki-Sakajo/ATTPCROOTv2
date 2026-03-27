#include <TFile.h>
#include <TTree.h>

void extractEvents(const char* inputFileName, const char* outputFileName, const char* treeName, Long64_t eventNumber) {
    // 入力ファイルを開く
    TFile* inputFile = new TFile(inputFileName, "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open input file: " << inputFileName << std::endl;
        return;
    }

    // 入力ファイルからTreeを取得
    TTree* inputTree = dynamic_cast<TTree*>(inputFile->Get(treeName));
    if (!inputTree) {
        std::cerr << "Error: Could not find tree with name " << treeName << " in input file." << std::endl;
        inputFile->Close();
        return;
    }

    // 出力ファイルを作成
    TFile* outputFile = new TFile(outputFileName, "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Could not create output file: " << outputFileName << std::endl;
        inputFile->Close();
        return;
    }

    // 入力Treeから指定されたイベントを抜き出して出力Treeに追加
    outputTree = inputTree->CloneTree(0);  // 0は空のTreeを作成

    if (eventNumber >= 0 && eventNumber < inputTree->GetEntries()) {
        inputTree->GetEntry(eventNumber);
        outputTree->Fill();
    }

    // 出力Treeを保存
    outputFile->Write();

    // ファイルを閉じる
    inputFile->Close();
    outputFile->Close();
}

int main() {
    const char* inputFileName = "input.root";
    const char* outputFileName = "output.root";
    const char* treeName = "yourTreeName";
    Long64_t eventNumberToExtract = 42;  // 抜き出したいイベントの番号

    extractEvents(inputFileName, outputFileName, treeName, eventNumberToExtract);

    return 0;
}
