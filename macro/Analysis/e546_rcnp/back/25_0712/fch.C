{
	TChain *fc = new TChain("tree");
	
	TString ss;

	for(int i = 28; i <= 109; i++){
		ss = TString::Format("~/fair_install/transfer_data/run_%04d.root", i);
		cout<<ss.Data()<<endl;
		fc->Add(ss.Data());	
	}	
}
