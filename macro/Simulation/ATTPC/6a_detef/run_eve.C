/*#include "TString.h"
#include "AtEventDrawTask.h"
#include "AtEventManager.h"

#include "FairLogger.h"
#include "FairParRootFileIo.h"
#include "FairRunAna.h"
*/

//void run_eve(int runNum = 1, TString OutputDataFile = "./data1/output.reco_display.root")
void run_eve(int subnum = 0)
{
  //   TString InputDataFile = TString::Format("./data/run_%04d.root", runNum);
  //  TString InputDataFile = TString::Format("./data1/output_digi_rcnp_sequential_hoylehoyle_80matm_"+to_string(subnum)+".root");
  TString infname="./data1/output_digi_rcnp_sequential_hoylehoyle_80matm_"+to_string(subnum)+".root";
  //  TString InputDataFile = TString::Format(infname);
  TString InputDataFile = infname;
   std::cout << "Opening: " << InputDataFile << std::endl;

   TString dir = getenv("VMCWORKDIR");
   //   TString geoFile = "ATTPC_v1.1_geomanager.root";
   TString geoFile = "rcnp_80matm_geomanager.root";
   //   TString mapFile = "e12014_pad_mapping.xml";
   TString mapFile = "e12014_pad_mapping.xml";

   TString OutputDataFile = "./data1/output.reco_display.root";
   TString InputDataPath = InputDataFile;
   TString OutputDataPath = OutputDataFile;
   TString GeoDataPath = dir + "/geometry/" + geoFile;
   TString mapDir = dir + "/scripts/" + mapFile;

   FairRunAna *fRun = new FairRunAna();
   FairRootFileSink *sink = new FairRootFileSink(OutputDataFile);
   FairFileSource *source = new FairFileSource(InputDataFile);
   fRun->SetSource(source);
   fRun->SetSink(sink);
   fRun->SetGeomFile(GeoDataPath);

   FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
   FairParRootFileIo *parIo1 = new FairParRootFileIo();
   // parIo1->open("param.dummy.root");
   rtdb->setFirstInput(parIo1);

   auto fMap = std::make_shared<AtTpcMap>();
   fMap->ParseXMLMap(mapDir.Data());
   AtViewerManager *eveMan = new AtViewerManager(fMap);

   auto tabMain = std::make_unique<AtTabMain>();
   tabMain->SetMultiHit(100); // Set the maximum number of multihits in the visualization
   eveMan->AddTab(std::move(tabMain));

   eveMan->Init();

   std::cout << "Finished init" << std::endl;
   // eveMan->RunEvent(27);
}
