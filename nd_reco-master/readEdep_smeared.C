#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TBits.h>
#include <TObjString.h>
#include <string>


using std::cout;
using std::endl;
void readEdep_smeared(TString filename){
  
  TFile file(filename);
  TTree * tree = (TTree *) file.Get("edep_smeared_tree");
  
  const int kNPmax = 300;

  int IEntry;
  int NPar;
  int NPrim;
  double RecoP4[kNPmax][4];
  double TrueP4[kNPmax][4];
  double RecoNuP4[4];
  double TrueNuP4[4];
  int  Pdg[kNPmax];
  int  TrackId[kNPmax];
  int  ParentId[kNPmax];
  int  TopParentId[kNPmax];
  double   Len[kNPmax];
  int  NXhit[kNPmax];
  int  NYhit[kNPmax];
  char      Info[kNPmax][10];
  //  TObjString *Info[kNPmax];
  //  std::string Info;
  TBranch *brIEntry=tree->GetBranch("IEntry");
  TBranch *brNPar= tree->GetBranch("NPar");
  TBranch *brNPrim= tree->GetBranch("NPrim");
  TBranch *brRecoP4= tree->GetBranch("RecoP4");
  TBranch *brTrueP4= tree->GetBranch("TrueP4");
  TBranch *brRecoNuP4= tree->GetBranch("RecoNuP4");
  TBranch *brTrueNuP4= tree->GetBranch("TrueNuP4");

  TBranch *brPdg=tree->GetBranch("Pdg");
  TBranch *brTrackId=tree->GetBranch("TrackId");
  TBranch *brParentId=tree->GetBranch("ParentId");
  TBranch *brTopParentId=tree->GetBranch("TopParentId");
  TBranch *brLen=tree->GetBranch("Len");
  TBranch *brNXhit=tree->GetBranch("NXhit");
  TBranch *brNYhit=tree->GetBranch("NYhit");
  TBranch *brInfo=tree->GetBranch("Info");
  
  
  brIEntry-> SetAddress ( &IEntry);
  brNPar-> SetAddress ( &NPar);
  brNPrim-> SetAddress ( &NPrim);
  brRecoP4-> SetAddress (RecoP4);
  brTrueP4-> SetAddress (TrueP4);
  brRecoNuP4-> SetAddress (RecoNuP4);
  brTrueNuP4-> SetAddress (TrueNuP4);
  brPdg-> SetAddress ( Pdg);
  brTrackId-> SetAddress ( TrackId);
  brParentId-> SetAddress ( ParentId);
  brTopParentId->SetAddress(TopParentId);
  brLen-> SetAddress ( Len);
  brNXhit-> SetAddress ( NXhit);
  brNYhit-> SetAddress ( NYhit);
  brInfo -> SetAddress ( Info);
  

  for(int i=0; i < 5; i++) {
    tree->GetEntry(i);
    printf("\n ----------------------------------------------------------------------------------");
    //  std::cout<<"info:"<<Info<<std::endl;
    printf("\n %d ientry:%d   nprim: %5d", i, IEntry,NPrim);
    printf("\n | %8.1f | %8.1f | %8.1f | %8.1f || %8.1f | %8.1f | %8.1f | %8.1f |", RecoNuP4[0], RecoNuP4[1],RecoNuP4[2], RecoNuP4[3], TrueNuP4[0],TrueNuP4[1],TrueNuP4[2], TrueNuP4[3]);
    printf("\n");
    printf("\n |   pdg  | trkid paren topPar|  recoPx    recoPy     recoPz     recoE  ||    truePx     truePy     truePz      trueE || length    nxhit   nyhit | INFO");
    for(int ip=0; ip<NPar; ip++) {
      std::string s;
      for(int j=0;j<10;j++){
	s+=Info[ip][j];
      }
      printf("\n | %6d | %3d | %3d | %3d | %8.1f | %8.1f | %8.1f | %8.1f ||  %8.1f | %8.1f | %8.1f | %8.1f || %7.1f | %5d | %5d | %s",
	     Pdg[ip], TrackId[ip], ParentId[ip], TopParentId[ip],
	     RecoP4[ip][0],RecoP4[ip][1],RecoP4[ip][2], RecoP4[ip][3],
	     TrueP4[ip][0],TrueP4[ip][1],TrueP4[ip][2],TrueP4[ip][3],
	     Len[ip], NXhit[ip], NYhit[ip], s.c_str()
	     );
      if(std::isnan(RecoP4[ip][0])) std::cout<<" nan exist"<<std::endl;
    }
  }
  printf("\n");

}
