//future improvement:
//NOTE: better external pi0 stats, beta stats
//NOTE: ecal energy smearing equation only works for em shower, but currently used for any particles failed STT smearing

// k0->k0L, primary will break
// k0L-> , primary not break
// lambda -> , primary will break ( only for gamma, for charged particle, their primary are still lambda)
// eta--> will break
// sigma0 --> will break
// sigma- --> will break for gamma(and other particles), for neutron not break

// decay: pi+ -> mu+
// decay: mu+ -> e+
// decay : kaon- -> e- + pi0
// decay: pi- -> gamma + gamma (hadron inelastic)



// TEST MACRO FOR USE WITH OLDER ROOT6.  DOESN"T WORK WHEN CLING KNOWS ABOUT
// THE VARIOUS CLASSES.
#include <TFile.h>
#include "TSystem.h"
#include <TGeoManager.h>
#include <iostream>
#include <TRandom3.h>
#include "TH1.h"
#include "TH2.h"
#include "TGeoTrd2.h"
#include "TGeoTube.h"
#include <TTree.h>
#include "TDatabasePDG.h"

#include "TG4Event.h"

#include <cstring>
#include <fstream>
#include <sstream>
#include <utility>
#include <functional>
#include <cassert>
#include <map>
#include <cstdlib>
#include <assert.h> 
using std::cout;
using std::endl;
TTree* gEDepSimTree = NULL;
TGeoManager *geo;
TG4Event* event;
double sigmas=200E-6; // m
double B=0.6; 
double x0=2.8; // m 
TRandom3 *ran;
TH2 *herr_dipAngle_stt;
TH2 *herr_dipAng100_stt;
TH2 *herr_thetaYZ100_stt;
TH2 *herr_pt_stt;
TH2 *herr_p_stt;
TH2 *herr_theta_stt;
TH1 *herr_E_ecal;
TH1 *herr_theta_ecal;
TH1 *herr_phi_ecal;
TH1 *herr_theta_N;
TH1 *herr_phi_N;
TH2 *herr_E_equa_N;
TH2 *herr_p_equa_N;
TH1 *herr_p_beta_N;
TH1 *herr_p_pi0;
TH1 *herr_theta_pi0;
TH1 *herr_phi_pi0;

TH1 *herr_nu_E;
double recoNuPx;
double recoNuPy;
double recoNuPz;
double trueNuE;
TFile *outf;
TFile *outTreeF;
TTree * tree;

int iEntry;
int targetpdg;
bool isHtarget;
int iChannel;
int debug;
TDatabasePDG *dbpdg;

TH2* hPi0_mom_recotrue;
TH2* hPi0_ang_recotrue;
TH1* hNeutron_ang_reso;
TH2* hNeutron_beta_recotrue_stt;
TH2* hNeutron_beta_recotrue_ecal;
//std::ofstream treefile;
std::map<int,std::pair<int,int> > sttMap;  // build only "the first" continuous chain by PrimaryId  -> for STT track building
//std::map<int,std::pair<int,int> > ecalMap;
//std::map<int,std::pair<int,int> > sttMap_prim;
//std::map<int,std::pair<int,int> > ecalMap_prim;
std::map<int, std::vector<int> > sttMap_prim2; // build ALL continuous chain by PrimaryId -> for pi0-decayed gamma finding 
std::map<int, std::vector<int> > ecalMap_prim2;  // build ALL continuous chain by PrimaryId -> for finding hits in ecal 
std::map<int, std::vector<int> > N_sttMap_contr_primIndex;  // build ALL continuous chain by contrib[0] -> for finding neutron all daughters' hits
std::map<int, std::vector<int> > N_ecalMap_contr_primIndex; // build ALL continuous chain by contrib[0] -> for finding neutron all daughters' hits

const int kNPmax = 300;
int         brIEntry;
int         brNPar;
int         brNPrim;
double      brRecoP4[kNPmax][4];
double      brTrueP4[kNPmax][4];
double  brRecoNuP4[4];
double  brTrueNuP4[4];
int         brPdg   [kNPmax];
int         brTrackId[kNPmax];
int         brParentId[kNPmax];
int         brTopParentId[kNPmax];
double      brLen[kNPmax];
int         brNXhit[kNPmax];
int         brNYhit[kNPmax];
char        brInfo[kNPmax][10];
int iFillPar;

void N_organizeHits_contr();
int findgammaPrimaryId(int trackid);
void smearPi0(int trackid);
void smearPar(int trackid, std::string name);
void findEvis_forCell(int starthit, int nhit, std::map<int, double> &cellId_Evis);
void findEvis_forCell(std::vector<int> allhits, std::map<int, double> &cellId_Evis);

bool inFV(double x, double y, double z){
  //  double centerX=0.;
  double centerY=-2384.73;  // mm
  double centerZ=23910; // mm

  if(abs(x)>1490) return false;
  double r=sqrt((y-centerY)*(y-centerY)+(z-centerZ)*(z-centerZ));
  if(r>1800) return false;
  return true;
}


void cleanBranch(){
  recoNuPx=0;
  recoNuPy=0;
  recoNuPz=0;
  brNPar=0;
  for(int i=0;i<kNPmax;i++){
    brRecoP4[i][0]=-999;
    brRecoP4[i][1]=-999;
    brRecoP4[i][2]=-999;
    brRecoP4[i][3]=-999;
    brTrueP4[i][0]=-999;
    brTrueP4[i][1]=-999;
    brTrueP4[i][2]=-999;
    brTrueP4[i][3]=-999;
    brPdg[i]=-999;
    brTrackId[i]=-999;
    brParentId[i]=-999;
    brTopParentId[i]=-999;
    brLen[i]=-999;
    brNXhit[i]=-999;
    brNYhit[i]=-999;
    // brInfo[i]=0;
    strcpy(brInfo[i],"");
  }  
}

void showAll(){
  std::cout<<"============================================="<<std::endl;
  for (std::vector<TG4Trajectory>::iterator
         t = event->Trajectories.begin();
       t != event->Trajectories.end(); ++t) {
    std::cout << "   Traj " << t->TrackId;
    std::cout << " " << t->ParentId;
    std::cout << " " << t->Name;
    std::cout << " " << t->Points.size();
    std::cout<< " E:"<<t->GetInitialMomentum().E();
    std::cout<<" beginpro:"<<t->Points.begin()->Process<<" "<<t->Points.begin()->Subprocess<<" endpro:"<<(t->Points.end()-1)->Process<<" "<<(t->Points.end()-1)->Subprocess;
    std::cout << std::endl;
  }
  for (auto d = event->SegmentDetectors.begin();
       d != event->SegmentDetectors.end(); ++d) {
    std::cout << "   det " << d->first;
    std::cout << " " << d->second.size();
    int count = 10;
    std::cout << " up to " << count << " segments";
    std::cout << std::endl;
    //    if(d->first!="Straw") continue;
    int i=0;
    for (std::vector<TG4HitSegment>::iterator
	   h = d->second.begin();
	 h != d->second.end();
	 ++h) {
      std::cout << "      "<<i;
      i++;
      std::cout << " P: " << h->PrimaryId<<" "<<h->Contrib[0];
      std::cout << " E: " << h->EnergyDeposit;
      std::cout << " S: " << h->SecondaryDeposit;
      std::cout << " C: " << h->Contrib.size()<<"->";
      for(unsigned long j=0;j<h->Contrib.size();j++){
	std::cout<<" "<<h->Contrib[j];
      }
      //      std::cout<<" name:"<<h->GetVolName();
      //            std::cout << " L: " << h->TrackLength;
      TLorentzVector mid= (h->Start+h->Stop)*0.5;
      TString name=geo->FindNode(mid.X(),mid.Y(),mid.Z())->GetName();
      std::cout<<" "<<name;
      std::cout<<" start:"<<h->Start.X()<<" "<<h->Start.Y()<<" "<<h->Start.Z()<<" "<<h->Start.T()<<" endT:"<<h->Stop.T();
      if((h+1)!= d->second.end() && (h+1)->Start.T()<h->Start.T()) std::cout<<"   !!!!!!! time reverted";
      std::cout<<std::endl;
    }
  }

  std::cout<<"============================================="<<std::endl;
}

int  findTopParent(int trackid){
  TG4Trajectory trk=event->Trajectories[trackid];
  while(trk.ParentId!=-1){
    trk=event->Trajectories[trk.ParentId];
  }
  return trk.TrackId;
}

void fill1Par2tree(double recoPx, double recoPy, double recoPz, int trackid, double len, int nXhit, int nYhit, const char info[10]){
  //  if(std::isnan(recoPx) || std::isnan(recoPy) || std::isnan(recoPz))  { std::cout<<"fill nan"; std::exit(EXIT_FAILURE);}
  if(recoPx==0 && recoPy==0 && recoPz==0) std::cout<<" recop=0 000000000000000000000000000000000000000000000, need check "<<std::endl;
  brRecoP4[iFillPar][0]=recoPx;
  brRecoP4[iFillPar][1]=recoPy;
  brRecoP4[iFillPar][2]=recoPz;
  int pdg=event->Trajectories[trackid].PDGCode;
  double m=dbpdg->GetParticle(pdg)->Mass()*1000;  // MeV
  brRecoP4[iFillPar][3]=sqrt(recoPx*recoPx+recoPy*recoPy+recoPz*recoPz+m*m);

  brTrueP4[iFillPar][0]=event->Trajectories[trackid].InitialMomentum.X();
  brTrueP4[iFillPar][1]=event->Trajectories[trackid].InitialMomentum.Y();
  brTrueP4[iFillPar][2]=event->Trajectories[trackid].InitialMomentum.Z();
  brTrueP4[iFillPar][3]=event->Trajectories[trackid].InitialMomentum.E();
  brTrackId[iFillPar]= trackid;
  brLen[iFillPar]= len;
  brNXhit[iFillPar]= nXhit;
  brNYhit[iFillPar]= nYhit;
  brParentId[iFillPar]= event->Trajectories[trackid].ParentId;
  brPdg[iFillPar]= pdg;
  brTopParentId[iFillPar]= findTopParent(trackid);
  strcpy(brInfo[iFillPar], info);
  iFillPar++;
  if(debug>2) std::cout<<"fill 1 par pdg:"<<pdg<<std::endl;

  recoNuPx+=recoPx;
  recoNuPy+=recoPy;
  recoNuPz+=recoPz;

}


class Node{
public:
  TG4Trajectory *Traj;
  Node *Parent;
  Node *FirstChild;
  Node *RightSibling;
  Node(TG4Trajectory *traj, Node *parent, Node *firstChild, Node *rightSibling): Traj(traj),Parent(parent),FirstChild(firstChild),RightSibling(rightSibling){}
};
Node *root;
std::vector<int> verlines;

Node *findNodeFast(int trackid){
  if(trackid==-1) return root;
  std::vector<int> parents;
  parents.push_back(trackid);
  TG4Trajectory *traj=&event->Trajectories[trackid];
  while(1){
    if(traj->ParentId==-1) break;
    parents.push_back(traj->ParentId);
    traj=&event->Trajectories[traj->ParentId];
  }
  Node *no=root->FirstChild;
  while(no){
    while(no){
      if(no->Traj->TrackId==parents.back()) { parents.pop_back(); break;}	
      no=no->RightSibling;
    }
    if(parents.empty()) return no;
    no=no->FirstChild;
  }
  return 0;
}
Node *findNode(Node *topNode, int trackid){
  if(trackid==-1) { /*std::cout<<"primary particle"<<std::endl; */return topNode;}
  if(topNode->Traj!=0) 
    if(trackid==topNode->Traj->TrackId) 
      return topNode;
  if(topNode->FirstChild==0) return 0;
  Node *next=topNode->FirstChild;
  while(1){    
    Node *m=findNode(next,trackid);
    if(m!=0) return m;
    if(next->RightSibling==0) break;
    next=next->RightSibling;
  }
  return 0;
}

void insertNode(Node *n){

  Node *parent=findNodeFast(n->Traj->ParentId);
  if(parent==0) std::cout<<"cannot find parent, something must be wrong"<<std::endl;
  n->Parent=parent;
  if(parent->FirstChild==0) {  parent->FirstChild=n; return;}
  Node *next=parent->FirstChild;
  while(1){
    if(next->RightSibling==0) { next->RightSibling=n; break;}
    next=next->RightSibling;    
  }
}

void makeTree(){
  root=new Node(0,0,0,0);
  for (unsigned long i=0;i<event->Trajectories.size();i++){
    insertNode(new Node(&(event->Trajectories[i]),0,0,0));
  }
}
bool findit(int i, std::vector<int> nums){
  for(unsigned long j=0;j<nums.size();j++){
    if(nums[j]==i) return true;
  }
  return false;
}


void drawNode(Node *n, int icol){
  if(icol==0) std::cout<<"_";
  else std::cout<<"__"<<n->Traj->Name<<"("<<n->Traj->TrackId<<")";
  if(n->FirstChild==0) return;
  Node *c=n->FirstChild;
  if(c->RightSibling!=0) { verlines.push_back(icol);}
  while(1){
    if(c->Parent->FirstChild!=c &&  c->RightSibling==0) {verlines.pop_back();}
    drawNode(c, icol+ 2+(c->Traj->Name).size()+ std::to_string(c->Traj->TrackId).size()+2);
    //    if(c->RightSibling==0) {std::cout<<"one pop"<<std::endl;verlines.pop_back();}
    if(c->RightSibling==0) break;
    c=c->RightSibling;
    std::cout<<std::endl;
    for(int i=0;i<=verlines.back();i++){
      if(findit(i,verlines)) std::cout<<"|";
      else std::cout<<" ";        
    }
  }

}

void drawTree(){
  verlines.clear();
  drawNode(root, 0);
  std::cout<<std::endl;
}
void dumpTree(){
  makeTree();
  drawTree();
}

void findEvis_forCell(std::vector<int> hitchains, std::map<int, double> &cellId_Evis){
  int id, cellID, planeID, modID;
  cellId_Evis.clear();
  assert(hitchains.size()%2==0);
  //  if(hitchains.size()%2!=0) std::exit(EXIT_FAILURE);
  for(unsigned int i=0;i<hitchains.size()/2;i++){
    for(int j = hitchains[2*i]; j <= hitchains[2*i+1]; j++){
      const TG4HitSegment& h = event->SegmentDetectors["ECAL"].at(j);
      double x = 0.5*(h.Start.X()+h.Stop.X());
      double y = 0.5*(h.Start.Y()+h.Stop.Y());
      double z = 0.5*(h.Start.Z()+h.Stop.Z());
      TGeoNode* node = geo->FindNode(x,y,z);
      TString slabstr = node->GetName();
      TString modstr=geo->GetMother()->GetName();
      //    std::cout<<"slabstr:"<<slabstr<<std::endl;
      if(slabstr.Contains("volECALActiveSlab") == true)
	{
	TObjArray* obj1 = slabstr.Tokenize("_");  //volECALActiveSlab_125_PV_0
	TObjArray* obj2 = modstr.Tokenize("_");  //ECAL_lv_PV_19

	int slabID;
	modID  = ((TObjString*) obj2->At(3))->GetString().Atoi();
	slabID = ((TObjString*) obj1->At(1))->GetString().Atoi();

	//    std::cout<<"modID:"<<modID<<" slabID:"<<slabID<<std::endl;
	delete obj1;
	delete obj2;
	// planeID==0 -> smallest slab
	// planeID==208 -> biggest slab
	planeID = slabID/40;

	if (planeID > 4) planeID = 4;
	double Pmaster[3]={x,y,z};
	double Plocal[3];
	geo->GetCurrentNavigator()->MasterToLocal(Pmaster,Plocal);
	//    std::cout<<"Plocal[0]:"<<Plocal[0]<<std::endl;  // along circular second longest
	//    std::cout<<"Plocal[1]:"<<Plocal[1]<<std::endl;  // along column longest
	//    std::cout<<"Plocal[2]:"<<Plocal[2]<<std::endl; //along radial, short
	TGeoTrd2* trd = (TGeoTrd2*) node->GetVolume()->GetShape();

	double dx1 = trd->GetDx1();  // shorter one along circumferential
	double dx2 = trd->GetDx2();  // longer one along circumferential
	double dz  = trd->GetDz();   // half thickness , along radial
	//	double dy1 = trd->GetDy1();  // along axial/fiber , same as dy2
	//	d1 = dy1 + Plocal[1];
	//	d2 = dy1 - Plocal[1];
	double dx = (dx2 - dx1) / dz * Plocal[2];
	double dis= Plocal[0]>0? (dx1+dx2)/2. + Plocal[0] - dx/2.: (dx1+dx2)/2. + Plocal[0] + dx/2.;
	double cellw = (dx1+dx2) / 12.;
	cellID = dis / cellw;
      }
    else if(slabstr.Contains("vol_endECALActiveSlab") == true)
      {

	TObjArray* obja = slabstr.Tokenize("_");
	int slabID;
	modID  = x>0?30:40;
	slabID = ((TObjString*) obja->At(1))->GetString().Atoi();

	delete obja;

	planeID = (208 - slabID)/40;

	if (planeID > 4) planeID = 4;

        double Pmaster[3]={x,y,z};
	double Plocal[3];
	//      std::cout<<"x:"<<x<<" y:"<<y<<" z:"<<z<<std::endl;
	geo->GetCurrentNavigator()->MasterToLocal(Pmaster,Plocal);

	TGeoTube* tub = (TGeoTube*) node->GetVolume()->GetShape();
	double rmax = tub->GetRmax();
	//      double dz  = tub->GetDz();
	// Plocal[0] : horizontal distance to tube center
	// Plocal[1] : vertical distance
	//	d1 = rmax * TMath::Sin(TMath::ACos(Plocal[0]/rmax)) - Plocal[1]; // d1 is shorter distance here
	//	d2 = rmax * TMath::Sin(TMath::ACos(Plocal[0]/rmax)) + Plocal[1];

	cellID = int((Plocal[0]/rmax + 1.) * 45);
      }
    else { continue;}
    id = cellID + 100 * planeID + 1000 * modID;
    //    std::cout<<"cellID:"<<cellID<<" planeID:"<<planeID<<" modID:"<<modID<<std::endl;
    if(cellId_Evis.find(id)!=cellId_Evis.end()) cellId_Evis[id]+= h.EnergyDeposit;
    else cellId_Evis[id]= h.EnergyDeposit;
    
    } // for 
  }

}

void smearFirstHitPosition(int starthit, double sigmaX, double sigmaY, double *pos_smear){
  sigmaX*=10.; //cm -> mm
  sigmaY*=10.; //cm --> mm
  //  std::cout<<"start smearFirstHitPosition starthit"<<starthit<<std::endl;
  const TG4HitSegment& h = event->SegmentDetectors["ECAL"].at(starthit);
  double x = 0.5*(h.Start.X()+h.Stop.X());
  double y = 0.5*(h.Start.Y()+h.Stop.Y());
  double z = 0.5*(h.Start.Z()+h.Stop.Z());
  pos_smear[0]=x;
  pos_smear[1]=y;
  pos_smear[2]=z;
  TGeoNode* node = geo->FindNode(x,y,z);
  double Pmaster[3]={x,y,z};
  double Plocal[3];

  geo->GetCurrentNavigator()->MasterToLocal(Pmaster,Plocal);
  ////// barrel ///////
  ////// barrel ///////

  TString slabstr = node->GetName();
  if(slabstr.Contains("volECALActiveSlab") == true) {
    double Plocal_smear[3]= { Plocal[0]+ran->Gaus(0,sigmaX), Plocal[1]+ran->Gaus(0,sigmaY), Plocal[2]};
    //////////////////////// boundary check 
    TGeoTrd2* trd = (TGeoTrd2*) node->GetVolume()->GetShape();
    double dx1 = trd->GetDx1();  // shorter one along circumferential
    double dx2 = trd->GetDx2();  // longer one along circumferential
    double dz  = trd->GetDz();   // half thickness , along radial
    double dy1 = trd->GetDy1();  // along axial/fiber , same as dy2
    //    std::cout<<"dx1:"<<dx1<<" dx2:"<<dx2<<" dz:"<<dz<<" dy1:"<<dy1<<std::endl;
    double maxX= (dx1+dx2)/2. + Plocal[2]*(dx2-dx1)/dz;
    if(Plocal_smear[0] > maxX) Plocal_smear[0]=maxX;
    else if(Plocal_smear[0] < -maxX) Plocal_smear[0]=-maxX;
    if(Plocal_smear[1]> dy1 ) Plocal_smear[1]= dy1;
    else   if(Plocal_smear[1]< -dy1 ) Plocal_smear[1]=-dy1;
    geo->GetCurrentNavigator()->LocalToMaster(Plocal_smear, pos_smear);
  }
  else if(slabstr.Contains("vol_endECALActiveSlab") == true){
    ////// need to know which simga to use for endcap ???????? ////////
    ////// now temperatly use same as barrel, sigma ??????? /////////
    double Plocal_smear[3]= { Plocal[0]+ran->Gaus(0,sigmaX), Plocal[1]+ran->Gaus(0,sigmaY), Plocal[2]};
    // boundary check 
    TGeoTube* tub = (TGeoTube*) node->GetVolume()->GetShape();
    double rmax = tub->GetRmax();
    double width=rmax*2./90.;
    double low=floor(Plocal[0]/width)*width;
    double high=ceil(Plocal[0]/width)*width;
    if(Plocal_smear[0]>high) Plocal_smear[0]=high;
    else     if(Plocal_smear[0]<low) Plocal_smear[0]=low;
    double ymax= sqrt(rmax*rmax+low*low);
    if(ymax<sqrt(rmax*rmax+high*high)) ymax=sqrt(rmax*rmax+high*high);
    if(Plocal_smear[1]>ymax) Plocal_smear[1]=ymax;
    else if(Plocal_smear[1] <-ymax) Plocal_smear[1]=-ymax;
    geo->GetCurrentNavigator()->LocalToMaster(Plocal_smear, pos_smear);
  }
  else { std::cout<<"this hit is neither in endcap or barrel of ecal, this must be very rare, simply return true "<<slabstr<<std::endl;}
  return;
  
}

bool isthisParent(int lowId, int highId){
  if(lowId==highId) return true;
  int parentId=event->Trajectories[lowId].ParentId;
  while(parentId!=-1){
    if(parentId==highId) return true;    
    parentId=event->Trajectories[parentId].ParentId;
  }
  return false;
}

void bruteforceFindHit_ecal(int trackid, int primid, std::vector<int> &allhits, bool havechild=true){
  allhits.clear();
  assert(primid!=-1);
  assert(trackid!=primid);
  std::vector<int> hitchains=ecalMap_prim2[primid];
  if(havechild){
    for(unsigned int i=0;i<hitchains.size()/2;i++){
      for(int j = hitchains[2*i]; j <= hitchains[2*i+1]; j++){
	if( !isthisParent(event->SegmentDetectors["ECAL"].at(j).Contrib[0], trackid)) { continue;}
	//      std::cout<<"j:"<<j<<" check1pass"<<std::endl;
	if(allhits.size()==0) {allhits.push_back(j);allhits.push_back(j);}
	else if(j-allhits.back()==1) allhits.back()=j;
	else
	  {allhits.push_back(j);allhits.push_back(j);}	
      }
    }
  }
  else{
    for(unsigned int i=0;i<hitchains.size()/2;i++){
      for(int j = hitchains[2*i]; j <= hitchains[2*i+1]; j++){
	if(event->SegmentDetectors["ECAL"].at(j).Contrib[0]!=trackid) { continue;}
	//      std::cout<<"j:"<<j<<" check1pass"<<std::endl;
        if(allhits.size()==0) {allhits.push_back(j);allhits.push_back(j);}
        else if(j-allhits.back()==1) allhits.back()=j;
        else
          {allhits.push_back(j);allhits.push_back(j);}
      }
    }
  }
  /*
  std::cout<<"bruteforce allhits pairs: ";
  for(auto a: allhits){
    std::cout<<" "<<a;
  }
  std::cout<<std::endl;
  */
}
int findPrimaryId(int trackid){
  
  int parentId=event->Trajectories[trackid].ParentId;
  if(event->Trajectories[parentId].Name=="gamma") return findgammaPrimaryId(parentId);
  else { std::cout<<"findPrimaryId check!!!!, trackid:"<<trackid<<" parent:"<<event->Trajectories[parentId].Name<<" iEntry:"<<iEntry<<std::endl; return parentId;}
}

bool smearPar_ecal(int trackid, int primaryId=-1){
  if(debug>=2) std::cout<<"start ecal smear for trackid:"<<trackid<<std::endl;
  //  5.7%/sqrt(E)  (GeV)
  int starthit;
  //  int nhit;  
  std::vector<int> allhits;
  std::map<int, double> cellId_Evis;
  std::string parentName=event->Trajectories[trackid].ParentId==-1?"empty":event->Trajectories[event->Trajectories[trackid].ParentId].Name;
  if(event->Trajectories[trackid].ParentId==-1 || primaryId==trackid ||  parentName=="pi0"){
    if(debug>=1) std::cout<<"didnot choose bruteforce since it is its own primaryId, trackid:"<<trackid<<" iEntry:"<<iEntry<<std::endl;
    if(ecalMap_prim2.find(trackid)==ecalMap_prim2.end()) return false;
    findEvis_forCell(ecalMap_prim2[trackid], cellId_Evis);
    starthit=ecalMap_prim2[trackid][0];
  }
  else {
    
    if(primaryId==-1) { primaryId= findPrimaryId(trackid); if(debug>=1) std::cout<<"@@@@ lookforprim ";}
    if(debug>=1) std::cout<<"---> start bruteforceFindHit for trackid: "<<trackid<<" "<<event->Trajectories[trackid].Name<<" havechild:"<<(findNodeFast(trackid)->FirstChild!=0)<<" prim:"<<primaryId<<std::endl;
    bruteforceFindHit_ecal(trackid, primaryId,  allhits, findNodeFast(trackid)->FirstChild!=0);
    //    std::cout<<"-------> after bruteforceFindHit"<<std::endl;
    if(allhits.size()==0) {  return false;}

    findEvis_forCell(allhits, cellId_Evis);
    starthit=allhits[0]; 
  }
  bool haveDetectableHit=false;
  for(auto cellE: cellId_Evis){
    if(cellE.second > 0.1) {
      haveDetectableHit=true;
      break;
    }
  }
  
  if(!haveDetectableHit) { 
    //    std::cout<<"ecal smear failed due to no detectable ecal hit"<<std::endl; 
    return false;}
  double E=event->Trajectories[trackid].InitialMomentum.E();
  double de2e= 0.057/sqrt(E/1000.);
  double E_smear= E*ran->Gaus(1,de2e); // MeV
  int pdg=event->Trajectories[trackid].PDGCode;
  double m=dbpdg->GetParticle(pdg)->Mass() * 1000; // MeV
  while(E_smear<=m){
    E_smear= E*ran->Gaus(1,de2e);
  }

  double P_smear=E_smear>m?sqrt(E_smear*E_smear-m*m):0; // MeV
  //  std::cout<<"E:"<<E<<" E_smear:"<<E_smear<<" P_smear:"<<P_smear<<std::endl;
  double p0                        =    0.0034793;
  double p1                        =   -0.0151348;
  double p2                        =      1.52382;
  double p3                        =      0.57308;  
  double q0                        =     -1.93866;
  double q1                        =      13.5211;
  // sigmaX fitting function : [0]*x*exp([1]*x+[2])+[3]"
  // sigma Z fitting function:  [0]*log(x)+[1]
  //sigmaX ,sigmaZ, plot extracted from paper <<A FLUKA simulation of the KLOE electromagnetic calorimeter>> 
  // fit on my own 

  double sigmaX= p0*E*exp(p1*E+p2)+p3;    // cm  along circumferential
  double sigmaZ= q0*log(E)+q1;  // cm  along fiber/axial 
  double Pos_smear[3];
  smearFirstHitPosition( starthit, sigmaX, sigmaZ, Pos_smear); // sigmaZ is actually sigmaY in local coordinate of Ttrd2

  TVector3 firstHitPos_smear(Pos_smear);
  TVector3 vtx=event->Primaries[0].Position.Vect(); 
  TVector3 dir=firstHitPos_smear-vtx;  
  TVector3 P3_smear=P_smear*dir.Unit();
  //  std::cout<<"ecal smear succeeded, fill it"<<std::endl;
  
  double theta=event->Trajectories[trackid].InitialMomentum.Theta();
  double phi=event->Trajectories[trackid].InitialMomentum.Phi();
  double theta_smear=dir.Theta();
  double phi_smear=dir.Phi();
  herr_E_ecal->Fill((E_smear-E)/E*100);
  herr_theta_ecal->Fill((theta_smear-theta)/theta*100);
  herr_phi_ecal->Fill((phi_smear-phi)/phi*100);
  
  fill1Par2tree(P3_smear.X(), P3_smear.Y(), P3_smear.Z(),  trackid, -999, -999, -999, "ecalsmear ");
  return true;
}

bool smearChargedPar_stt(int trackid){
  //  std::cout<<"stt smearing start --------------"<<trackid<<std::endl;
  // Ptran RMS from the equation, 
  // the angle between Ptran and Px smear by PDG multiscattering-RMS-equation, in which Px is decided.  
  // angle between Py and Pz is smeared by Roberto-provided equation with multiple-scattering(second) term replaced by PDG one
  //  std::cout<<"start smearChargedPar_stt:"<<trackid<<std::endl;
  int nYhit=0;
  int nXhit=0;
  double Lyz=0;
  double L=0;
  double Lx=0;
  //EDEP FACT: for hitsegment, the primary Id will be its top parent's track id, unless if there's "decay" process happens, the decayed daughters will start to become top parent
  //EDEP FACT: same track's higsegment most time are id-connected,form a long chain,if later broke, only a few hits left,  thoese later hits are wierd hits, themselves are not time-connected, space-connected, better disregard.   for the long chain, most time they are time-forwarding, but the last a few hits may time-reverse, these hits are also wierd hits. 

  // EDEP FACT: for pion0->gamma, gamma will use its own trackid as primaryid, but gamma->e+ e-, e+ e- will use gammas trackid as their primary ID

  if(sttMap.find(trackid)==sttMap.end()) {
    //    std::cout<<"zero stt hit, stt smear fail!"<<std::endl;
    return false;
  }

  TLorentzVector prePos, postPos;
  double dy,dz,dx;
  unsigned int ihit=sttMap[trackid].first;
  int nhit=sttMap[trackid].second;
  //  if(ihit>0) { if(event->SegmentDetectors["Straw"].at(ihit-1).Contrib[0]==trackid) { std::cout<<" wrong"; std::exit(EXIT_FAILURE);}}
  //  if(event->SegmentDetectors["Straw"].size()>(ihit+nhit) ) {if(event->SegmentDetectors["Straw"].at(ihit+nhit).Contrib[0]==trackid) {std::cout<<" wrong"; std::exit(EXIT_FAILURE);}}
  
  TG4HitSegment  h= event->SegmentDetectors["Straw"].at(ihit);
  unsigned int i=(ihit+1);
  if( (event->SegmentDetectors["Straw"].at(ihit+1).Start.T()<h.Start.T() || event->SegmentDetectors["Straw"].at(ihit+1).Stop.T()<h.Stop.T() ) && h.Contrib.size()>1) 
    { h= event->SegmentDetectors["Straw"].at(ihit+1); i=(ihit+2);}
  //  const TG4HitSegment& hseg = ev->SegmentDetectors["Straw"].at(j);
  TLorentzVector mid= (h.Start+h.Stop)*0.5;
  prePos=mid;
  TString name=geo->FindNode(mid.X(),mid.Y(),mid.Z())->GetName();
  if(name.Contains("horizontal")) nYhit++;
  else nXhit++;
  for( ;i<(ihit+nhit);i++){

    h=event->SegmentDetectors["Straw"].at(i);
    postPos= (h.Start+h.Stop)*0.5;
    //    if(h.Contrib[0]!=trackid) { std::cout<<" !!!!!!!!!!!!!!!!!!!!! wrong trackid"<<trackid<<" ihit"<<ihit<<" nhit:"<<nhit<<std::endl; showAll();  std::exit(EXIT_FAILURE);}
    //    if(h.Contrib.size()>1) {  std::cout<<" contribution more than 2 tracks ihit:"<<ihit<<" nhit:"<<nhit<<" i:"<<i<<std::endl;  continue;}
    if(postPos.T()< prePos.T()) { // the time reversed hits are usually bad hits, break it
      //      if((i-ihit)*1.0<=0.5*nhit) std::cout<<"trackid:"<<trackid<<" time reverse, cut! ihit:"<<ihit<<" nhit+ihit:"<<nhit+ihit<<" i:"<<i<<std::endl;
      if((i-ihit)*1.0<0.5*nhit && nhit>15) std::cout<<"thistrack has >15hits but break befor half,mayneed lookinto!!!nhit: "<<nhit<<" i-ihit:"<<i-ihit<<" trackid:"<<trackid<<" iEntry:"<<iEntry<<std::endl;
      break;}
    if(h.EnergyDeposit<250E-6) continue;
    //    postPos= (h.Start+h.Stop)*0.5;
    dx= postPos.X()-prePos.X();
    dy= postPos.Y()-prePos.Y();
    dz= postPos.Z()-prePos.Z();
    name=geo->FindNode(postPos.X(),postPos.Y(),postPos.Z())->GetName();
    if(name.Contains("horizontal")) nYhit++;
    else nXhit++;
    Lyz+= sqrt(dy*dy+dz*dz);
    L+= sqrt(dx*dx+dy*dy+dz*dz);    
    prePos=postPos;
  }

  if(nYhit<4)  { 
    //    std::cout<<"stt smear failed "<<nYhit<<std::endl;
    return false;
  } // fill1Par2tree(-999,-999,-999, trackid, L, nXhit, nYhit, "sttFail") ; return false;}
  //  std::cout<<"nYhit:"<<nYhit<<" nXhit:"<<nXhit<<" postT:"<<postPos.T()<<std::endl;
  /*
  if(h==event->SegmentDetectors["Straw"].end()) return false;
  h++;
  std::cout<<" broken points time:";
  for( ;h != event->SegmentDetectors["Straw"].end(); h++){
    if(h->Contrib[0] ==trackid ) std::cout<<" "<<h->Start.T();
  }
  std::cout<<std::endl;
  */
  TVector3 initP=event->Trajectories[trackid].InitialMomentum.Vect();
  double Pt=sqrt(initP.Y()*initP.Y()+initP.Z()*initP.Z());
  double Px=initP.X();
  double P=initP.Mag();
  double dipAng=atan(Px/Pt);
  double thetaYZ=TMath::ATan2(initP.Y(),initP.Z());
  //  std::cout<<"Lyz:"<<Lyz<<" P:"<<P<<std::endl;
  Lx=sqrt(L*L-Lyz*Lyz);
  L/=1000;
  Lyz/=1000.; //mm-> m
  Lx/=1000.; //mm-> m
  Pt/=1000.;  // Mev-->GeV
  Px/=1000.;
  P/=1000.;
  //  double dPt2Pt=sqrt(pow(sigmas*Pt/0.3/B/Lyz/Lyz*sqrt(720./(nYhit+4)),2)+pow(0.045/B/sqrt(Lyz*x0),2));
  double dPt2Pt=sqrt(pow(sigmas*Pt/0.3/B/L/L*sqrt(720./(nYhit+4)),2)+pow(0.045/B/sqrt(L*x0),2));
  double Pt_smear=Pt*ran->Gaus(1, dPt2Pt);
  while(Pt_smear<=0){
    Pt_smear=Pt*ran->Gaus(1, dPt2Pt);
  }
  double sigma_dipAng=13.6E-3/P*sqrt(L/x0)*(1+0.038*log(L/x0)); // from pdg
  //  double sigma_dipAng2=13.6E-3/P*sqrt(Lx/x0)*(1+0.038*log(Lx/x0)); // from pdg
  //  double sigma_dipAng3=13.6E-3/P*sqrt(Lyz/x0)*(1+0.038*log(Lyz/x0)); // from pdg
  double dipAng_smear=dipAng + ran->Gaus(0,sigma_dipAng);
  //  double dipAng_smear2=dipAng + ran->Gaus(0,sigma_dipAng2);
  //  double dipAng_smear3=dipAng + ran->Gaus(0,sigma_dipAng3);
  double Px_smear=Pt*tan(dipAng_smear);
  double sigma_thetaYZ=13.6E-3/P*sqrt(Lyz/x0)*(1+0.038*log(Lyz/x0)); 
  double thetaYZ_smear=thetaYZ + ran->Gaus(0,sigma_thetaYZ);
  double Py_smear=Pt_smear*sin(thetaYZ_smear);
  double Pz_smear=Pt_smear*cos(thetaYZ_smear);
  double P_smear=sqrt(Pt_smear*Pt_smear+Px_smear*Px_smear);
  //  double Px_smear2=Px + Pt*tan(ran->Gaus(0,sigma_thetaYZ));
  double namecode;
  if(event->Trajectories[trackid].Name=="mu+" || event->Trajectories[trackid].Name=="mu-")
    namecode=0;
  else if(event->Trajectories[trackid].Name=="proton" || event->Trajectories[trackid].Name=="anti_proton")
    namecode=1;
  else if(event->Trajectories[trackid].Name=="pi+" || event->Trajectories[trackid].Name=="pi-")
    namecode=2;
  else if(event->Trajectories[trackid].Name=="e+" || event->Trajectories[trackid].Name=="e-")
    namecode=3;
  else if(event->Trajectories[trackid].Name=="kaon+" || event->Trajectories[trackid].Name=="kaon-")
    namecode=4;
  else 
    namecode=5;
  
  herr_dipAngle_stt->Fill(dipAng_smear-dipAng, namecode);
  herr_pt_stt->Fill((Pt_smear-Pt)/Pt*100,namecode);
  herr_p_stt->Fill((P_smear-P)/P*100,namecode);
  herr_dipAng100_stt->Fill((dipAng_smear-dipAng)/dipAng*100, namecode);
  herr_thetaYZ100_stt->Fill((thetaYZ_smear-thetaYZ)/thetaYZ*100,namecode);
  
  //  double theta=initP.Theta();
  // TVector3 P3smear(Px_smear,Py_smear,Pz_smear);
  ///// constant smearing  2mrad for thetaX and thetaY, 5% for P
  /*
  if(false){
    double thetaX=atan(initP.X()/initP.Z());
    double thetaY=atan(initP.Y()/initP.Z());
    double thetaX_smear=thetaX+ ran->Gaus(0,0.002);
    double thetaY_smear=thetaY+ ran->Gaus(0,0.002);
    double Psmear=initP.Mag()*ran->Gaus(1,0.05); //MeV
    double mod=sqrt(pow(tan(thetaX_smear),2)+pow(tan(thetaY_smear),2)+1);
    if(initP.Z()<0) mod*=-1;
    Px_smear=tan(thetaX_smear)/mod*Psmear/1000.;
    Py_smear=tan(thetaY_smear)/mod*Psmear/1000.;
    Pz_smear=1/mod*Psmear/1000.;
    herr_theta_stt->Fill(thetaX_smear-thetaX, namecode);
  }
  */
  //  std::cout<<"stt smear succeed, L"<<L<<" nXhit:"<<nXhit<<" nYhit:"<<nYhit<<std::endl;
  fill1Par2tree(Px_smear*1000., Py_smear*1000., Pz_smear*1000., trackid, L*1000, nXhit, nYhit, "sttsmear  "); // always use MeV to fill
  
  return true;
}

bool smearChargedPar(int trackid){  
  bool sttSmear=false;
  bool ecalSmear=false;
  sttSmear=smearChargedPar_stt(trackid);

  if(!sttSmear) {
    //    std::cout<<"stt smear fail, will start ecal smear"<<std::endl;
    ecalSmear=smearPar_ecal(trackid); 
  }
  return (sttSmear || ecalSmear);  
}

int findgammaPrimaryId(int trackid){
  int parentId=event->Trajectories[trackid].ParentId;
  if(parentId==-1) return trackid;
  if(event->Trajectories[parentId].Name=="gamma") return findgammaPrimaryId(parentId);
  if(event->Trajectories[parentId].Name=="pi0") return trackid;
  if(event->Trajectories[parentId].Name=="lambda" || event->Trajectories[parentId].Name=="eta" || event->Trajectories[parentId].Name=="sigma0" || event->Trajectories[parentId].Name=="anti_lambda" || event->Trajectories[parentId].Name=="sigma-" || event->Trajectories[parentId].Name=="sigma+") return trackid;
  if(event->Trajectories[parentId].Name=="kaon0L" || event->Trajectories[parentId].Name=="kaon0S") return parentId;
  if(event->Trajectories[parentId].ParentId==-1) { std::cout<<" check this gamma parent(also top):"<<event->Trajectories[parentId].Name<<" gid:"<<trackid<<" iEntry:"<<iEntry<<std::endl;; return parentId;}
  int grandid=event->Trajectories[parentId].ParentId;
  if((event->Trajectories[parentId].Name=="e+" || event->Trajectories[parentId].Name=="e-") && event->Trajectories[grandid].Name=="gamma") return findgammaPrimaryId(grandid);
  std::cout<<"parent:"<<event->Trajectories[parentId].Name<<" grand:"<<event->Trajectories[grandid].Name<<" this need to check trackid:"<<trackid<<" ientry:"<<iEntry<<std::endl;

  return trackid;

}

void smearRemnantGamma(int trackid, int primaryId=-1){
  
  if(primaryId==-1)
    primaryId=findgammaPrimaryId(trackid);

  smearPar_ecal(trackid, primaryId); // this remnant gamma very very not possible to have any stt hits, so just do ecal smearing

}

void smearGamma(int trackid){
  // most time you smear gamma by its daughters one by one, but you don't want to go deep in this loop
  // gamma->gamma->gamma (the third gamma here, we won't smear its daughters but treat it as a remnant gamma
  // pi0->gamma->gamma( do the same thing for the second gamma) 
  //  std::cout<<"start to smear gamma:"<<trackid<<std::endl;
  int parentId=event->Trajectories[trackid].ParentId;
  int grandId=(parentId==-1)?-1:event->Trajectories[parentId].ParentId;
  if(parentId!=-1 && grandId!=-1 && event->Trajectories[parentId].Name=="gamma"){
    if(event->Trajectories[grandId].Name=="gamma") { smearRemnantGamma(trackid); return;}
    if(event->Trajectories[grandId].Name=="pi0") { smearRemnantGamma(trackid, parentId); return;}
    if(event->Trajectories[grandId].Name==("lambda" || "anti_lambda")) { smearRemnantGamma(trackid, parentId); return;}
    if(event->Trajectories[grandId].Name=="kaon0L") { smearRemnantGamma(trackid, grandId); return;}
    if(event->Trajectories[grandId].Name=="kaon0S") { smearRemnantGamma(trackid, grandId); return;}
  }


  Node *no=findNodeFast(trackid);
  no=no->FirstChild;
  if(no==0) { smearRemnantGamma(trackid); return;}
  
  while(no){
    //    std::cout<<"start to smear 1 gamma daughter: "<<no->Traj->TrackId<<"  "<<no->Traj->Name<<std::endl;
    if(no->Traj->Name!="neutron") 
      smearPar(no->Traj->TrackId, no->Traj->Name);
    no=no->RightSibling;
  }
}

void smearDaughters(int trackid){
  //  std::cout<<"start to smearDaughters"<<std::endl;
  Node *no=findNodeFast(trackid);
  no=no->FirstChild;
  while(no){
    if(no->Traj->PDGCode>999999) ; // if pdgcode is too large, it will make dbpdg->GetParticle()->Charge() bad!
    else if(dbpdg->GetParticle(no->Traj->PDGCode)->Charge()!=0)
      smearChargedPar_stt(no->Traj->TrackId); // // for these daughters, should only do STT smear, if only ecal hits, even not possible to recognize them
    else if(no->Traj->Name=="neutron" || no->Traj->Name=="anti_neutron") ;
    else if(no->Traj->Name=="gamma") smearGamma(no->Traj->TrackId);
    else if(no->Traj->Name=="pi0") smearDaughters(no->Traj->TrackId); // don't  use smearPi0 here
    else if(no->Traj->Name=="kaon0L" || no->Traj->Name=="kaon0S" || no->Traj->Name=="lambda" || no->Traj->Name=="eta") { smearDaughters(no->Traj->TrackId);} 
    else { 
      std::cout<<"a 0 charge daughter particle try to smear, will smear daughters!!  Name:"<<no->Traj->Name<<" id:"<<no->Traj->TrackId<<" parentID:"<<trackid<<" ientry"<<iEntry<<std::endl;
      smearDaughters(no->Traj->TrackId);
    }    
    no=no->RightSibling;
  }
}

void smearSTT_or_Daughters(int trackid){ 
  
  bool sttSucceed=smearChargedPar_stt(trackid);
  if(sttSucceed) return;
  else 
    smearDaughters(trackid);

}
void smearPi0_external(int trackid){
  double P=event->Trajectories[trackid].InitialMomentum.P()/1000.;
  double Phi=event->Trajectories[trackid].InitialMomentum.Phi(); // -pi , pi 
  double Theta=event->Trajectories[trackid].InitialMomentum.Theta(); // 0 - pi
  double lowP=0.;
  double highP=4.; // GeV
  int  nPbin=100;
  double lowAng=0;
  double highAng=180;
  int nAngbin=180;
  int iPbin=TMath::CeilNint(P/(highP-lowP)*nPbin);
  int iThetabin=TMath::CeilNint(Theta/(highAng-lowAng)*nAngbin);
  int iPhibin=TMath::CeilNint(abs(Phi)/(highAng-lowAng)*nAngbin);

  if(iPbin>nPbin) std::cout<<"iPbin -------> out of hPi0_mom_recotrue range P:"<<P<<std::endl;
  //  std::cout<<" get random P start"<<"iPbin:"<<iPbin<<std::endl;
  double P_smear= P<4?hPi0_mom_recotrue->ProjectionY("", iPbin,iPbin)->GetRandom():P; // GeV
  if(iPbin>37){ // the stats is too little, have to do in  this way, right now, treat the reco mean same as true mean
    double rms;
    if(iPbin<=39)
      rms=hPi0_mom_recotrue->ProjectionY("",38,39)->GetRMS();
    else if(iPbin<=41)
      rms=hPi0_mom_recotrue->ProjectionY("",40,41)->GetRMS();
    else if(iPbin<=45)
      rms=hPi0_mom_recotrue->ProjectionY("",42,45)->GetRMS();
    else if(iPbin<=50)
      rms=hPi0_mom_recotrue->ProjectionY("",46,50)->GetRMS();
    else if(iPbin<=56)
      rms=hPi0_mom_recotrue->ProjectionY("",51,56)->GetRMS();
    else if(iPbin<=69)
      rms=hPi0_mom_recotrue->ProjectionY("",57,69)->GetRMS();
    else
      rms=hPi0_mom_recotrue->ProjectionY("",70,100)->GetRMS();
    P_smear=P+ran->Gaus(0,rms);
    while(P_smear<=0) {
      P_smear=P+ran->Gaus(0,rms);
    }
  }

  double Theta_smear= hPi0_ang_recotrue->ProjectionY("", iThetabin,iThetabin)->GetRandom();
  //  std::cout<<"get phi"<<std::endl;
  double Phi_smear= hPi0_ang_recotrue->ProjectionY("", iPhibin,iPhibin)->GetRandom();
  //  std::cout<<" get random end ---<"<<std::endl;
  if(Phi<0) Phi_smear*=-1.;

  double Pz_smear=P_smear*cos(Theta_smear);
  double Px_smear=P_smear*sin(Theta_smear)*cos(Phi_smear);
  double Py_smear=P_smear*sin(Theta_smear)*sin(Phi_smear);
  
  fill1Par2tree(Px_smear*1000, Py_smear*1000, Pz_smear*1000 , trackid, -999, -999, -999, "pi0extern ");
  herr_p_pi0->Fill((P_smear-P)/P*100);
  herr_theta_pi0->Fill((Theta_smear-Theta)/Theta*100);
  herr_phi_pi0->Fill((Phi_smear-Phi)/Phi*100);
  
}

void smearPi0(int trackid){
  // only two mode: gamma + gamma / gamma e+ e-, for each mode, no one will miss
  //for the two modes, all the daughters  are always connected, trackid connected.

  // if both gamma convert in ECAL, use external theta true2reco to smear angle, external E true2reco to smear energy
  // if one gamma convert in STT, use charged-track way to reconstruct this gamma, the other gamma convert in ECAL, use 5.7%/sqrt(E) to smear E, use sigmaZ, sigmaX from KLOE ECAL Paper to smear angle
  
  Node *no=findNodeFast(trackid);
  no=no->FirstChild;
  int nChild=0;
  std::string names[3];
  int trackids[3];
  while(no){
    names[nChild]= no->Traj->Name;
    trackids[nChild]=no->Traj->TrackId;
    no=no->RightSibling;
    nChild++;
  }
  if(nChild==2) {
    assert(names[0]=="gamma" && names[1]=="gamma");
    if(sttMap_prim2.find(trackids[0])==sttMap_prim2.end() && sttMap_prim2.find(trackids[1])==sttMap_prim2.end()){ // <--- here need to use primaryId_organized map
      //      std::cout<<"no any hit for pi0 two gammas, use external way to smear pi0"<<std::endl;
      smearPi0_external(trackid);
    }
    else{
      //      std::cout<<" ==pi0==> smear pi0 first gamma: id:"<<trackids[0]<<std::endl;
      smearGamma(trackids[0]);
      //      std::cout<<" ==pi0==> smear pi0 second gamma: id:"<<trackids[1]<<std::endl;
      smearGamma(trackids[1]);
    }
  }
  else if(nChild==3) {
    //    std::cout<<"  ==pi0==> three children for pi0, smear one by one"<<std::endl;
    smearPar(trackids[0], names[0]);
    smearPar(trackids[1], names[1]);
    smearPar(trackids[2], names[2]);
  }
  else std::cout<<" ==pi0==> more than 3 children from pi0, check!!!!"<<std::endl;
}

bool smearN_byEquation(double &psmear, int trackid){  // only for antinumu events with 1 neutron produced 
  if(debug>=2) std::cout<<"--------------------------------------smear neutron by equation -----------------"<<std::endl;
  if(iFillPar<1) return false;
  if(abs(brPdg[0])!=13 && abs(brPdg[0])!=11)  {std::cout<<"^^^^^^^^^^^^^ fail neutron smearing by equation since no lepton reconstructed, iEntry: "<<iEntry<<std::endl; return false;}
  
  const double mpr = dbpdg->GetParticle(2212)->Mass()*1000;
  const double mmu = dbpdg->GetParticle(brPdg[0])->Mass()*1000; // actually it could be electron
  const double mn = dbpdg->GetParticle(2112)->Mass()*1000;

  TLorentzVector p4hadreco(0,0,0,0);
  TLorentzVector p4mureco(brRecoP4[0][0],brRecoP4[0][1],brRecoP4[0][2],brRecoP4[0][3]);
  for(int i=1;i<iFillPar;i++){
    //    if(iChannel==0) std::cout<<"i:"<<i<<" brPdg[i]:"<<brPdg[i]<<std::endl;
    p4hadreco+=TLorentzVector(brRecoP4[i][0],brRecoP4[i][1],brRecoP4[i][2],brRecoP4[i][3]);
  }
  p4mureco.RotateX(-0.101);
  p4hadreco.RotateX(-0.101);

  double en = 0.5*( mmu*mmu + pow(p4hadreco.M(),2) + mpr*mpr - mn*mn - 2*mpr*(p4mureco.E() + p4hadreco.E()) +
  		    2*p4mureco*p4hadreco)/(p4mureco.E() + p4hadreco.E() - p4mureco.Pz() - p4hadreco.Pz() - mpr);
  en = en - p4mureco.E() - p4hadreco.E() + mpr;

  if(en<mn) return false;
  double E=event->Trajectories[trackid].InitialMomentum.E();
  double P=event->Trajectories[trackid].InitialMomentum.P();
  psmear=sqrt(en*en-mn*mn);
  if(abs(brPdg[0])==13)  herr_E_equa_N->Fill((en-E)/E*100, iChannel);
  if(abs(brPdg[0])==13)  herr_p_equa_N->Fill((psmear-P)/P*100, iChannel);
  
    /*
  if(iChannel==0) {
  TLorentzVector p4nu(StdHepP4[0][0],StdHepP4[0][1],StdHepP4[0][2],StdHepP4[0][3]);
  TLorentzVector p4n(event->Trajectories[1].InitialMomentum);
    p4nu.RotateX(-0.101);
    p4nu*=1000;
    p4n.RotateX(-0.101);
    std::cout<<"mpr:"<<mpr<<" mmu:"<<mmu<<" mn:"<<mn<<std::endl;
    std::cout<<"targetpdg:"<<targetpdg<<" neutrinopdg:"<<StdHepPdg[0]<<" brPdg[0]:"<<brPdg[0]<<std::endl;
    std::cout<<"nu px:"<<StdHepP4[0][0]*1000<<" py:"<<StdHepP4[0][1]*1000<<" pz:"<<StdHepP4[0][2]*1000<<" E:"<<StdHepP4[0][3]*1000<<std::endl;
    std::cout<<" true mupx:"<<brTrueP4[0][0]<<" py:"<<brTrueP4[0][1]<<" pz:"<<brTrueP4[0][2]<<" E:"<<brTrueP4[0][3]<<std::endl;
    std::cout<<" true neutron->: "; event->Trajectories[1].InitialMomentum.Print();
    std::cout<<"after transform"<<std::endl;
    std::cout<<"neutrino now:"; p4nu.Print();
    std::cout<<"mu now:"; p4mureco.Print();
    std::cout<<"neutron now:"; p4n.Print();
    p4hadreco.Print();
    std::cout<<"E:"<<E<<" en:"<<en<<" err:"<<(en-E)/E*100<<std::endl;
    if((en-E)/E*100>1) std::cout<<"!!!!!!!!!!!!!!!!!!!!!!! too large check !!!!!!!!!! "<<iEntry<<std::endl;
  }
  */

  return true;
}

bool smearNeutron(int trackid){
  if(event->Trajectories[trackid].ParentId!=-1)  {
    std::cout<<" about to smear a non-primary neutron, stop right now, check why it happens"<<std::endl;
    return false;
  }
  N_organizeHits_contr();

  bool isSTTdetectable=false;
  bool isECALdetectable=false;

  if(N_sttMap_contr_primIndex.find(trackid)!=N_sttMap_contr_primIndex.end()){
    const std::vector<int> &vec= N_sttMap_contr_primIndex[trackid];
    for(unsigned int i=0;i<vec.size()/2;i++){
      for(int j=vec[2*i];j<=vec[2*i+1];j++){
	if(event->SegmentDetectors["Straw"].at(j).EnergyDeposit>250E-6) {isSTTdetectable=true; break;}
      }
      if(isSTTdetectable) break;
    }
    //    if(!isSTTdetectable) std::cout<<"this neutron has straw hits but none of them are detectable, check!!!"<<std::endl;
  }
  if(!isSTTdetectable) {
    if(N_ecalMap_contr_primIndex.find(trackid)!=N_ecalMap_contr_primIndex.end()){
      std::map<int, double> cellId_Evis;
      findEvis_forCell(N_ecalMap_contr_primIndex[trackid], cellId_Evis);
      for(auto cellE: cellId_Evis){
	if(cellE.second > 0.1) { isECALdetectable=true;	  break;}
      }      
    }
  }
  if( !isSTTdetectable && !isECALdetectable) {
    //    std::cout<<"smear neutron fail due to no detectable stt hit and ecal hit"<<std::endl;
    return false;
  }
  
  
  double P_smear;
  double P=event->Trajectories[trackid].InitialMomentum.P();
  bool PequationSmearSucceed=false;
  if(isHtarget) PequationSmearSucceed=smearN_byEquation(P_smear, trackid);

  if(!PequationSmearSucceed || !isHtarget)  {
    double E=event->Trajectories[trackid].InitialMomentum.E();
    //    double m=event->Trajectories[trackid].InitialMomentum.Mag();
    double beta=P/E;
    int iTrueBetaBin=TMath::CeilNint(beta/(1.-0.)*100);
    if(iTrueBetaBin==0) iTrueBetaBin=1; // <--- this solves static neutron decay problem ( neutron with P=0 will decay right at vertex to e-, proton) 
    double beta_smear=isSTTdetectable?hNeutron_beta_recotrue_stt->ProjectionY("",iTrueBetaBin,iTrueBetaBin)->GetRandom(): hNeutron_beta_recotrue_ecal->ProjectionY("",iTrueBetaBin,iTrueBetaBin)->GetRandom();
    //NOTE:  don't use this (when beta_smear>1, P_smear will become nan)  ----->  P_smear=m*beta_smear/sqrt(1-beta_smear*beta_smear);    
    P_smear=E*beta_smear;  // this sometimes will make P_smear > E true, but that's okay
    //    if(beta_smear>1) std::cout<<" beta_smear>1"<<std::endl;
    herr_p_beta_N->Fill((P_smear-P)/P*100);

    //    std::cout<<"beta:"<<beta<<" beta_smear:"<<beta_smear<<" P_smear:"<<P_smear<<std::endl;
  }
  /*
  double Phi=event->Trajectories[trackid].InitialMomentum.Phi(); // -pi , pi
  double Theta=event->Trajectories[trackid].InitialMomentum.Theta(); // 0 - pi
  double Phi_smear=Phi*(1+hNeutron_ang_reso->GetRandom());
  double Theta_smear=Theta*(1+hNeutron_ang_reso->GetRandom());
  double Pz_smear=P_smear*cos(Theta_smear);
  double Px_smear=P_smear*sin(Theta_smear)*cos(Phi_smear);
  double Py_smear=P_smear*sin(Theta_smear)*sin(Phi_smear);
  */
  double thetaX=atan(event->Trajectories[trackid].InitialMomentum.X()/event->Trajectories[trackid].InitialMomentum.Z());
  double thetaY=atan(event->Trajectories[trackid].InitialMomentum.Y()/event->Trajectories[trackid].InitialMomentum.Z());
  double thetaX_smear=thetaX*(1+hNeutron_ang_reso->GetRandom()/sqrt(2.));
  double thetaY_smear=thetaY*(1+hNeutron_ang_reso->GetRandom()/sqrt(2.));
  if(thetaX_smear > TMath::Pi()/2) thetaX_smear = TMath::Pi() - thetaX_smear;
  if(thetaX_smear < -TMath::Pi()/2) thetaX_smear = -TMath::Pi()- thetaX_smear;
  if(thetaY_smear > TMath::Pi()/2) thetaY_smear = TMath::Pi() - thetaY_smear;
  if(thetaY_smear < -TMath::Pi()/2) thetaY_smear = -TMath::Pi()- thetaY_smear;
  double Pz_smear = P_smear/sqrt(1 + tan(thetaX_smear)*tan(thetaX_smear)+tan(thetaY_smear)*tan(thetaY_smear));
  double Px_smear = Pz_smear*tan(thetaX_smear);
  double Py_smear = Pz_smear*tan(thetaY_smear);
  
  if(PequationSmearSucceed) fill1Par2tree(Px_smear, Py_smear, Pz_smear , trackid, -999, -999, -999, "NsmearEqua");
  else fill1Par2tree(Px_smear, Py_smear, Pz_smear , trackid, -999, -999, -999, "NsmearBeta");

  double Theta=event->Trajectories[trackid].InitialMomentum.Theta();
  double theta_Smear=atan(sqrt(Px_smear*Px_smear+Py_smear*Py_smear)/Pz_smear);
  herr_theta_N->Fill((theta_Smear-Theta)/Theta);
  
  //  herr_phi_N->Fill((Phi_smear-Phi)/Phi*100);
  
  return true;
}

void organizeHits(){
  sttMap.clear();
  //  ecalMap.clear();
  int pretrackid;
  int posttrackid;
  int nhit;
  int istart;

  if(event->SegmentDetectors["Straw"].size()>0){    
    pretrackid=event->SegmentDetectors["Straw"].begin()->Contrib[0];
    nhit=1;
    istart=0;
    posttrackid=pretrackid;
    for(unsigned int i=1; i<event->SegmentDetectors["Straw"].size(); i++){
      posttrackid=event->SegmentDetectors["Straw"].at(i).Contrib[0];
      //      std::cout<<"posttrackid:"<<posttrackid<<" nhit:"<<nhit<<std::endl;
      if(posttrackid==pretrackid) { nhit++;continue;}
      if(sttMap.find(pretrackid) ==sttMap.end()) 
	sttMap[pretrackid]= std::make_pair(istart,nhit);
      if(sttMap.find(posttrackid) ==sttMap.end()){
	nhit=1;
	istart=i;
      }    
      pretrackid=posttrackid;
    }
    if(sttMap.find(posttrackid) ==sttMap.end())
      sttMap[posttrackid]= std::make_pair(istart,nhit);
  }
  for(auto it = sttMap.begin(); it != sttMap.end(); ) {
    if(it->second.second <4)
      it = sttMap.erase(it);
    else
      ++it;
  }

  /*
  if(event->SegmentDetectors["ECAL"].size()>0){
    pretrackid=event->SegmentDetectors["ECAL"].begin()->Contrib[0];
    nhit=1;
    istart=0;
    posttrackid=pretrackid;
    for(unsigned long i=1; i<event->SegmentDetectors["ECAL"].size(); i++){
      posttrackid=event->SegmentDetectors["ECAL"].at(i).Contrib[0];
      //      std::cout<<"posttrackid:"<<posttrackid<<" nhit:"<<nhit<<std::endl;
      if(posttrackid==pretrackid) { nhit++;continue;}
      if(ecalMap.find(pretrackid) ==ecalMap.end())
        ecalMap[pretrackid]= std::make_pair(istart,nhit);
      if(ecalMap.find(posttrackid) ==ecalMap.end()){
        nhit=1;
        istart=i;
      }
      pretrackid=posttrackid;
    }
    if(ecalMap.find(posttrackid) ==ecalMap.end())
      ecalMap[posttrackid]= std::make_pair(istart,nhit);
  }
  */
}
/*
void organizeHits_prim(){
  sttMap_prim.clear();
  ecalMap_prim.clear();
  int pretrackid;
  int posttrackid;
  int nhit;
  int istart;
  if(event->SegmentDetectors["Straw"].size()>0){
    pretrackid=event->SegmentDetectors["Straw"].begin()->PrimaryId;
    nhit=1;
    istart=0;
    posttrackid=pretrackid;
    for(unsigned int i=1; i<event->SegmentDetectors["Straw"].size(); i++){
      posttrackid=event->SegmentDetectors["Straw"].at(i).PrimaryId;
      //      std::cout<<"posttrackid:"<<posttrackid<<" nhit:"<<nhit<<std::endl;
      if(posttrackid==pretrackid) { nhit++;continue;}
      if(sttMap_prim.find(pretrackid) ==sttMap_prim.end())
        sttMap_prim[pretrackid]= std::make_pair(istart,nhit);
      if(sttMap_prim.find(posttrackid) ==sttMap_prim.end()){
        nhit=1;
        istart=i;
      }
      pretrackid=posttrackid;
    }
    if(sttMap_prim.find(posttrackid) ==sttMap_prim.end())
      sttMap_prim[posttrackid]= std::make_pair(istart,nhit);
  }

  if(event->SegmentDetectors["ECAL"].size()>0){
    pretrackid=event->SegmentDetectors["ECAL"].begin()->PrimaryId;
    nhit=1;
    istart=0;
    posttrackid=pretrackid;
    for(unsigned long i=1; i<event->SegmentDetectors["ECAL"].size(); i++){
      posttrackid=event->SegmentDetectors["ECAL"].at(i).PrimaryId;
      //      std::cout<<"posttrackid:"<<posttrackid<<" nhit:"<<nhit<<std::endl;
      if(posttrackid==pretrackid) { nhit++;continue;}
      if(ecalMap_prim.find(pretrackid) ==ecalMap_prim.end())
        ecalMap_prim[pretrackid]= std::make_pair(istart,nhit);
      if(ecalMap_prim.find(posttrackid) ==ecalMap_prim.end()){
        nhit=1;
        istart=i;
      }
      pretrackid=posttrackid;
    }
    if(ecalMap_prim.find(posttrackid) ==ecalMap_prim.end())
      ecalMap_prim[posttrackid]= std::make_pair(istart,nhit);
  } 
}
*/
void organizeHits_prim2(){
  sttMap_prim2.clear();
  ecalMap_prim2.clear();
  int trackid;
  for(unsigned int i=0; i<event->SegmentDetectors["Straw"].size(); i++){
    trackid=event->SegmentDetectors["Straw"].at(i).PrimaryId;
    if(sttMap_prim2.find(trackid)==sttMap_prim2.end()) { sttMap_prim2[trackid].push_back(i); sttMap_prim2[trackid].push_back(i); continue;}
    if(i-sttMap_prim2[trackid].back()==1)
      sttMap_prim2[trackid].back()=i;
    else
      {sttMap_prim2[trackid].push_back(i); sttMap_prim2[trackid].push_back(i);}    
  } // for
  
  for(unsigned int i=0; i<event->SegmentDetectors["ECAL"].size(); i++){
    trackid=event->SegmentDetectors["ECAL"].at(i).PrimaryId;
    if(ecalMap_prim2.find(trackid)==ecalMap_prim2.end()) { ecalMap_prim2[trackid].push_back(i); ecalMap_prim2[trackid].push_back(i);continue;}
    if(i-ecalMap_prim2[trackid].back()==1)
      ecalMap_prim2[trackid].back()=i;
    else
      {ecalMap_prim2[trackid].push_back(i); ecalMap_prim2[trackid].push_back(i);}
  } // for
}

bool isTopNeutron(int trackid, int& topId){
  TG4Trajectory trk=event->Trajectories[trackid];

  while(trk.ParentId!=-1){
    trk=event->Trajectories[trk.ParentId];
  }
  topId=trk.TrackId;
  return (trk.Name=="neutron");
}

void N_organizeHits_contr(){
  std::map<int, std::vector<int> >  sttMap_temp;
  std::map<int, std::vector<int> >  ecalMap_temp;
  N_sttMap_contr_primIndex.clear();
  N_ecalMap_contr_primIndex.clear();
  int trackid;
  for(unsigned int i=0; i<event->SegmentDetectors["Straw"].size(); i++){
    trackid=event->SegmentDetectors["Straw"].at(i).Contrib[0];
    if(sttMap_temp.find(trackid)==sttMap_temp.end()) { sttMap_temp[trackid].push_back(i); sttMap_temp[trackid].push_back(i); continue;}
    if(i-sttMap_temp[trackid].back()==1)
      sttMap_temp[trackid].back()=i;
    else
      {sttMap_temp[trackid].push_back(i); sttMap_temp[trackid].push_back(i);}
  } // for

  for(unsigned int i=0; i<event->SegmentDetectors["ECAL"].size(); i++){
    trackid=event->SegmentDetectors["ECAL"].at(i).Contrib[0];
    if(ecalMap_temp.find(trackid)==ecalMap_temp.end()) { ecalMap_temp[trackid].push_back(i); ecalMap_temp[trackid].push_back(i);continue;}
    if(i-ecalMap_temp[trackid].back()==1)
      ecalMap_temp[trackid].back()=i;
    else
      {ecalMap_temp[trackid].push_back(i); ecalMap_temp[trackid].push_back(i);}
  } // for
  int topid;
  if(!sttMap_temp.empty()){
    for(auto it = sttMap_temp.begin(); it != sttMap_temp.end(); it++ ) {
      if(!isTopNeutron(it->first,topid)) continue;
      else
	N_sttMap_contr_primIndex[topid].insert(N_sttMap_contr_primIndex[topid].end(),it->second.begin(),it->second.end());
      
    }
  }
  if(!ecalMap_temp.empty()){
    for(auto it = ecalMap_temp.begin(); it != ecalMap_temp.end(); it++) {
      if(!isTopNeutron(it->first,topid)) continue;
      else
	N_ecalMap_contr_primIndex[topid].insert(N_ecalMap_contr_primIndex[topid].end(),it->second.begin(),it->second.end());      
    }
  }
  /*
  std::cout<<"-------N_sttMap_contr_primIndex map show:"<<std::endl;
  for(auto pp: N_sttMap_contr_primIndex){
    std::cout<<pp.first<<" ->";
    assert(pp.second.size()%2==0);
    for(auto ii:pp.second)
      std::cout<<" "<<ii;
    std::cout<<std::endl;
  }
  std::cout<<"-------N_ecalMap_contr_primIndex  map show:"<<std::endl;
  for(auto pp: N_ecalMap_contr_primIndex){
    std::cout<<pp.first<<" ->";
    assert(pp.second.size()%2==0);
    for(auto ii:pp.second)
      std::cout<<" "<<ii;
    std::cout<<std::endl;
  }
  */
}

void showHitMap(){
  std::cout<<"-------stt map show:"<<std::endl;
  for(auto pp:sttMap){
    std::cout<<pp.first<<" "<<pp.second.first<<" "<<pp.second.second<<std::endl;
  }
  /*
  std::cout<<"-------ecal map show:"<<std::endl;
  for(auto pp:ecalMap){
    std::cout<<pp.first<<" "<<pp.second.first<<" "<<pp.second.second<<std::endl;
  }
  */
  std::cout<<"-------stt prim2 map show:"<<std::endl;
  for(auto pp:sttMap_prim2){
    std::cout<<pp.first<<" ->";
    assert(pp.second.size()%2==0);
    for(auto ii:pp.second)
      std::cout<<" "<<ii;
    std::cout<<std::endl;
  }
  std::cout<<"-------ecal prim2 map show:"<<std::endl;
  for(auto pp:ecalMap_prim2){
    std::cout<<pp.first<<" ->";
    assert(pp.second.size()%2==0);
    for(auto ii:pp.second)
      std::cout<<" "<<ii;
    std::cout<<std::endl;
  }
  
}
void smearPar(int trackid, std::string name){

  if(name=="mu-" || name=="mu+" || name=="proton" || name=="pi+" || name=="pi-") smearChargedPar(trackid);
  else if(name=="e-" || name=="e+") smearChargedPar(trackid);
  else if(name=="kaon+" || name=="kaon-") smearChargedPar(trackid);
  else if(name=="pi0") smearPi0(trackid);
  else if(name=="neutron") smearNeutron(trackid);
  else if(name=="gamma") smearGamma(trackid);  
  else if(name=="kaon0" || name=="anti_kaon0" || name=="kaon0L") smearDaughters(trackid);
  else if(name=="lambda" || name=="anti_lambda")  smearDaughters(trackid);   // lambda travel very short distance (5cm) then decay to pi- and proton
  else if(name=="sigma+" || name=="sigma-") smearSTT_or_Daughters(trackid); // sigma+/- can create medium track then decay 
  else if(name=="sigma0") smearDaughters(trackid);  // sigma0 decay right at vertex to 1 lambda and 1 gamma
  else if(name=="eta") smearDaughters(trackid); //eta decay right there to 3 pi0s or other modes
  else if(name=="anti_proton") smearChargedPar(trackid);
  else if(name=="nu_e"|| name=="anti_nu_e" || name=="anti_nu_mu" || name=="nu_mu" || name=="C12" || name=="O16" || name=="Ar40" || name=="deuteron") return;
  else if(name=="Fe56" || name=="Al27" || name=="Pb207") return;
  else if(name=="anti_neutron") smearNeutron(trackid);
  else {
    std::cout<<"--->####################################### unknown par:"<<name<<"  iEntry:"<<iEntry<<" trackid:"<<trackid<<std::endl;
    std::cout<<"--->####################################### unknown par:"<<name<<std::endl;
    std::cout<<"--->####################################### unknown par:"<<name<<std::endl;
  }
  //  else std::cout<<"--->unknown par:"<<name<<std::endl;
}

bool checkPrimNeutron(bool &NeutronAtLast, int &Neutron_trackId){

  Neutron_trackId=-1;
  NeutronAtLast=false;
  int nPrim= event->Primaries.begin()->Particles.size();
  int i=0;
  for (;i<nPrim;i++){    
    if(event->Primaries.begin()->Particles[i].Name=="neutron")
      { Neutron_trackId=i; NeutronAtLast=(i==(nPrim-1)); return true;}
  }
  return false;

}

void smearEvent(){
  cleanBranch();
  //  treefile.open(Form("treemap%d.txt",i));
  //  std::cout<<" ############################################################## new event ####################################  "<<i<<std::endl;
  if(debug>=3)  showAll();
  organizeHits();
  //  organizeHits_prim();
  organizeHits_prim2();
  //  showHitMap();
  makeTree();
  if(debug>=2) drawTree();
  //  treefile.close();
  iFillPar=0;
  int nPrim= event->Primaries.begin()->Particles.size();
  bool havePrimNeutron=false;
  bool NeutronAtLast=false;
  int neutron_trackid;
  if(isHtarget) havePrimNeutron=checkPrimNeutron(NeutronAtLast, neutron_trackid);
  if(isHtarget && havePrimNeutron && (!NeutronAtLast)) {
    for (int i=0;i<nPrim;i++){      
      if(event->Primaries.begin()->Particles[i].TrackId==neutron_trackid) continue;
      if(debug>=2) std::cout<<" +++++++++++++++++++++++++++++++++++++++++++++++++++++smearPar: "<<event->Primaries.begin()->Particles[i].Name<<" id: "<<event->Primaries.begin()->Particles[i].TrackId<<std::endl;
      smearPar(event->Primaries.begin()->Particles[i].TrackId , event->Primaries.begin()->Particles[i].Name);
    }
    if(debug>=2)   std::cout<<"smear neutron at the very end since its a numubar hydrogen target event"<<std::endl;
    smearNeutron(neutron_trackid);
  }
  else {
    for (int i=0;i<nPrim;i++){
      if(debug>=2) std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++smearPar: "<<event->Primaries.begin()->Particles[i].Name<<" id: "<<event->Primaries.begin()->Particles[i].TrackId<<std::endl;
      smearPar(event->Primaries.begin()->Particles[i].TrackId , event->Primaries.begin()->Particles[i].Name);
    }
  }
  //  std::cout<<"iFillPar:"<<iFillPar<<std::endl;
  if(iFillPar==0) return;
  brNPar=iFillPar;
  brNPrim=nPrim;
  brRecoNuP4[0]=recoNuPx;
  brRecoNuP4[1]=recoNuPy;
  brRecoNuP4[2]=recoNuPz;
  double recoNuE=sqrt(recoNuPx*recoNuPx+recoNuPy*recoNuPy+recoNuPz*recoNuPz);
  brRecoNuP4[3]=recoNuE;
  tree->Fill();

  herr_nu_E->Fill((recoNuE-trueNuE)/trueNuE);
  //  std::cout<<"recoNuE:"<<recoNuE<<" trueNuE:"<<trueNuE<<std::endl;
}

std::vector<std::string> makefilelist(std::string st,int Nfilelist=0){
  // if Nfilelist is 0, then the default means input all the lines/files  i have.
  std::vector<std::string> files;
  files.clear();

  std::ifstream filelist(st.c_str());
  //count how many lines in the filelistlist
  int num=std::count(std::istreambuf_iterator<char>(filelist), std::istreambuf_iterator<char>(), '\n');

  if(Nfilelist==0)
    Nfilelist=num;
  else if(Nfilelist>num)
    std::cout<<"we don't have that many filelist, please change the  number."<<std::endl;
  //the following close and reopen is important, cause if you didn't close, right now the pointer is at the end of the filelistlist
  filelist.close();
  filelist.open(st.c_str());

  std::string onefile;
  for(int i=0;i<Nfilelist;i++){
    getline(filelist,onefile);
    files.push_back(onefile);
  }
  filelist.close();
  return files;
}//////////////////////////////////////////////////////////////////


int main(int argc, char *argv[]){
  if(argc<4) { std::cout<<"at least three arguments needed: smeartype(0: regular, 1:sttfv) , input, output, [debug1/2/3],[startentry],[nentry]"<<std::endl;return 0;}
  debug=-1;
  bool sttFVsmear=std::atoi(argv[1]);
  std::cout<<"sttFVsmear:"<<sttFVsmear<<std::endl;
  int testStartEntry=-1;
  int testNEntry=-1;
  if(argc>=5){
    if(std::strcmp(argv[4],"debug1")==0) debug=1;
    else if(std::strcmp(argv[4],"debug2")==0) debug=2;
    else if(std::strcmp(argv[4],"debug3")==0) debug=3;
  }
  if(argc>=6) { testStartEntry=std::atoi(argv[5]);}
  if(argc==7) { testNEntry=std::atoi(argv[6]);}
  std::cout<<" %%%%%%%%%%%%%%%%%%%%%%%% debug level ########################   :  "<<debug<<std::endl;

  //  outf=new TFile("outf.root","recreate");
  TFile *fpi0=new TFile("data/Histograms_1pi0_complete.root");
  TFile *fneutron_beta=new TFile("data/RecoVsTrue_Beta_RealCal_pdg2112_20190315_192325_ThrStt2.500000e-07.root");
  hPi0_mom_recotrue=(TH2*)fpi0->Get("h_mom_recotrue");
  hPi0_ang_recotrue=(TH2*)fpi0->Get("h_arctg_recotrue");
  hNeutron_beta_recotrue_stt=(TH2*)fneutron_beta->Get("Beta STT");
  hNeutron_beta_recotrue_ecal=(TH2*)fneutron_beta->Get("Beta calorimeter");
  /*
    TFile *fneutron_ang=new TFile("data/Neutron_H_Angle_qeok_smearingok_All.root");
    TH2 *hNeutron_ang_reso_2D=(TH2*)fneutron_ang->Get("isto_res_stt");
    hNeutron_ang_reso_2D->Add((TH2*)fneutron_ang->Get("isto_res_calbarrel"));
    hNeutron_ang_reso=hNeutron_ang_reso_2D->ProjectionY("hNeutron_ang_reso",1,30);
  */
  //  TFile *fneutron_ang=new TFile("data/plotres.root");
  //  hNeutron_ang_reso=(TH1*)fneutron_ang->Get("STT Resolution");
  //  hNeutron_ang_reso->Add((TH1*)fneutron_ang->Get("Calorimeter Resolution"));
  TFile *fneutron_ang=new TFile("data/plotres.root");
  hNeutron_ang_reso=(TH1*)fneutron_ang->Get("hRes1DTot");
  std::cout<<"hNeutron_ang_reso rms"<<hNeutron_ang_reso->GetRMS()<<std::endl;
  
  hNeutron_ang_reso->Smooth();
  hPi0_mom_recotrue->Smooth();
  hPi0_mom_recotrue->Smooth();
  hPi0_ang_recotrue->Smooth();
  hPi0_ang_recotrue->Smooth();
  hNeutron_beta_recotrue_stt->Smooth();
  hNeutron_beta_recotrue_ecal->Smooth();
  
  outTreeF=new TFile(argv[3],"recreate");

  TBranch * brEvtCode;
  TObjString* EvtCode = 0;
  TBranch * brStdHepPdg=0;
  double  StdHepP4[kNPmax][4];
  int  StdHepPdg[kNPmax];
  TTree *rootrackerTree;

  tree = new TTree("edep_smeared_tree"," reconstructed on edep sim");
  tree->Branch("IEntry", &brIEntry, "IEntry/I");
  tree->Branch("NPar", &brNPar, "NPar/I");
  tree->Branch("NPrim", &brNPrim, "NPrim/I");
  tree->Branch("RecoP4",    brRecoP4,        "RecoP4[NPar][4]/D"); 
  tree->Branch("TrueP4",    brTrueP4,        "TrueP4[NPar][4]/D");
  tree->Branch("RecoNuP4",    brRecoNuP4,        "RecoNuP4[4]/D");
  tree->Branch("TrueNuP4",    brTrueNuP4,        "TrueNuP4[4]/D");
  tree->Branch("Pdg",    brPdg,        "Pdg[NPar]/I");
  tree->Branch("TrackId",    brTrackId,        "TrackId[NPar]/I");
  tree->Branch("ParentId",    brParentId,        "ParentId[NPar]/I");
  tree->Branch("TopParentId",    brTopParentId,        "TopParentId[NPar]/I");
  tree->Branch("Len",    brLen,        "Len[NPar]/D");
  tree->Branch("NXhit",    brNXhit,        "NXhit[NPar]/I");
  tree->Branch("NYhit",    brNYhit,        "NYhit[NPar]/I");
  tree->Branch("Info",    brInfo,        "Info[NPar][10]/C");

  herr_dipAngle_stt=new TH2F("herr_dipAngle_stt","",200,-0.05,0.05,10,0,10); // rad
  herr_dipAng100_stt=new TH2F("herr_dipAng100_stt","",100,-30,30,10,0,10); // percent
  herr_thetaYZ100_stt=new TH2F("herr_thetaYZ100_stt","",100,-30,30,10,0,10); // percent
  herr_theta_stt=new TH2F("herr_theta_stt","",200,-0.05,0.05,10,0,10);
  herr_pt_stt=new TH2F("herr_pt_stt","",100, -30,30,10,0,10); // percent
  herr_p_stt=new TH2F("herr_p_stt","",100, -30,30,10,0,10); //percent
  herr_E_ecal=new TH1F("herr_E_ecal","",100,-30,30); 
  herr_theta_ecal=new TH1F("herr_theta_ecal","",100,-30,30);
  herr_phi_ecal=new TH1F("herr_phi_ecal","",100,-30,30);
  herr_theta_N=new TH1F("herr_theta_N","",100,-30,30);
  herr_phi_N=new TH1F("herr_phi_N","",100,-30,30);
  herr_E_equa_N=new TH2F("herr_E_equa_N","",100,-50,50, 5, 0,5);
  herr_p_equa_N=new TH2F("herr_p_equa_N","",100,-50,50, 5, 0,5);
  herr_p_beta_N=new TH1F("herr_p_beta_N","",100,-50,50);
  herr_p_pi0=new TH1F("herr_p_pi0","",100,-50,50);
  herr_theta_pi0=new TH1F("herr_theta_pi0","",100,-50,50);
  herr_phi_pi0=new TH1F("herr_phi_pi0","",100,-50,50);
  herr_nu_E=new TH1F("herr_nu_E","",100,-1,1);

  ran=new TRandom3(0); // 0 will always give different result when you recreate it
  //  ran->SetSeed(3722147861);
  std::cout<<"seed:"<<ran->GetSeed()<<std::endl;
  dbpdg = new TDatabasePDG();
  gSystem->Load("libGeom");  
  
  std::cout<<"****************** input file:"<<argv[2]<<std::endl;
  gFile=new TFile(argv[2]);
  geo = (TGeoManager*) gFile->Get("EDepSimGeometry");
  gEDepSimTree = (TTree*) gFile->Get("EDepSimEvents");
  gEDepSimTree->SetBranchAddress("Event",&event);    
  rootrackerTree=(TTree*) gFile->Get("DetSimPassThru/gRooTracker");
  brStdHepPdg    = rootrackerTree -> GetBranch ("StdHepPdg");
  TBranch * brStdHepP4        = rootrackerTree -> GetBranch ("StdHepP4");
  brEvtCode      = rootrackerTree -> GetBranch ("EvtCode");
  brStdHepPdg -> SetAddress (StdHepPdg);
  brEvtCode   -> SetAddress (&EvtCode);
  brStdHepP4        -> SetAddress (  StdHepP4 );

  outTreeF->cd();

  int nEntry=gEDepSimTree->GetEntries();
  int rootrackerEntry=rootrackerTree->GetEntries();
  if(nEntry!=rootrackerEntry) {std::cout<<"----->----->not same entries"<<std::endl; return 1;}
  std::cout<<"------------------------------------------------how many entries in this file:"<<nEntry<<std::endl;
  int startEntry=(testStartEntry!=-1)?testStartEntry:0;
  int endplusEntry=(testNEntry!=-1)?(startEntry+testNEntry):nEntry;
  for(int i=startEntry;i<endplusEntry;i++){
    if(i%100==0) std::cout<<"=====> ientry:"<<i<<std::endl;
    if(debug>=1) std::cout<<" ############################################################## new event ####################################  "<<i<<std::endl;
    iEntry=i;
    gEDepSimTree->GetEntry(i);
    if(sttFVsmear) {
      TLorentzVector vtx=event->Primaries.begin()->GetPosition();
      if(!inFV(vtx.X(),vtx.Y(),vtx.Z())) continue;
    }
    isHtarget=false;
    brIEntry=i;
    rootrackerTree->GetEntry(i);
    const char * code=EvtCode->String().Data();
    if(std::strstr(code,"DIS"))
      iChannel=2;
    else if(std::strstr(code,"RES"))
      iChannel=1;
    else if(std::strstr(code,"QES"))
      iChannel=0;
    else if(std::strstr(code,"COH"))
      iChannel=3;
    else
      { std::cout<<"--code --->"<<code<<std::endl;std::exit(EXIT_FAILURE);}

    targetpdg=StdHepPdg[1];
    trueNuE = StdHepP4[0][3]*1000.;
    brTrueNuP4[0]=StdHepP4[0][0]*1000.;
    brTrueNuP4[1]=StdHepP4[0][1]*1000.;
    brTrueNuP4[2]=StdHepP4[0][2]*1000.;
    brTrueNuP4[3]=StdHepP4[0][3]*1000.;

    if (StdHepPdg[1]==2212)  isHtarget=true;
    //    if(StdHepPdg[0]!=14 && StdHepPdg[0]!=-14 && StdHepPdg[1]==2212) std::cout<<" electron neutrino + htarget"<<" StdHepPdg[1]:"<<StdHepPdg[1]<<std::endl;
    //    std::cout<<"nupx:"<<StdHepP4[0][0]<<" nupy:"<<StdHepP4[0][1]<<" nupz:"<<StdHepP4[0][2]<<" nue:"<<StdHepP4[0][3]<<std::endl;
    smearEvent();
    
  }    
  std::cout<<"close files"<<std::endl;
  
  tree->Write();
  herr_dipAngle_stt->Write();
  herr_dipAng100_stt->Write();
  herr_thetaYZ100_stt->Write();
  herr_pt_stt->Write();
  herr_p_stt->Write();
  herr_theta_stt->Write();
  herr_E_ecal->Write();
  herr_theta_ecal->Write();
  herr_phi_ecal->Write();
  herr_theta_N->Write();
  herr_phi_N->Write();
  herr_E_equa_N->Write();
  herr_p_equa_N->Write();
  herr_p_beta_N->Write();
  herr_p_pi0->Write();
  herr_theta_pi0->Write();
  herr_phi_pi0->Write();
  herr_nu_E->Write();
  //  outTreeF->Write();
  outTreeF->Close();

}
/*
 --------------------------------------------------------------------------------------------------------------------------
 | Idx | Ist |    PDG     | Rescat |   Mother  |  Daughter |   Px   |   Py   |   Pz   |   E    |   x    |   y    |    z   |
 |     |     |            |        |           |           |(GeV/c) |(GeV/c) |(GeV/c) | (GeV)  |  (fm)  |  (fm)  |   (fm) |

 |   0 |   0 |        -14 | -280947920 |  -1 |  -1 |   4 |   4 |  0.003 | -0.289 |  2.791 |  2.806 |  1.184 |  1.936 |  1.051 |
 |   1 |   0 | 1000060120 |  32766 |  -1 |  -1 |   2 |   3 |  0.000 |  0.000 |  0.000 | 11.175 |  0.000 |  0.000 |  0.000 |
 |   2 |  11 |       2212 | -649084800 |   1 |  -1 |   5 |   5 |  0.113 | -0.021 | -0.154 |  0.921 |  1.184 |  1.936 |  1.051 |
 |   3 |   2 | 1000050110 |  32666 |   1 |  -1 |  16 |  16 | -0.113 |  0.021 |  0.154 | 10.254 |  0.000 |  0.000 |  0.000 |
 |   4 |   1 |        -13 |      4 |   0 |  -1 |  -1 |  -1 | -0.089 | -0.177 |  0.000 |  0.225 |  1.184 |  1.936 |  1.051 |

*/
