

double createBeamTime(TRandom *ran){
  int nbatch=6;
  int nbunch=86; //
  int nEmpty1Bunch=2; //
  //  double frequency=53.1E6;  // per second
  double time1Bunch=19; // ns
  //  double time1BunchRMS=1; //ns
  int ibatch=ran->Integer(nbatch);
  int ibunch=ran->Integer(nbunch-nEmpty1Bunch);
  double t=(ibatch*nbunch+ibunch)*time1Bunch + ran->Gaus(0,1.5);
  //  hBeamT->Fill(t);
  return t; // ns
}
