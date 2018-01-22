#define Mixing_cxx
#include "Mixing.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <valarray>
void Mixing:: Eccentricity()
{
//   In a ROOT session, you can do:
//      root> .L Mixing.C
//      root> Mixing t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch


   if (fChain == 0) return;

     Long64_t nentries = fChain->GetEntriesFast();
     Long64_t nbytes = 0, nb = 0;

     TH1D *ecc1DHisto = new TH1D("ecc1DHisto", "Eccentricity", 100,0,1.2);
     //     TH1D *lin_ecc1DHisto = new TH1D("lin_ecc1DHisto", "Linear Eccentricity", 300,0,6);
     TH2D *m2DHisto = new TH2D("m2DHisto", "m20 vs m02", 40, 0 , 20, 40, 0 , 20);
     double Ecc;
     //     double LinEcc;

     for (Long64_t jentry=0; jentry<nentries;jentry++) {
       Long64_t ientry = LoadTree(jentry);
       if (ientry < 0) break;

	 nb = fChain->GetEntry(jentry);   nbytes += nb;	
	 for (int n=0; n< ncluster; n++){  //initially forgot to loop over ncluster in each event
	   if ((cluster_e[n] > 6) && (cluster_ncell[n]>1)){	   
	   if ( (cluster_lambda_square[n][0] > 0) && (cluster_lambda_square[n][0]>0) ){ //how to make sense of 0's?
	     Ecc = sqrt (1 - ((pow(cluster_lambda_square[n][0],2)) / (pow(cluster_lambda_square[n][0],2))) ); //m02 is major
	     ecc1DHisto->Fill(Ecc); // not finding eccentricity if either axis yields zero

	     //  LinEcc = Ecc*(cluster_lambda_square[n][0],2);
	     //  lin_ecc1DHisto->Fill(LinEcc);

	   }

	   m2DHisto->Fill(cluster_lambda_square[n][0], cluster_lambda_square[n][0]);

	 }

	 // if (Cut(ientry) < 0) continue;
       }
     }

     
     TCanvas *e1D = new TCanvas("e1D", "Eccentricity", 800,800);
     ecc1DHisto->Draw();


     TCanvas *c2D = new TCanvas("c2D","MajorMinor",800,800);
     m2DHisto->GetXaxis()->SetTitle("m02");
     m2DHisto->GetYaxis()->SetTitle("m20");
     gStyle->SetOptStat("");
     c2D->SetLogz();
     m2DHisto->Draw("COLZ");

     // TCanvas *l1D = new TCanvas("l1D", "LinearEccentricity", 800,800);
     // lin_ecc1DHisto->Draw();

}


void Mixing::Correlation()
{
  if (fChain == 0) return;

     Long64_t nentries = fChain->GetEntriesFast();
     Long64_t nbytes = 0, nb = 0;

     TH2D *Corr = new TH2D("Correlation", "Inclusive gamma - hadron Correlation", 60,-M_PI/2,3*M_PI/2, 34, -1.7, 1.7); 
                                                                                  //range -pi/2 to 3pi/2
     Corr->Sumw2();

     double DeltaPhi;
     double Cluster_rads;
     double Track_rads;

     double DeltaEta;
     double Cluster_eta;
     double Track_eta;
     
     for (Long64_t jentry=0; jentry<nentries;jentry++) {
       Long64_t ientry = LoadTree(jentry);
       if (ientry < 0) break;
       cout<<jentry<<endl;
       fChain->GetEntry(jentry);  
       for (int n=0; n < ncluster; n++){
	 if ((cluster_e[n] > 6) && (cluster_ncell[n]>1)){
	     
	   for (int m=0; m < ntrack; m++){
	     if(track_pt[m] > 1){ //Cut on 1GeV Tracks     
	   
	       DeltaPhi = (cluster_phi[n] - track_phi[m]);
	       if (DeltaPhi < -M_PI/2){DeltaPhi= DeltaPhi + 2*M_PI;}  //if less then -pi/2 add 2pi
	       if (DeltaPhi > 3*M_PI/2){DeltaPhi = DeltaPhi -2*M_PI;}
	       DeltaEta = (cluster_eta[n] - track_eta[m]);
	       Corr->Fill(DeltaPhi,DeltaEta);
	     }
	   }
	 }
       }
     }
     
     TCanvas *Corr2D = new TCanvas("Corr2D", "Inc_g-h_Correlation", 800,800);
     Corr->SetMinimum(0.); //set Y-Axis minimum to 0 when drawing 
     Corr->Draw("SURF1");
     TFile *MyFile = new TFile("g_h_Correlation.root","RECREATE");
     Corr->Write();
     MyFile->Print(); //write to file so I don't have to run this again

}


void Mixing::H_HCorrelation()
{
  if (fChain == 0) return;

     Long64_t nentries = fChain->GetEntriesFast();
     Long64_t nbytes = 0, nb = 0;

     TH2D *Corr = new TH2D("Correlation", "hadron - hadron Correlation", 60,-M_PI/2,3*M_PI/2, 34, -1.7, 1.7); //range -pi/2 to 3pi/2
     Corr->Sumw2();

     double DeltaPhi;
     double DeltaEta;
     
     for (Long64_t jentry=0; jentry<nentries;jentry++) { //normally go until 'nentries', shortening to text mixing
       Long64_t ientry = LoadTree(jentry);
       if (ientry < 0) break;
       cout<<jentry<<endl;
       fChain->GetEntry(jentry);  
       for (int n=0; n < ntrack; n++){
	 if (track_pt[n] > 1){
	     
	   for (int m=0; m < ntrack; m++){
	     if((track_pt[m] > 1) && (n != m)){ //Cut on 1GeV Tracks     
	   
	       DeltaPhi = (track_phi[n] - track_phi[m]);
	       if (DeltaPhi < -M_PI/2){DeltaPhi= DeltaPhi + 2*M_PI;}  //if less then -pi/2 add 2pi
	       if (DeltaPhi > 3*M_PI/2){DeltaPhi = DeltaPhi -2*M_PI;}
	       DeltaEta = (track_eta[n] - track_eta[m]);
	       Corr->Fill(DeltaPhi,DeltaEta);
	     }
	   }
	 }
       }
     }
     
     TCanvas *Corr2D = new TCanvas("Corr2D", "Hadron-Hadron_Correlation", 800,800);
     Corr->SetMinimum(0.); //set Y-Axis minimum to 0 when drawing 
     Corr->Draw("SURF1");
     TFile *MyFile = new TFile("h_h_Correlation.root","RECREATE");
     Corr->Write();
     MyFile->Print(); //write to file so I don't have to run this again
}

void Mixing::divide_hadron(){
  TFile *corr = new TFile("h_h_Correlation.root"); 
  TH2F * h1 = (TH2F*)corr->Get("Correlation"); 
  h1->Sumw2(); 

  TFile *mix = new TFile("h_h_MixedEvents.root");
  TH2F * h2 = (TH2F*)mix->Get("Correlation"); 
  h1->Divide(h2);
  h1->Draw("SURF2");
}


void Mixing::divide_gammah(){
  TFile *corr = new TFile("g_h_Correlation.root"); 
  TH2F * h1 = (TH2F*)corr->Get("Correlation"); 
  h1->Sumw2(); 

  TFile *mix = new TFile("g_h_MixedEvents.root");
  TH2F * h2 = (TH2F*)mix->Get("Correlation"); 
  h1->Divide(h2);
  h1->Draw("SURF2");
}

//need to normalize and also take projections before ratios

void Mixing::CentralityList(){
  
  //This outpus an ordered list of events according to multiplicity. 
  //This will be usefull for centrality binning, as well as event mixing in the usual way

  if (fChain ==0) return;
 
  Long64_t nentries = fChain->GetEntriesFast();
  vector <pair <Long64_t,Float_t> > CentList;// Fill with global entry number and V0M
  //is there a better alterniatve to std::pair?

  for (Long64_t jentry=0; jentry < nentries; jentry++) {//normally stop at nentries
    Long64_t ientry = LoadTree(jentry);  if (ientry < 0) break;
    fChain->GetEntry(jentry); 
    Float_t V0M = 0;                                    
    for (int i = 0; i < 64;i++) //64 channnels for v0
      {
	V0M += multiplicity_v0[i];
      }
    cout<<V0M<<endl;
    if (V0M > 0) CentList.push_back(make_pair(jentry, V0M));
    
  }
  
  cout<<"Multiplicity checks done, ordering according to V0M"<<endl;
  
  //the following sorts according to multiplicty, and fill a vector with global entry
  sort(CentList.begin(), CentList.end(),[](const pair<Long64_t,Float_t> &left, const pair<Long64_t,Float_t> &right) 
       {
	 int a;
	 return(left.second < right.second); //low multiplicity filled first
	 cout<<a<<endl;
	 a++;
       });
  
  cout<<"Sort Complete"<<endl;
  
  for(int i = 0; i < CentList.size(); i++)
    {cout << "(" << CentList[i].first << ", " << CentList[i].second << ")" << endl;}
  
  FILE *out = fopen("CentralityList.txt", "w");
  if(!out) { cout<<"Can't open file"<<endl; }
  
  for(int i=0; i != CentList.size();++i )  {
    fprintf(out, "%llu", CentList[i].first);
    fprintf(out, "\n");}
   
 fclose(out);
}


void Mixing::SimpleMix(){ //   mixes EMCal triggered events. Tirgger bias seen in the form of
                           //   near and away side riges in ∆phi. Not flat as in EMCal-MinBias mixing.
                           //   should ask for trigger information in NTuple.
  if (fChain == 0) return;

  //need to read from text file
  Long64_t nentries = fChain->GetEntriesFast();
  vector <Long64_t> EntryVector;
  
  
  ifstream inputFile("CentralityList.txt");
  if (inputFile){
    Long64_t line;
    while(inputFile >> line){
      EntryVector.push_back(line); //read into a vector line by line
    }
  }

  cout<<"vector filled"<<endl;

  TH2D *MixCorr = new TH2D("Correlation", "Mixed Events", 60,-M_PI/2,3*M_PI/2, 34, -1.7, 1.7);
  MixCorr->Sumw2();
  
  double DeltaPhi;
  double Cluster_rads;
  double Track_rads;
  
  double DeltaEta;
  double Cluster_eta;
  double Track_eta;

  double_t PoolSize = 10;
 
  for (Long64_t v=0; v<20;v++){//loop over 5% centrality bins
    cout<<v<<endl;
    
    for (Long64_t j = v*EntryVector.size()/20; j < (v+1)*EntryVector.size()/20; j++) { //select event j within c-bin
      Long64_t ientry = LoadTree(EntryVector[j]); 
      if (ientry < 0) break;
      
      fChain->GetEntry(EntryVector[j]);
      Float_t E1Cluster_e[ncluster], E1Cluster_phi[ncluster], E1Cluster_eta[ncluster], E1Cluster_ncell[ncluster];
      copy(cluster_e, cluster_e+ncluster, E1Cluster_e);
      copy(cluster_phi, cluster_phi+ncluster, E1Cluster_phi);
      copy(cluster_eta, cluster_eta+ncluster, E1Cluster_eta);
      copy(cluster_ncell, cluster_ncell+ncluster, E1Cluster_ncell);
      ULong64_t NClusters = ncluster;

      cout<<j<<endl;      

      int k= v*EntryVector.size()/20; //start entry k in the same centrality bin as entry j
      double_t MixCounter = 0;      
      
      for (int k = j+1; (k<(v+1)*(EntryVector.size()/20) && MixCounter < PoolSize); k++){
	
	fChain->GetEntry(EntryVector[k]);     
	Float_t E2Track_pt[ntrack], E2Track_phi[ntrack], E2Track_eta[ntrack];
	copy(track_pt, track_pt+ntrack, E2Track_pt);
	copy(track_phi, track_phi+ntrack, E2Track_phi);
	copy(track_eta, track_eta+ntrack, E2Track_eta);
	ULong64_t NTracks = ntrack;

	  for (int n=0; n < ncluster; n++){ //get phi and eta for tracks in event j
	   if ((E1Cluster_e[n] > 6) && (E1Cluster_ncell[n]>1)){
	      	     
	     for (int m=0; m < NTracks; m++){
	       if(E2Track_pt[m] > 1){ 
		 //Cut on 1GeV Tracks

		 DeltaPhi = (E1Cluster_phi[n] - E2Track_phi[m]);
		 if (DeltaPhi < -M_PI/2){DeltaPhi= DeltaPhi + 2*M_PI;}
		 if (DeltaPhi > 3*M_PI/2){DeltaPhi = DeltaPhi - 2*M_PI;}
		 
		 DeltaEta = (E1Cluster_eta[n] - E2Track_eta[m]);
		 
		 MixCorr->Fill(DeltaPhi,DeltaEta,PoolSize);
		 
// 		 if((abs(DeltaPhi)<0.2)&&(abs(DeltaEta)<0.2)){
// 		   cout<<DeltaPhi<<"    "<<DeltaEta<<"    "<<E1Cluster_e[n]<<"    "<<E2Track_eta[m]<<endl;
// 		 }
		 
	       }

	     }//mtracks

	   }//cluster energy check

	  }//nclusters
	  MixCounter++;
	
      }//k 
      
    }//j
    
  }//v
  
  TCanvas *Corr1D = new TCanvas("Corr1D", "Correlation", 800,800);
  MixCorr->SetMinimum(0.);  
  MixCorr->Draw("SURF1");
  TFile *MyFile = new TFile("g_h_MixedEvents.root","RECREATE");
  MixCorr->Write();
  MyFile->Print(); 
}

void Mixing::HadronMix(){ //   mixes EMCal triggered events, but now only track Phi and eta from different events
                           //   same EMCal Bias as SimpleMix()

  if (fChain == 0) return;

  //need to read from text file
  Long64_t nentries = fChain->GetEntriesFast();
  vector <Long64_t> EntryVector; //just read entry #'s from centrality list
  
  
  ifstream inputFile("CentralityList.txt");
  if (inputFile){
    Long64_t line;
    while(inputFile >> line){
      EntryVector.push_back(line); //read into a vector line by line
    }
  }

  cout<<"vector filled"<<endl;
  
  TH2D *MixCorr = new TH2D("Correlation", "Mixed Events", 60,-M_PI/2,3*M_PI/2, 34, -1.7, 1.7); //range -pi/2 to 3pi/2
  MixCorr->Sumw2();
  
  double DeltaPhi;
  double Track_radsj;
  double Track_radsk;
  
  double DeltaEta;
  double Track_etaj;
  double Track_etak;
  double ptj;
  double ptk;

  double_t PoolSize = 10; 
 
  for (Long64_t v=0; v<20;v++){//loop over 5% centrality bins, maybe skip 90-100% centrality (peripherals)
    cout<<v<<endl;
    
    for (Long64_t j = v*EntryVector.size()/20; j < (v+1)*EntryVector.size()/20; j++) {//loop events j within c-bin
      Long64_t ientry = LoadTree(j); if (ientry < 0) break;

      fChain->GetEntry(EntryVector[j]);
      Float_t E1Track_pt[ntrack],E1Track_phi[ntrack],E1Track_eta[ntrack];
      copy(track_pt, track_pt+ntrack, E1Track_pt); //I think copy results in memory leak
      copy(track_phi, track_phi+ntrack, E1Track_phi);
      copy(track_eta, track_eta+ntrack, E1Track_eta);
      ULong64_t NTracks1 = ntrack;
      
      cout<<j<<endl;

      int k= v*EntryVector.size()/20;
      int MixCounter = 0;

      for (int k = v*EntryVector.size()/20; (k<(v+1)*EntryVector.size()/20 && MixCounter < PoolSize); k++){

	fChain->GetEntry(EntryVector[k]);
	Float_t E2Track_pt[ntrack],E2Track_phi[ntrack],E2Track_eta[ntrack];
	copy(track_pt, track_pt+ntrack, E2Track_pt);
	copy(track_phi, track_phi+ntrack, E2Track_phi); //copy to vectors b/c two instances
	copy(track_eta, track_eta+ntrack, E2Track_eta); //of GetEntry() simply didn't work
	ULong64_t NTracks2 = ntrack;	
	
	if (k != j) {
	  
	  for(int n = 0; n<NTracks1; n++){
	    if(E1Track_pt[n] > 1){
	      
	      for (int m=0; m < NTracks2; m++){
		if(E2Track_pt[m] > 1){
		  
		  DeltaPhi = (E1Track_phi[n] - E2Track_phi[m]);
		  if (DeltaPhi < -M_PI/2){DeltaPhi= DeltaPhi + 2*M_PI;}  //if less then -pi/2 add 2pi
		  if (DeltaPhi > 3*M_PI/2){DeltaPhi = DeltaPhi - 2*M_PI;}
		  
		  DeltaEta = E1Track_eta[n] - E2Track_eta[m];
		  
		  MixCorr->Fill(DeltaPhi,DeltaEta,PoolSize);
		  
		  // 		    if((abs(DeltaPhi)<0.2)&&(abs(DeltaEta)<0.2)){
		  // 		    cout<<DeltaPhi<<"    "<<DeltaEta<<"    "<<ptj<<"    "<<ptk<<endl;
		  // 		    }
		} //E2 pt check
		
	      }//mtracks
	      
	    } // E1 track pt check
	    
	  } //ntracks
	    
	    MixCounter++; //only increment mixing counter if mix actually occurs (k!=j)
	    
	}//k!=j
	  
      }//k
	
    }//j
      
  }//v
    
    TCanvas *HCorr2D = new TCanvas("HCorr1D", "Hadron-Hadron", 800,800);
    MixCorr->SetMinimum(0.); //set Y-Axis minimum to 0 when drawing 
    MixCorr->Draw("SURF1");
    TFile *MyFile = new TFile("h_h_MixedEvents.root","RECREATE");
    MixCorr->Write();
    MyFile->Print(); //write to file so I don't have to run this again
}
  
vector <Long64_t> Mixing::entries_after_cuts(Long64_t start, int size){
    
    //used for loops over PoolA and B, IMPORTANT: Pool B starts where A ends
    //this is method can change if I implement different cuts
    //I can change this to simply read from centralitylist.txt, which has the V0M check...
    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t an_entry = start;
    int counter = 0;
    vector <Long64_t> vec;
    while ((counter < size)&&(an_entry<nentries)){
      fChain->GetEntry(an_entry);
      if((centrality_v0m>0) && (centrality_v0m > 0)) //the cut
	{
	  //	  if (cluster_m20[an_entry] <0.3)
	  vec.push_back(an_entry);
	  counter++;
	}
    an_entry++;//iterate regardless of if statement 
  }
  return vec;//return a vector of entries (long64_t) that survive cuts I defined in this method
}
  //this is mostly for when I use small pools to check things. 


void Mixing::PTable(){
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();

  //The first element is the event #
  //The second element is the metric
  vector <vector <pair <Long64_t,Float_t> > > mPreference, wPreference;  
  Float_t D;


  //was using 100 for a while
  int group_size = 500;
  vector <Long64_t> Avector, Bvector;
  Avector = entries_after_cuts(0,group_size); 
  Bvector = entries_after_cuts(Avector.back()+1,Avector.size()); //want Pool A&B to be same size
  
  for (Long64_t A=0; A < Avector.size(); A++) {
    Long64_t ientry = LoadTree(Avector[A]); 
    if (ientry < 0) break;
    
    fChain->GetEntry(Avector[A]);    
    
    valarray<Float_t> multip (multiplicity_v0,64);
    Float_t V0M_i = multip.sum();
    cout<<Avector[A]<<" "<<"A"<<endl;
    vector <pair <Long64_t,Float_t> > temprow; 
    temprow.push_back(make_pair(Avector[A], 0));
    //replace with GetMetric() when full NTuple is available


    for (Long64_t B=0; B < Bvector.size(); B++) {
      ientry = LoadTree(Bvector[B]);
      if (ientry < 0) break;
      
      fChain->GetEntry(Bvector[B]); 
      valarray<Float_t> multip (multiplicity_v0,64);
      Float_t V0M_j = multip.sum();
      D = abs(V0M_i - V0M_j);
      temprow.push_back(make_pair(Bvector[B], D));

    }//loop PoolB
  
    mPreference.push_back(temprow);
    
  }//loop of poolA

  //sorting
  for (int i = 0; i !=mPreference.size(); i++){//Loop over rows
    sort( mPreference[i].begin()+1, mPreference[i].end(),[](const pair<Long64_t,Float_t> &left, 
							    const pair<Long64_t,Float_t> &right) 

	  { return(left.second < right.second); } );    //sort each row by "D"
  }//sort from begin()+1 because first entry in each row is the suitor, not a prefernce


  //writing
  FILE *outm = fopen("mPTable.txt", "w");
  if(!outm) { cout<<"Can't open file"<<endl; }
  
  for(int i=0; i !=  mPreference.size();++i )  {
    
    for (int j=0; j !=  mPreference[i].size(); j++){
      fprintf(outm, "%llu",  mPreference[i][j].first); 
      fprintf(outm, "\t");
    }
    fprintf(outm, "\n");//new line for each row
  }
  
  fclose(outm);
  
  cout<<"Table is"<< mPreference.size()<<"by"<< mPreference[0].size()<<endl;

  cout << "----------Male Preference Table Done------------" <<endl;
  
  //Now I fill the partner table. Pool B is now the outer loop, poolA is inner loop
  
  for (Long64_t B=0; B<Bvector.size(); B++) {
    Long64_t ientry = LoadTree(Bvector[B]); 
    if (ientry < 0) break;
    
    fChain->GetEntry(Bvector[B]);
    cout<<Bvector[B]<<" "<<"B"<<endl;
    vector <pair <Long64_t,Float_t> > temprow; 
    Float_t V0M_i = (centrality_v0m+centrality_v0m)/2;
    temprow.push_back(make_pair(Bvector[B], 0)); 

    //first entry is now "bride", then suitors

    for (Long64_t A = 0; A < Avector.size(); A++) {
      
      fChain->GetEntry(Avector[A]); 
      Float_t V0M_j = (centrality_v0m+centrality_v0m)/2;
      D = abs(V0M_i - V0M_j);
      temprow.push_back(make_pair(Avector[A], D)); 
      
    } //loop A

    wPreference.push_back(temprow);
    
  } //loop B
 
 
  //sorting 
  for (int i = 0; i !=wPreference.size(); i++){//Loop over rows
    sort(wPreference[i].begin()+1, wPreference[i].end(),[](const pair<Long64_t,Float_t> &left, 
							   const pair<Long64_t,Float_t> &right) 
	 {  return(left.second < right.second); } );  
  }
  
  
  // writing 
  FILE *outw = fopen("wPTable.txt", "w");
  if(!outw) { cout<<"Can't open file"<<endl; }
  
  for(int i=0; i != wPreference.size();++i )  {//loop over rows     
    for (int j=0; j != wPreference[i].size(); j++){//loop over elements in row //
      fprintf(outw, "%llu", wPreference[i][j].first); //write global entry number to file
      fprintf(outw, "\t");
    }
    fprintf(outw, "\n");//new line for each row
  }
  
  fclose(outw);
  cout<<"Table is"<<wPreference.size()<<"by"<<wPreference[0].size()<<endl;
  
}

//will have to make a class that obtains "D" once NTupple is updated
//A pressing question at the moment is how do we weight the similarity criteria
 
void Mixing::PairingTest()
{
  vector <vector <Long64_t> > mprefer = GetPreferences("mPTable.txt");
  vector <vector <Long64_t> > wprefer = GetPreferences("wPTable.txt");
  StableMarriage(mprefer,wprefer);
  
}

//get preference table from text file
vector <vector <Long64_t> > Mixing::GetPreferences(const string filename)
{
  vector <vector <Long64_t> > PreferenceVector;
  vector <Long64_t> row; 
  Long64_t a;
  
  istringstream lin;

  ifstream infile(filename);  

  for (string line; std::getline(infile, line); ) {
    lin.clear();
    lin.str(line);
    while (lin >> a) 
      {
	row.push_back(a);
      }
    PreferenceVector.push_back(row); 
    row.clear();
  }
  
  return  PreferenceVector;
}


bool Mixing::Hot_or_Not(vector <Long64_t> prefer, Long64_t first, Long64_t second){
  
  for (int i = 0; i < prefer.size(); i++) 
    { 
      if (prefer[i] == first)
        return true; 
      if (prefer[i] == second)
        return false;
    }  
  return false;  // avoid warning "control may reach end of non-void function"
}

void Mixing::StableMarriage(vector <vector <Long64_t> > men_data, 
			     vector <vector <Long64_t> > women_data) {

  map<Long64_t, vector<Long64_t> > men_pref, women_pref, engaged;
  deque <Long64_t>  bachelors;
  map <Long64_t,vector <Long64_t> > men_pairings;

  Long64_t poolsize = 100;
  
  // init data structures
  for (int i = 0; i < men_data.size(); ++i) // person
    {
      for (int j = 1; j < men_data[0].size(); ++j) // preference
        {
	  men_pref[men_data[i][0]].push_back(men_data[i][j]);
	  women_pref[women_data[i][0]].push_back(women_data[i][j]);
        }
      vector <Long64_t> emptyvec;
      men_pairings.insert(make_pair(men_data[i][0],emptyvec));
      bachelors.push_back(men_data[i][0]);
      cout<<i<<endl;
    }
  
  cout<<"Data Initialized"<<endl;

  while (!bachelors.empty()) {

    Long64_t &suitor = bachelors.front();
    bool unfullfilled = false;
                              
    while ((men_pairings.at(suitor)).size() < poolsize && !unfullfilled){ 
      
      if ((men_pref[suitor]).size() == 0)
	{
	  bachelors.pop_front();
	  cout<<suitor<<"is out of brides to ask!"<<endl;
	  unfullfilled = true;
	  continue;
	}

      vector <Long64_t> &preflist = men_pref[suitor];
      
      Long64_t suitors_engagements = (men_pairings.at(suitor)).size(); 
      //required exit condition for when a suitor rejected by all remaining brides
      
      //      for (vector <Long64_t>::iterator it = preflist.begin(); (it != preflist.end() && (men_pairings.at(suitor)).size() < poolsize); it++)
      for (int n = 0; n<preflist.size(); n++)
	{//something wrong with exit condition...

	  Long64_t &bride = preflist[n];

	  if ( (engaged[bride]).size() < poolsize)
	    {
	      (engaged[bride]).push_back(suitor);
	      (men_pairings.at(suitor)).push_back(bride);
	      
	      cout <<"suitor  " << suitor << "  counter  is  "
		   << (men_pairings.at(suitor)).size()<<"\n"<<"\n";
	    }
	  
	  else {
	    Long64_t worst_groom = Least_Liked(women_pref[bride], engaged[bride]);
	    
	    if (Hot_or_Not(women_pref[bride], suitor, worst_groom)) {
	      //current suitor is better than her worst groom
	      
	      cout << "\t" << bride << " dumped " << worst_groom << " for " << suitor << "\n";
	      
	      if ((men_pairings.at(worst_groom)).size() == poolsize) bachelors.push_back(worst_groom);
	      
	     (men_pairings.at(worst_groom)).erase(remove((men_pairings.at(worst_groom)).begin(), 
	      (men_pairings.at(worst_groom)).end(), bride),(men_pairings.at(worst_groom)).end());
	    //remove bride from worst_groom match list
	   
	
	     //	      (men_pref[worst_groom]).erase(remove((men_pref[worst_groom]).begin(),
	     //   (men_pref[worst_groom]).end(), bride),(men_pref[worst_groom]).end());
	      //do the same for the preflist
	     //the above is uneeded if i remove from preflist on first passthrough

	     //	      cout<<worst_groom<<" has "<<(men_pref[worst_groom]).size()<<" left on his preflist"<<endl;

	      cout<<"\t" << "suitor "<< worst_groom <<" has "
		  <<(men_pairings.at(worst_groom)).size()<<" brides"<<endl;

	      //So this is a key difference from the original gale-shapely algorithm
	      //every suitor is no longer guranteed to be fullfilled. 

	      //So condition for "stable" matching is rarely realized
	      //To this end, "stability" is hardcoded my removing the women 
	      //the dumped worst groom from from HIS plist.
	      

	      
	      (engaged[bride]).erase(remove((engaged[bride]).begin(), 
			    (engaged[bride]).end(), worst_groom),(engaged[bride]).end()); 
	      //remove groom from brides matchlist

	      //FIXME: This is ugly. please remember to change from map of vectors to
	      //to map of unordered SETs. Including maybe the queue, but hopefully not. Also should
	      //maybe make men_pref a set as well. Because removing anything from  the queue is rare
	      // I think i should leave it as a queue. The men_pairings and engaged structures should
	      //be changed.
	      
	      
	      (engaged[bride]).push_back(suitor);
	      
	      (men_pairings.at(suitor)).push_back(bride);
	      
	      cout<<"\t" << "suitor "<< suitor <<" has "
		  <<(men_pairings.at(suitor)).size()<<" brides"<<endl<<endl;
	      
	    }//hot_or_not

	  
	  }//else (bride is full)
	  
	  (men_pref[suitor]).erase(remove((men_pref[suitor]).begin(),
      	      (men_pref[suitor]).end(), bride),(men_pref[suitor]).end());

	  //	   cout<<suitor<<" now has "<<(men_pref[suitor]).size()<<" left on his preflist (ev) "<<endl;


	  if ((men_pairings.at(suitor)).size() == poolsize) 
	    {
	      bachelors.pop_front(); 
	      cout<<endl<<"break"<<endl;
	      break; //haha it's fucking break!!!!! of course. 
	    }
	  

	  // I no longer care about perfect stability, then there is no need for him
	  //to propose to this bride again! He will stay with this bride as long as
	  //he is not a worst_groom that gets kicked out!
	  

	   //Neither Manipulating the queue 
	}//preflist  

      //If the suitor has no new wife after going through his entire preflist
      if(suitors_engagements == (men_pairings.at(suitor)).size())
	{
	  cout<<bachelors.front();
	  bachelors.pop_front();
	  cout <<endl << suitor<<" is unfullfilled. Queue now has " << bachelors.size()<< " suitors" <<endl;
	  unfullfilled = true;
	  
	}
      
    } //men_pairings size
      //   if (suitor ==108) break; //works (uh, not really), suitors have set of unique wives
    
    
  }//bachelors queue
    

  cout << "Engagements:\n";

  FILE *outm = fopen("Engagements.txt", "w");
  if(!outm) { cout<<"Can't open file"<<endl; }

  for(map<Long64_t,vector <Long64_t> >::iterator iter = men_pairings.begin(); iter != men_pairings.end(); ++iter)
    { 
      fprintf(outm, "%llu", iter->first);
      fprintf(outm, "\t");

      for (int i = 0; i != (iter->second).size(); i++)
	{
	  fprintf(outm, "%llu", (iter->second)[i]);
	  fprintf(outm, "\t");
	}
      fprintf(outm, "\n");
    }

  //each suitor has the right number of brides!!!!
  //but repeats a bride several times

  fclose(outm);



}


Long64_t Mixing::Least_Liked(vector <Long64_t> &bride_pref, vector <Long64_t> &brides_men)
{
  //  for (vector <Long64_t>::reverse_iterator rit = bride_pref.rbegin();                       
  //    rit != bride_pref.rend(); rit ++)   {                                                   


  sort(brides_men.begin(), brides_men.end(),[&bride_pref](const Long64_t &left,
                                                          const Long64_t &right)
       {
         for (int i = 0; i < bride_pref.size(); i++)
           {
             if (bride_pref[i] == left)
               return true; //if the first element comes first in pref, do nothing              

             if (bride_pref[i] == right)
               return false; //if the second element comes first, std::sort will swap           
           }
         //this will sort the bride's current grooms from most to least liked                   

         return true;
       } );

  return brides_men.back();
}



void Mixing::pi_zero_mass()
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  TH1D *InvMass = new TH1D("Correlation", "m_{#gamma#gamma} 5-13 GeV 50-90 Centrality ", 300,0,0.75);
  TH1D *OpnAngle = new TH1D("OpeningAngle","#theta 5-13 GeV 50-90 Centrality",72,0,M_PI);

  TLorentzVector Lnvector;
  TLorentzVector Livector;


  vector <Long64_t> EntryVector;

  ifstream inputFile("CentralityList.txt");
  if (inputFile){
    Long64_t line;
    while(inputFile >> line){
     EntryVector.push_back(line); //read into a vector line by line                                                                 
    }
  }
  cout<<"EntryVector Filled"<<endl;
  // for (int e = 0; e < Lowest_M; e++){
  // Long64_t jentry = EntryVector[e];
  // int Lowest_M = EntryVector.size()/2; //100-50% centrality
  //  cout<<Lowest_M<<endl;

  Long64_t multentry = EntryVector[(EntryVector.size()/10)];
  fChain->GetEntry(multentry);
  Long64_t mult = (centrality_v0m+centrality_v0m)/2;
  cout<<"Multiplicity "<<mult<<endl;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {    
    Long64_t ientry = LoadTree(jentry);
    cout<<jentry<<endl;
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;   
    //    cout<<"trigger?"<<mixed_event<<endl;
    if(((centrality_v0m+centrality_v0m)/2 < 1325) && ((centrality_v0m+centrality_v0m)/2 > 87)){//80-100% (383) 50+ (1325), 0-90% (87)
      for (int n=0; n< ncluster; n++){
	if (cluster_ncell[n]>1){
	  if (cluster_pt[n] > 0.7){
	    if ((cluster_lambda_square[n][0]< 0.3) && (cluster_lambda_square[n][0] > 0.1)){ //Lambda_20 cut
	      for(int i=0; i<ncluster; i++){
		if (i != n) {
		  if (cluster_pt[i] > 0.7){
		    Float_t assym = abs(((cluster_e[i] - cluster_e[n])/(cluster_e[i] + cluster_e[n])));
		    if ((assym < 0.8) && (assym>0)) { //asymmetry cut
		      //			if(abs(cluster_phi[n]-cluster_phi[i]) > 0.03){
		      //			  if (abs(cluster_eta[n]-cluster_eta[i]) > 0.02){
		      Lnvector.SetPtEtaPhiE(cluster_pt[n], cluster_eta[n], cluster_phi[n], cluster_e[n]);
		      Livector.SetPtEtaPhiE(cluster_pt[i], cluster_eta[i], cluster_phi[i], cluster_e[i]);
		      Float_t Sum_pT = (Lnvector + Livector).Pt();
		      Double_t angl = (Lnvector.Vect()).Angle(Livector.Vect());
		      OpnAngle->Fill(angl);
		      if (angl>0.017){
			if ((Sum_pT > 5) && (Sum_pT < 13)){
			  InvMass->Fill((Lnvector+Livector).M());
			  // }
			  //	}
			}
		      }
		    }
 		  }
 		}
 	      }
 	    }
 	  } 
 	}
      }   
    }
  }  
  
  TCanvas *Mass = new TCanvas("Mass", "piMass", 800,800);
  InvMass->Draw(); 
  //  OpnAngle->Draw();
  TFile *MyFile = new TFile("pi_mass_50_90_Cent.root","RECREATE");
  InvMass->Write();
  //FIXME change root file names back
}


void Mixing::mass_mixing()
{
  TH1D *InvMass = new TH1D("Correlation", "Mass", 300,0,0.75);
  TH1D *OpnAngle = new TH1D("OpeningAngle","#theta 5-13 GeV 50-90 Centrality",72,0,M_PI);

  vector <vector <Long64_t> > Engagements= GetPreferences("Engagements.txt");

  TLorentzVector Lnvector;
  TLorentzVector Livector;

  double DeltaPhi;
  double Cluster_radsj;
  double Cluster_radsk;

  double DeltaEta;
  double Cluster_etaj;
  double Cluster_etak;
  double ptj;
  double ptk;
  
  int counter = 0;
  double_t PoolSize = 10;

  for (int i = 0; i < Engagements.size(); ++i) // row enumerated by glabal entry number                                                      
    {

      //mix engagements[i][0] with events from inner loop
      //do something like 50 tracks for each event for speed
      //I would much rather just use GetEntry() when needed, but switching events is non-trivial. Copy is used at the moment. 
      
      Long64_t ientry = LoadTree(Engagements[i][0]); if (ientry < 0) break;
      
      
      fChain->GetEntry(Engagements[i][0]);
      if(((centrality_v0m+centrality_v0m)/2 < 1325) && ((centrality_v0m+centrality_v0m)/2 > 87)) {//80-100% (383) 50+ (1325), be >87 for 0-90%
	  //	  Float_t Mult1 = centrality_v0m+centrality_v0m/2;
	  //	  cout<<Mult1<<endl;
	  //	  cout<<(centrality_v0m+centrality_v0m)/2<<endl;
	  Float_t E1Cluster_pt[ncluster],E1Cluster_phi[ncluster],E1Cluster_eta[ncluster];
	  Float_t E1Cluster_e[ncluster],E1Cluster_m02[ncluster][2];
	  copy(cluster_pt, cluster_pt+ncluster, E1Cluster_pt);
	  copy(cluster_phi, cluster_phi+ncluster, E1Cluster_phi);
	  copy(cluster_eta, cluster_eta+ncluster, E1Cluster_eta);
	  copy(cluster_e, cluster_e+ncluster, E1Cluster_e);
	  //	  copy(cluster_m02, cluster_m02+ncluster, E1Cluster_m02);
	  copy(&cluster_lambda_square[0][0], &cluster_lambda_square[0][0]+ncluster*2,&E1Cluster_m02[0][0]);
	  ULong64_t NClusters1 = ncluster;

	  cout<<counter<<endl;
	  counter++;	  
	  
	  for (int j = 1; j < Engagements[i].size(); ++j) { // events paired with above event
	    fChain->GetEntry(Engagements[i][j]);
	    if((centrality_v0m+centrality_v0m)/2 < 1325){
	      //  cout<<(centrality_v0m+centrality_v0m)/2<<endl;
	      Float_t E2Cluster_pt[ncluster],E2Cluster_phi[ncluster],E2Cluster_eta[ncluster];
	      Float_t E2Cluster_e[ncluster],E2Cluster_m02[ncluster][2];
	      copy(cluster_pt, cluster_pt+ncluster, E2Cluster_pt);
	      copy(cluster_phi, cluster_phi+ncluster, E2Cluster_phi); //copy to vectors b/c two instances 
	      copy(cluster_eta, cluster_eta+ncluster, E2Cluster_eta); //of GetEntry() simply didn't work
	      copy(cluster_e, cluster_e+ncluster, E2Cluster_e);
	      //	      copy(cluster_m02, cluster_m02+ncluster, E2Cluster_m02);
	      copy(&cluster_lambda_square[0][0], &cluster_lambda_square[0][0]+ncluster*2,&E2Cluster_m02[0][0]);
	      ULong64_t Nclusters2 = ncluster;
	      

	      for (int n=0; n<NClusters1; n++){
		//if (cluster_ncell[n]>1){
		  if (E1Cluster_pt[n] > 0.7){
		    if ((E1Cluster_m02[n][0]< 0.3) && (E1Cluster_m02[n][0] > 0.1)){ //Lambda_20 cut
			for(int i=0; i<Nclusters2; i++){
			  if (i != n) {
			    if (E2Cluster_pt[i] > 0.7){
			      if ((E2Cluster_m02[i][0]< 0.3) && (E2Cluster_m02[i][0] > 0.1)){
				//if (((E1Cluster_pt[n]+E2Cluster_pt[i]) > 5) && ((E1Cluster_pt[n]+E2Cluster_pt[i]) < 13)){ //pT 5-13 peak
				  Float_t assym = abs(((E2Cluster_e[i] - E1Cluster_e[n])/(E2Cluster_e[i] + E1Cluster_e[n])));
				  if ((assym < 0.8) && (assym>0)) { //asymmetry cut
				    // if(abs(E1Cluster_phi[n]-E2Cluster_phi[i]) > 0.03){ //PHI CUT
				    // if (abs(E1Cluster_eta[n]-E2Cluster_eta[i]) > 0.02){ //ETA CUT
				    Lnvector.SetPtEtaPhiE(E1Cluster_pt[n], E1Cluster_eta[n], E1Cluster_phi[n], E1Cluster_e[n]);
				    Livector.SetPtEtaPhiE(E2Cluster_pt[i], E2Cluster_eta[i], E2Cluster_phi[i], E2Cluster_e[i]);
				    Float_t Sum_pT = (Lnvector + Livector).Pt();
				    Double_t angl = (Lnvector.Vect()).Angle(Livector.Vect());
				    OpnAngle->Fill(angl);
				    if (angl>0.017){
				      if ((Sum_pT > 5) && (Sum_pT < 13)){
					InvMass->Fill((Lnvector+Livector).M());

				      // }
				      // }
				      }
				    }	      //might want to implement delta eta and delta phi cuts 0.02 cut in deltaEta and 0.03 in deltaPhi
				  }           // the above is soley to remove charged tracks. need to implement OPENENING angle cuts (Lorentz vector)
			      }
			    }
			  }
			}
		    }
		  }
	      }
	    }
	  }
      }
    }
  



  TFile *MyFile = new TFile("Mass_mix_50_90.root","RECREATE");
  InvMass->Write();

  TCanvas *Mass = new TCanvas("Mass", "piMass", 800,800);
  InvMass->Draw();

  TCanvas *angl = new TCanvas("OpnAngl", "#theta", 800,800);
  OpnAngle->Draw();

}

void Mixing::divide_mass(){
  TCanvas *Masses = new TCanvas("Masses", "piMass", 800,800);
  TFile *corr = new TFile("pi_mass_Lorrentz_50_90.root");
  TH1F * h1 = (TH1F*)corr->Get("Correlation");
 

  TFile *mix = new TFile("Mass_mix_50_90.root");
  TH1F * h2 = (TH1F*)mix->Get("Correlation");
  h2->Scale(6.50513/8);//without eta and phi cuts, 1.6132011375 with rebinnign 4
  //swiching to lorentz vector, 6.50513 base scale
  //Lorentz and Opening angle, the scale is 6.4799
  h1->SetTitle(" Galey-Shapley Mixing for m_{#gamma#gamma} Background; m_{#gamma#gamma} [GeV] ; counts"); 
  //  h2->Scale(0.8099875);
  h2->Rebin(8);
  h2->SetLineColor(2);
  h2->SetLineStyle(1);  
  h1->Draw();
  h2->Draw("same");
//  h1->SetTitle("Ratio;  m_{#gamma#gamma} [GeV] ; counts");
//  h1->Divide(h2);
//  h1->Draw();
  auto legend = new TLegend(0.1,0.7,0.3,0.75);
  // legend->SetHeader("Mixing and Same Event","C"); // option "C" allows to center the header
  legend->AddEntry(h1,"Same Event","l");
  legend->AddEntry(h2,"Mixed Event","l");
  legend->Draw();
}

void Mixing::Simple_mass_mixing()
{
  TH1D *InvMass = new TH1D("Correlation", "Mass", 300,0,0.75);

  //  vector <vector <Long64_t> > Engagements= GetPreferences("Engagements.txt");

  TLorentzVector Lnvector;
  TLorentzVector Livector;

  double DeltaPhi;
  double Cluster_radsj;
  double Cluster_radsk;

  double DeltaEta;
  double Cluster_etaj;
  double Cluster_etak;
  double ptj;
  double ptk;

  int counter = 0;
  double_t PoolSize = 10;

  vector <Long64_t> EntryVector;

  ifstream inputFile("CentralityList.txt");
  if (inputFile){
    Long64_t line;
    while(inputFile >> line){
      EntryVector.push_back(line); //read into a vector line by line                                                                                              
    }
  }
  cout<<"EntryVector Filled"<<endl;

  for (Long64_t i=0; i<EntryVector.size();i++) {
    Long64_t ientry = LoadTree(EntryVector[i]); if (ientry < 0) break;
    
      fChain->GetEntry(EntryVector[i]);
      if((centrality_v0m+centrality_v0m)/2 < 1325){//80-100% centrality (383) 50+ (1325), be >87 for 0-90%                                                                                                                                                                      
	Float_t E1Cluster_pt[ncluster],E1Cluster_phi[ncluster],E1Cluster_eta[ncluster];
	Float_t E1Cluster_e[ncluster],E1Cluster_m02[ncluster][2];
	copy(cluster_pt, cluster_pt+ncluster, E1Cluster_pt);
	copy(cluster_phi, cluster_phi+ncluster, E1Cluster_phi);
	copy(cluster_eta, cluster_eta+ncluster, E1Cluster_eta);
	copy(cluster_e, cluster_e+ncluster, E1Cluster_e);
	//	copy(cluster_m02, cluster_m02+ncluster, E1Cluster_m02);
	copy(&cluster_lambda_square[0][0], &cluster_lambda_square[0][0]+ncluster*2,&E1Cluster_m02[0][0]);
	ULong64_t NClusters1 = ncluster;

	cout<<counter<<endl;
	counter++;

	for (int j = 0; j < EntryVector.size(); j++) { // events paired with above event
	  if (j!=i){
	    fChain->GetEntry(EntryVector[j]);
	    if((centrality_v0m+centrality_v0m)/2 < 1325){
	      //  cout<<(centrality_v0m+centrality_v0m)/2<<endl;                                                                                                                                                                                                                
	      Float_t E2Cluster_pt[ncluster],E2Cluster_phi[ncluster],E2Cluster_eta[ncluster];
	      Float_t E2Cluster_e[ncluster],E2Cluster_m02[ncluster][2];
	      copy(cluster_pt, cluster_pt+ncluster, E2Cluster_pt);
	      copy(cluster_phi, cluster_phi+ncluster, E2Cluster_phi); //copy to vectors b/c two instances                                                                                                                                                                                 
	      copy(cluster_eta, cluster_eta+ncluster, E2Cluster_eta); //of GetEntry() simply didn't work                                                                                                                                                                                  
	      copy(cluster_e, cluster_e+ncluster, E2Cluster_e);
	      //	      copy(cluster_m02, cluster_m02+ncluster, E2Cluster_m02);
	      copy(&cluster_lambda_square[0][0], &cluster_lambda_square[0][0]+ncluster*2,&E2Cluster_m02[0][0]);
	      ULong64_t Nclusters2 = ncluster;
	      
	      
	      for (int n=0; n<NClusters1; n ++){
		//if (cluster_ncell[n]>1){                                                                                                                                                                                                                                                
		if (E1Cluster_pt[n] > 0.7){
		  if ((E1Cluster_m02[n][0]< 0.3) && (E1Cluster_m02[n][0] > 0.1)){ //Lambda_20 cut                                                                                                                                                                                             
		    for(int i=0; i<Nclusters2; i++){
		      if (i != n) {
			if (E2Cluster_pt[i] > 0.7){
			  if ((E2Cluster_m02[i][0]< 0.3) && (E2Cluster_m02[i][0] > 0.1)){
			    if (((E1Cluster_pt[n]+E2Cluster_pt[i]) > 5) && ((E1Cluster_pt[n]+E2Cluster_pt[i]) < 13)){ //pT 5-13 peak                                                                                                                                                  
			      Float_t assym = abs(((E2Cluster_e[i] - E1Cluster_e[n])/(E2Cluster_e[i] + E1Cluster_e[n])));
			      if ((assym < 0.8) && (assym>0)) { //asymmetry cut                                                                                                                                                                                                       
				//    if(abs(E1Cluster_phi[n]-E2Cluster_phi[i]) > 0.03){                                                                                                                                                                                              
				// if (abs(E1Cluster_eta[n]-E2Cluster_eta[i]) > 0.02){                                                                                                                                                                                                
				Lnvector.SetPtEtaPhiE(E1Cluster_pt[n], E1Cluster_eta[n], E1Cluster_phi[n], E1Cluster_e[n]);
				Livector.SetPtEtaPhiE(E2Cluster_pt[i], E2Cluster_eta[i], E2Cluster_phi[i], E2Cluster_e[i]);
				InvMass->Fill((Lnvector+Livector).M());
				// }                                                                                                                                                                                                                                              
				//          }                                                                                                                                                                                                                                     
			      }           //need to implement delta eta and delta phi cuts 0.02 cut in deltaEta and 0.03 in deltaPhi
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
  }




  TFile *MyFile = new TFile("Mass_mix.root","RECREATE");
  InvMass->Write();

  TCanvas *Mass = new TCanvas("Mass", "piMass", 800,800);
  InvMass->Draw();

}






void Mixing::GaleHadronMix(){

  if (fChain == 0) return;

  //need to read from text file                                                                                                                                                                                                                                                                                                        
  Long64_t nentries = fChain->GetEntriesFast();
  vector <Long64_t> Entryvector; //Just Read Entry #'S From centrality list                                                                                                                                                                       
  vector <vector <Long64_t> > Engagements= GetPreferences("Engagements.txt");
  //  for (vector <vector<Long64_t> >::iterator it = Engagements.begin(); it != Engagements.end(), it++)
  //    {
  //      cout<<it[0][0]<<endl;
  //    }

 


  cout<<"vector filled"<<endl;

  TH2D *MixHadron = new TH2D("Hadron Correlation", "Hadron Mixed Events (GS)", 60,-M_PI/2,3*M_PI/2, 34, -1.7, 1.7); //range -pi/2 to 3pi/2 

  MixHadron->Sumw2();

  double DeltaPhi;
  double Track_radsj;
  double Track_radsk;

  double DeltaEta;
  double Track_etaj;
  double Track_etak;
  double ptj;
  double ptk;

  double_t PoolSize = 10;

  for (int i = 0; i < Engagements.size(); ++i) // row enumerated by glabal entry number
    {

      //mix engagements[i][0] with events from inner loop
      //do something like 50 tracks for each event for speed
   
      
      
      
      Long64_t ientry = LoadTree(Engagements[i][0]); if (ientry < 0) break;
      
      fChain->GetEntry(Engagements[i][0]);
      Float_t E1Track_pt[ntrack],E1Track_phi[ntrack],E1Track_eta[ntrack];
      copy(track_pt, track_pt+ntrack, E1Track_pt); //I think copy results in memory leak                                                                                                                                                                                                                                               
      copy(track_phi, track_phi+ntrack, E1Track_phi);
      copy(track_eta, track_eta+ntrack, E1Track_eta);
      ULong64_t NTracks1 = ntrack;
      
      cout<<Engagements[i][0]<<endl;
      
      int MixCounter = 0;
      
      
      for (int j = 1; j < Engagements[i].size(); ++j) // events paired with above event                                                            
	{
	    
	    
	  fChain->GetEntry(Engagements[i][j]);
	  Float_t E2Track_pt[ntrack],E2Track_phi[ntrack],E2Track_eta[ntrack];
	  copy(track_pt, track_pt+ntrack, E2Track_pt);
	  copy(track_phi, track_phi+ntrack, E2Track_phi); //copy to vectors b/c two instances                                                                                                                                                                                                                                            
	  copy(track_eta, track_eta+ntrack, E2Track_eta); //of GetEntry() simply didn't work                                                                                                                                                                                                                                             
	  ULong64_t NTracks2 = ntrack;
	    
	    
	  for(int n = 0; n<NTracks1; n++){//in the future wich larger mixing pools, use n+=50
	    if(E1Track_pt[n] > 1){
	            
	      for (int m=0; m < NTracks2; m++){
		if(E2Track_pt[m] > 1){
		    
		  DeltaPhi = (E1Track_phi[n] - E2Track_phi[m]);
		  if (DeltaPhi < -M_PI/2){DeltaPhi= DeltaPhi + 2*M_PI;}  //if less then -pi/2 add 2pi                                                                                                                                                                                                                                  
		  if (DeltaPhi > 3*M_PI/2){DeltaPhi = DeltaPhi - 2*M_PI;}
		    
		  DeltaEta = E1Track_eta[n] - E2Track_eta[m];
		    
		  MixHadron->Fill(DeltaPhi,DeltaEta,PoolSize);
		    
		  //                if((abs(DeltaPhi)<0.2)&&(abs(DeltaEta)<0.2)){                                                                                                                                                                                                                                                      
		  //                cout<<DeltaPhi<<"    "<<DeltaEta<<"    "<<ptj<<"    "<<ptk<<endl;                                                                                                                                                                                                                                  
		  //                }                                                                                                                                                                                                                                                                                                  
		} //E2 pt check                                                                                                                                                                                                                                                                                                        
		
	      }//mtracks                                                                                                                                                                                                                                                                                                               
	            
	    } // E1 track pt check                                                                                                                                                                                                                                                                                                     
	        
	  } //ntracks                                                                                                                                                                                                                                                                                                                  
	    
	}//k                                                                                                                                                                                                                                                                                                                             
      
    }//j                                                                                                                                                                                                                                                                                                                               
 
 
 
  TCanvas *HCorr2D = new TCanvas("HCorr1D", "Hadron-Hadron", 800,800);
  MixHadron->SetMinimum(0.); //set Y-Axis minimum to 0 when drawing                                                                                                                                                                                                                                                                    
  MixHadron->Draw("SURF1");
  TFile *MyFile = new TFile("h_h_MixedEvents.root","RECREATE");
  MixHadron->Write();
  MyFile->Print(); //write to file so I don't have to run this again                                                                                                                                                                                                                                                                 
}











