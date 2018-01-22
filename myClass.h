//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed May 31 15:47:31 2017 by ROOT version 6.08/00
// from TTree _tree_event/
// found on file: /project/projectdirs/alice/ylai/000245683.root
//////////////////////////////////////////////////////////

#ifndef myClass_h
#define myClass_h


#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class myClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run_number;
   Char_t          mixed_event;
   Float_t         multiplicity_v0a[4];
   Float_t         multiplicity_v0c[4];
   Float_t         centrality[41];
   Bool_t          has_misalignment_matrix;
   Int_t           eg_ntrial;
   Float_t         eg_perp_hat;
   Float_t         eg_cross_section;
   Double_t        primary_vertex[3];
   ULong64_t       ncluster;
   Float_t         cluster_e[2298];   //[ncluster]
   Float_t         cluster_pt[2298];   //[ncluster]
   Float_t         cluster_eta[2298];   //[ncluster]
   Float_t         cluster_phi[2298];   //[ncluster]
   Float_t         cluster_m02[2298];   //[ncluster]
   Float_t         cluster_m20[2298];   //[ncluster]
   Float_t         cluster_tof[2298];   //[ncluster]
   Int_t           cluster_ncell[2298];   //[ncluster]
   ULong64_t       ntrack;
   Float_t         track_e[3980];   //[ntrack]
   Float_t         track_pt[3980];   //[ntrack]
   Float_t         track_eta[3980];   //[ntrack]
   Float_t         track_phi[3980];   //[ntrack]
   ULong64_t       nmc_truth;
   Float_t         mc_truth_e[1];   //[nmc_truth]
   Float_t         mc_truth_pt[1];   //[nmc_truth]
   Float_t         mc_truth_eta[1];   //[nmc_truth]
   Float_t         mc_truth_phi[1];   //[nmc_truth]
   Int_t           mc_truth_label[1];   //[nmc_truth]
   Int_t           mc_truth_pdg_id[1];   //[nmc_truth]
   Int_t           mc_truth_status[1];   //[nmc_truth]
   ULong64_t       njet;
   Float_t         jet_e_raw[41];   //[njet]
   Float_t         jet_e[41];   //[njet]
   Float_t         jet_e_charged[41];   //[njet]
   Float_t         jet_pt_raw[41];   //[njet]
   Float_t         jet_pt[41];   //[njet]
   Float_t         jet_pt_charged[41];   //[njet]
   Float_t         jet_eta[41];   //[njet]
   Float_t         jet_phi[41];   //[njet]
   Int_t           jet_truth_index_z_truth[41][2];   //[njet]
   Float_t         jet_truth_z_truth[41][2];   //[njet]
   Int_t           jet_truth_index_z_reco[41][2];   //[njet]
   Float_t         jet_truth_z_reco[41][2];   //[njet]

   // List of branches
   TBranch        *b_run_number;   //!
   TBranch        *b_mixed_event;   //!
   TBranch        *b_multiplicity_v0a;   //!
   TBranch        *b_multiplicity_v0c;   //!
   TBranch        *b_centrality;   //!
   TBranch        *b_has_misalignment_matrix;   //!
   TBranch        *b_eg_ntrial;   //!
   TBranch        *b_eg_perp_hat;   //!
   TBranch        *b_eg_cross_section;   //!
   TBranch        *b_primary_vertex;   //!
   TBranch        *b_ncluster;   //!
   TBranch        *b_cluster_e;   //!
   TBranch        *b_cluster_pt;   //!
   TBranch        *b_cluster_eta;   //!
   TBranch        *b_cluster_phi;   //!
   TBranch        *b_cluster_m02;   //!
   TBranch        *b_cluster_m20;   //!
   TBranch        *b_cluster_tof;   //!
   TBranch        *b_cluster_ncell;   //!
   TBranch        *b_ntrack;   //!
   TBranch        *b_track_e;   //!
   TBranch        *b_track_pt;   //!
   TBranch        *b_track_eta;   //!
   TBranch        *b_track_phi;   //!
   TBranch        *b_nmc_truth;   //!
   TBranch        *b_mc_truth_e;   //!
   TBranch        *b_mc_truth_pt;   //!
   TBranch        *b_mc_truth_eta;   //!
   TBranch        *b_mc_truth_phi;   //!
   TBranch        *b_mc_truth_label;   //!
   TBranch        *b_mc_truth_pdg_id;   //!
   TBranch        *b_mc_truth_status;   //!
   TBranch        *b_njet;   //!
   TBranch        *b_jet_e_raw;   //!
   TBranch        *b_jet_e;   //!
   TBranch        *b_jet_e_charged;   //!
   TBranch        *b_jet_pt_raw;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_pt_charged;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_truth_index_z_truth;   //!
   TBranch        *b_jet_truth_z_truth;   //!
   TBranch        *b_jet_truth_index_z_reco;   //!
   TBranch        *b_jet_truth_z_reco;   //!

   myClass(TTree *tree=0);
   virtual ~myClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Eccentricity();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     Correlation();
   virtual void     H_HCorrelation();
   virtual vector   <Long64_t> entries_after_cuts(Long64_t start, int size);

   virtual void     PTable();
   virtual void     CentralityList();
   virtual void     SimpleMix();
   virtual void     HadronMix();
   virtual void     divide_hadron();
   virtual void     divide_gammah();
   
   virtual void     PairingTest(); 
   virtual vector   <vector <Long64_t> > GetPreferences(const string filename);
   virtual bool     Hot_or_Not(vector <Long64_t> prefer, Long64_t first, Long64_t second);

   virtual void     StableMarriage(vector <vector <Long64_t> > men_data, 
				   vector <vector <Long64_t> > women_data);

   virtual Long64_t Least_Liked(vector <Long64_t> &bride_pref,
                                vector <Long64_t> &brides_men);

   virtual void     pi_zero_mass();

   virtual void     GaleHadronMix();
   virtual void     mass_mixing();
   virtual void     divide_mass();
   virtual void     Simple_mass_mixing();

};

#endif

#ifdef myClass_cxx
myClass::myClass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("000245683.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("000245683.root");
      }
      f->GetObject("_tree_event",tree);

   }
   Init(tree);
}

myClass::~myClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t myClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t myClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   } 
   return centry;
}

void myClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run_number", &run_number, &b_run_number);
   fChain->SetBranchAddress("mixed_event", &mixed_event, &b_mixed_event);
   fChain->SetBranchAddress("multiplicity_v0a", multiplicity_v0a, &b_multiplicity_v0a);
   fChain->SetBranchAddress("multiplicity_v0c", multiplicity_v0c, &b_multiplicity_v0c);
   fChain->SetBranchAddress("centrality", centrality, &b_centrality);
   fChain->SetBranchAddress("has_misalignment_matrix", &has_misalignment_matrix, &b_has_misalignment_matrix);
   fChain->SetBranchAddress("eg_ntrial", &eg_ntrial, &b_eg_ntrial);
   fChain->SetBranchAddress("eg_perp_hat", &eg_perp_hat, &b_eg_perp_hat);
   fChain->SetBranchAddress("eg_cross_section", &eg_cross_section, &b_eg_cross_section);
   fChain->SetBranchAddress("primary_vertex", primary_vertex, &b_primary_vertex);
   fChain->SetBranchAddress("ncluster", &ncluster, &b_ncluster);
   fChain->SetBranchAddress("cluster_e", cluster_e, &b_cluster_e);
   fChain->SetBranchAddress("cluster_pt", cluster_pt, &b_cluster_pt);
   fChain->SetBranchAddress("cluster_eta", cluster_eta, &b_cluster_eta);
   fChain->SetBranchAddress("cluster_phi", cluster_phi, &b_cluster_phi);
   fChain->SetBranchAddress("cluster_m02", cluster_m02, &b_cluster_m02);
   fChain->SetBranchAddress("cluster_m20", cluster_m20, &b_cluster_m20);
   fChain->SetBranchAddress("cluster_tof", cluster_tof, &b_cluster_tof);
   fChain->SetBranchAddress("cluster_ncell", cluster_ncell, &b_cluster_ncell);
   fChain->SetBranchAddress("ntrack", &ntrack, &b_ntrack);
   fChain->SetBranchAddress("track_e", track_e, &b_track_e);
   fChain->SetBranchAddress("track_pt", track_pt, &b_track_pt);
   fChain->SetBranchAddress("track_eta", track_eta, &b_track_eta);
   fChain->SetBranchAddress("track_phi", track_phi, &b_track_phi);
   fChain->SetBranchAddress("nmc_truth", &nmc_truth, &b_nmc_truth);
   fChain->SetBranchAddress("mc_truth_e", &mc_truth_e, &b_mc_truth_e);
   fChain->SetBranchAddress("mc_truth_pt", &mc_truth_pt, &b_mc_truth_pt);
   fChain->SetBranchAddress("mc_truth_eta", &mc_truth_eta, &b_mc_truth_eta);
   fChain->SetBranchAddress("mc_truth_phi", &mc_truth_phi, &b_mc_truth_phi);
   fChain->SetBranchAddress("mc_truth_label", &mc_truth_label, &b_mc_truth_label);
   fChain->SetBranchAddress("mc_truth_pdg_id", &mc_truth_pdg_id, &b_mc_truth_pdg_id);
   fChain->SetBranchAddress("mc_truth_status", &mc_truth_status, &b_mc_truth_status);
   fChain->SetBranchAddress("njet", &njet, &b_njet);
   fChain->SetBranchAddress("jet_e_raw", jet_e_raw, &b_jet_e_raw);
   fChain->SetBranchAddress("jet_e", jet_e, &b_jet_e);
   fChain->SetBranchAddress("jet_e_charged", jet_e_charged, &b_jet_e_charged);
   fChain->SetBranchAddress("jet_pt_raw", jet_pt_raw, &b_jet_pt_raw);
   fChain->SetBranchAddress("jet_pt", jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_pt_charged", jet_pt_charged, &b_jet_pt_charged);
   fChain->SetBranchAddress("jet_eta", jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_truth_index_z_truth", jet_truth_index_z_truth, &b_jet_truth_index_z_truth);
   fChain->SetBranchAddress("jet_truth_z_truth", jet_truth_z_truth, &b_jet_truth_z_truth);
   fChain->SetBranchAddress("jet_truth_index_z_reco", jet_truth_index_z_reco, &b_jet_truth_index_z_reco);
   fChain->SetBranchAddress("jet_truth_z_reco", jet_truth_z_reco, &b_jet_truth_z_reco);
   Notify();
}

Bool_t myClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
} 

void myClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t myClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef myClass_cxx
