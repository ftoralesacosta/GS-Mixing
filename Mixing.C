#define Mixing_cxx
#include "Mixing.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <valarray>

void Mixing::CentralityList(){
  
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

  vector <vector <pair <Long64_t,Float_t> > > mPreference, wPreference;  
  Float_t D;


  //was using 100 for a while
  int group_size = 1000;
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

  Long64_t poolsize = 20;
  
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

	  if ((men_pairings.at(suitor)).size() == poolsize) 
	    {
	      bachelors.pop_front(); 
	      cout<<endl<<"break"<<endl;
	      break; //haha it's fucking break!!!!! of course. 
	    }
	  
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


