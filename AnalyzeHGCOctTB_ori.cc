#define AnalyzeHGCOctTB_cxx
#include <iostream>
#include <vector>
#include <cstring>
#include "AnalyzeHGCOctTB.h"
#include "Math/SMatrix.h"
#include "Math/SVector.h"
#include <math.h>
#include<TF1.h>
using namespace std;



// chip 3022,44,3028




int main(int argc, char* argv[])//, int argvv[])
{

  if (argc < 5) {
    cerr << "Please give 5 arguments " << "runList " << " " << "outputFileName" << " " << "dataset" <<" " << "configuration" 
	 <<" " << "energy" << endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];
  const char *config          = argv[4];
  const char *energy = argv[5];
  AnalyzeHGCOctTB hgcOctTB(inputFileList, outFileName, data, config, energy);
  cout << "dataset " << data << " " << endl;
  cout << "configuration " << config << " " << endl;
  cout << "chi2 method " << energy << " " << endl;
 int chi2_method = atoi(energy);
  hgcOctTB.EventLoop(data,energy);
  return 0;
}

void AnalyzeHGCOctTB::EventLoop(const char *data, const char *energy) {


  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  //Long64_t nentries3 = fChain3->GetEntriesFast();
  Long64_t hgc_jentry = 0;

  cout << "nentries " << nentries << endl;
  cout << "Analyzing dataset " << data << " " << endl;

  Long64_t nbytes = 0, nb = 0;
  Long64_t nbytes2 = 0, nb2 = 0;
  Long64_t nbytes3 = 0, nb3 = 0;

  

  int decade = 0;
  int ahc_zeroHit = 0;



  bool DEBUG = false;
  //  bool DEBUG = true;
  int TOTAL_ACTIVE_LAYER = -1;
  int EE_LAYER = -1;
  int FH_LAYER = -1;
  int AH_LAYER = -1;
  if(!strcmp(conf_,"alpha") || !strcmp(conf_,"config1")) {
    TOTAL_ACTIVE_LAYER = 79;
    EE_LAYER = 28;
    FH_LAYER = 12;
    AH_LAYER = 39;
  }
  else if(!strcmp(conf_,"bravo") || !strcmp(conf_,"config2")){
    TOTAL_ACTIVE_LAYER = 78;
    EE_LAYER = 28;
    FH_LAYER = 11;
    AH_LAYER = 39;
  }
  else if(!strcmp(conf_,"charlie") || !strcmp(conf_,"config3")) {
    TOTAL_ACTIVE_LAYER = 59;
    EE_LAYER = 8;
    FH_LAYER = 12;
    AH_LAYER = 39;
    
  }
  else {
    cout<<"ERROR: Unknown configuration!!!!"<<endl;
    return;
  }
  int eve_count = 0;
  //counter
  int nHgc=0, nAhc=0, nrechits=0;
  float offset = 159.4; //A constant to be added to AHC layer z positions
  //cout<<energy<<endl;
  int chi2_method = atoi(energy);
  float FH_AH_relative_scale = 0.4;
  float alpha_ = FH_AH_relative_scale;
  float EE_scale = 94.624; //MIPs per Ge
  float FH_AH_scale = 12.788; //MIPs per GeV
 float ee_rescaling = 1.035;
  float fh_rescaling = 1.095;
  float ah_rescaling = 1.095;
  if(!strcmp(data,"data"))
    {
      ee_rescaling = 1;
      fh_rescaling = 1;
      ah_rescaling = 1;
    }
  cout<<"sim rescaling"<<"\t"<<ee_rescaling<<"\t"<<fh_rescaling<<"\t"<<ah_rescaling<<endl;
  
  std::ofstream file_ecal_ene_frac;
  file_ecal_ene_frac.open("Results/ecalEneFrac_true_pred_ss.txt",ios::out);

  int Elist[85] =  {10, 14 , 18 , 22 , 26 , 30,  34 , 38 , 42,  46,  50,  54,  58,  62,  66,  70,  74,  78,
                     82,  86,  90,  94,  98, 102, 106, 110, 114 ,118, 122, 126, 130, 134, 138, 142 ,146, 150,
                     154 ,158 ,162 ,166 ,170 ,174 ,178, 182, 186, 190 ,194 ,198, 202, 206, 210, 214 ,218, 222,
                  226, 230, 234 ,238, 242 ,246 ,250, 254, 258 ,262 ,266, 270, 274, 278, 282, 286, 290, 294,
                     298, 302 ,306, 310, 314, 318 ,322, 326, 330, 334 ,338 ,342 ,346};


  //to select the evnts within 2sigma region

  char* name1 = new char[1000];
  //reading chi2 weights
  sprintf(name1,"./txt_maps/updated_2022Maps/SSEE_2sigmaRange_trimAhcal_weights_21April22.txt");
  std::fstream file_EH_;
  file_EH_.open(name1,ios::in);
  sprintf(name1,"./txt_maps/updated_2022Maps/MipsEE_2sigmaRange_trimAhcal_weights_21April22.txt");
  std::fstream file_H_;
  file_H_.open(name1,ios::in);
  std::vector<map<int,double>> en_range_EE_xmin_vec;
  std::vector<map<int,double>> en_range_EE_xmax_vec;
  std::vector<map<int,double>> en_range_FH_xmin_vec;
  std::vector<map<int,double>> en_range_FH_xmax_vec;
  if(!file_H_.is_open())
    {
      std::cout << " file not opened" << std::endl;
    }
  else
    {
      int energy;
      double min_, max_;
      while (file_H_ >> energy >>min_>>max_) {
	// int energy;
	// float min_, max_;
	// file_H_ >> energy >>min_>>max_;
	std::map<int,double> en_range_FH_xmin;
	std::map<int,double> en_range_FH_xmax;
	en_range_FH_xmin[energy]=min_;
	en_range_FH_xmax[energy]=max_;
	en_range_FH_xmin_vec.push_back(en_range_FH_xmin);
	en_range_FH_xmax_vec.push_back(en_range_FH_xmax);
      }
    }
  if(!file_EH_.is_open())
    {
      std::cout << " file not opened" << std::endl;
    }
  else
    {
      int energy;
      double min_, max_;
      while (file_EH_ >>energy >>min_>>max_) {
	// int energy;
	// float min_, max_;
	// file_EH_ >>energy >>min_>>max_;
	std::map<int,double> en_range_FH_xmin;
	std::map<int,double> en_range_FH_xmax;
	en_range_FH_xmin[energy]=min_;
	en_range_FH_xmax[energy]=max_;
	en_range_EE_xmin_vec.push_back(en_range_FH_xmin);
	en_range_EE_xmax_vec.push_back(en_range_FH_xmax);
      }
    }

  // sprintf(name1,"./parametres_weights_v1.txt");
  // std::fstream file_p;
  // file_p.open(name1,ios::in);
  // std::vector<float> p0;
  // std::vector<float> p1;
  // std::vector<float> p2;
  // if(!file_p.is_open())
  //   {
  //     std::cout << " file not opened" << std::endl;
  //   }
  // else
  //   {
  //     float P0,P1,P2;
  //     while (file_p >>P0>>P1) {
  // 	p0.push_back(P0);
  // 	p1.push_back(P1);
  // 	// p2.push_back(P2);
  //     }
  //   }

  //fit function for 5M my fit parameters
//  TF1* f_EH_trimAhcal_w1 = new TF1("f_EH_trimAhcal_w1","sqrt([0]*[0]+[1]*[1]/x)+[2]/x+[3]", 5, 355);
//  f_EH_trimAhcal_w1->FixParameter(0,6.80016e-06);
//  f_EH_trimAhcal_w1->FixParameter(1,2.76128);
//  f_EH_trimAhcal_w1->FixParameter(2,2.03839);
//  f_EH_trimAhcal_w1->FixParameter(3,1.12455);
//   
//  TF1* f_EH_trimAhcal_w2 = new TF1("f_EH_trimAhcal_w2","sqrt([0]*[0]+[1]*[1]/x)+exp(-[2]*x+[3])", 5, 355);
//  f_EH_trimAhcal_w2->FixParameter(0,-2.04555e-06);
//  f_EH_trimAhcal_w2->FixParameter(1,7.60685);
//  f_EH_trimAhcal_w2->FixParameter(2,-4.12003);
//  f_EH_trimAhcal_w2->FixParameter(3,0.837049);
//  
//  TF1* f_EH_trimAhcal_w3= new TF1("f_EH_trimAhcal_w3","sqrt([0]*[0]+[1]*[1]/x)+[2]/x+[3]", 5, 355);
//  f_EH_trimAhcal_w3->FixParameter(0,4.03006e-07);
//  f_EH_trimAhcal_w3->FixParameter(1,0.839464);
//  f_EH_trimAhcal_w3->FixParameter(2,1.07043);
//  f_EH_trimAhcal_w3->FixParameter(3,0);
//
//  TF1* f_H_trimAhcal_w1 = new TF1("f_H_trimAhcal_w1","sqrt([0]*[0]+[1]/x)", 5, 355);
//  f_H_trimAhcal_w1->FixParameter(0,0);
//  f_H_trimAhcal_w1->FixParameter(1,0);
//  
//  TF1* f_H_trimAhcal_w2= new TF1("f_H_trimAhcal_w2","sqrt([0]*[0]+[1]/x)+[2]/x+[3]", 5, 355);
//  f_H_trimAhcal_w2->FixParameter(0,4.62966e-07);
//  f_H_trimAhcal_w2->FixParameter(1,1.38589);
//  f_H_trimAhcal_w2->FixParameter(2,0.414604);
//  f_H_trimAhcal_w2->FixParameter(3,0.787929);
//
//  TF1* f_H_trimAhcal_w3= new TF1("f_H_trimAhcal_w3","sqrt([0]*[0]+[1]/x)+[2]", 5, 355);
//  f_H_trimAhcal_w3->FixParameter(0,6.82505e-07);
//  f_H_trimAhcal_w3->FixParameter(1,1.79926);
//  f_H_trimAhcal_w3->FixParameter(2,0.978932);

// new updated 10M fit parameters
  TF1* f_EH_trimAhcal_w1 = new TF1("f_EH_trimAhcal_w1","sqrt([0]*[0]+[1]/x)+[2]/x+[3]", 5, 355);
  f_EH_trimAhcal_w1->FixParameter(0, -4.38488e-07 );
  f_EH_trimAhcal_w1->FixParameter(1, 4.06297);
  f_EH_trimAhcal_w1->FixParameter(2, 0.906989);
  f_EH_trimAhcal_w1->FixParameter(3, 1.124);

  TF1* f_EH_trimAhcal_w2 = new TF1("f_EH_trimAhcal_w2","sqrt([0]*[0]+[1]/x)+[2]/x+[3]", 5, 355);
  f_EH_trimAhcal_w2->FixParameter(0, 2.62734e-06);
  f_EH_trimAhcal_w2->FixParameter(1, 6.64859);
  f_EH_trimAhcal_w2->FixParameter(2, -3.83775);
  f_EH_trimAhcal_w2->FixParameter(3, 0.840613);

  TF1* f_EH_trimAhcal_w3= new TF1("f_EH_trimAhcal_w3","sqrt([0]*[0]+[1]/x)+[2]", 5, 355);
  f_EH_trimAhcal_w3->FixParameter(0, -2.23239e-06);
  f_EH_trimAhcal_w3->FixParameter(1, 0.596684);
  f_EH_trimAhcal_w3->FixParameter(2, 1.08795);
//  f_EH_trimAhcal_w3->FixParameter(3,0);

  TF1* f_H_trimAhcal_w1 = new TF1("f_H_trimAhcal_w1","sqrt([0]*[0]+[1]/x)", 5, 355);
  f_H_trimAhcal_w1->FixParameter(0,0);
  f_H_trimAhcal_w1->FixParameter(1,0);

  TF1* f_H_trimAhcal_w2= new TF1("f_H_trimAhcal_w2","sqrt([0]*[0]+[1]/x)+[2]/x+[3]", 5, 355);
  f_H_trimAhcal_w2->FixParameter(0, 8.73881e-07);
  f_H_trimAhcal_w2->FixParameter(1, 1.50499);
  f_H_trimAhcal_w2->FixParameter(2, -0.0352752);
  f_H_trimAhcal_w2->FixParameter(3, 0.790406);

  TF1* f_H_trimAhcal_w3= new TF1("f_H_trimAhcal_w3","sqrt([0]*[0]+[1]/x)+[2]", 5, 355);
  f_H_trimAhcal_w3->FixParameter(0, -1.40785e-07);
  f_H_trimAhcal_w3->FixParameter(1,  1.80381);
  f_H_trimAhcal_w3->FixParameter(2, 0.966127);


  // function for full ahcal scenarios
  TF1* f_EH_w1 = new TF1("f_EH_w1","sqrt([0]*[0]+[1]*[1]/x)+[2]/x+[3]", 5, 355);
  f_EH_w1->FixParameter(0,2.62110e-06);
  f_EH_w1->FixParameter(1,1.64811);
  f_EH_w1->FixParameter(2,1.18442);
  f_EH_w1->FixParameter(3,1.14035);

  TF1* f_EH_w2 = new TF1("f_EH_w2","sqrt([0]*[0]+[1]/x)+exp(-[2]*x+[3])", 5, 355);
  f_EH_w2->FixParameter(0,-3.15604e-06);//p0.at(4));
  f_EH_w2->FixParameter(1,1.27254e+00);//p1.at(4));
  f_EH_w2->FixParameter(2,1.55104e-04);//p2.at(4));
  f_EH_w2->FixParameter(3,-5.58071e-02);
   
  TF1* f_EH_w3 = new TF1("f_EH_w3","sqrt([0]*[0]+[1]*[1]/x)+[2]/x+[3]", 5, 355);
  f_EH_w3->FixParameter(0,-1.35815e-06);//p0.at(4));                                                                                       
  f_EH_w3->FixParameter(1,1.9127);//p1.at(4));                                                                                             
  f_EH_w3->FixParameter(2,-3.19441);//p2.at(4));                                                                                            
  f_EH_w3->FixParameter(3,1.01901);


  TF1* f_H_w1 = new TF1("f_H_w1","sqrt([0]*[0]+[1]/x)", 5, 355);
  f_H_w1->FixParameter(0,0);//p0.at());
  f_H_w1->FixParameter(1,0);//.at(6));
  TF1* f_H_w2= new TF1("f_H_w2","sqrt([0]*[0]+[1]/x)+[2]/x", 5, 355);
  f_H_w2->FixParameter(0,8.47109e-01);
  f_H_w2->FixParameter(1,2.38987e+01);
  f_H_w2->FixParameter(2,-6.49213e+00);
 

  TF1* f_H_w3= new TF1("f_H_w3","sqrt([0]*[0]+[1]*[1]/x)+[2]/x+[3]", 5, 355);
  f_H_w3->FixParameter(0,-6.04190e-07);
  f_H_w3->FixParameter(1,2.16114e+00);
  f_H_w3->FixParameter(2,-3.15370e+00);
  f_H_w3->FixParameter(3,9.17243e-01);

  TF1* fnct_EH_w1 = new TF1("fnct_EH_w1","sqrt([0]*[0]+[1]*[1]/x)", 5, 355);
  fnct_EH_w1->FixParameter(0,1.189897);
  fnct_EH_w1->FixParameter(1,5.215187);
      
  TF1* fnct_EH_w2 = new TF1("fnct_EH_w2","sqrt([0]*[0]+[1]*[1]/x)", 5, 355);
  fnct_EH_w2->FixParameter(0,1.004924);
  fnct_EH_w2->FixParameter(1,3.352913);
  TF1* fnct_EH_w3 = new TF1("fnct_EH_w3","sqrt([0]*[0]+[1]*[1]/x)", 5,355);
  fnct_EH_w3->FixParameter(0,1.025529);
  fnct_EH_w3->FixParameter(1,3.602263);
      
  TF1* fnct_H_w1 = new TF1("fnct_H_w1","sqrt([0]*[0]+[1]*[1]/x)", 5, 355);
  fnct_H_w1->FixParameter(0,0.0);
  fnct_H_w1->FixParameter(1,0.0);
  TF1* fnct_H_w2 = new TF1("fnct_H_w2","sqrt([0]*[0]+[1]*[1]/x)", 5, 355);
  fnct_H_w2->FixParameter(0,0.878231);
  fnct_H_w2->FixParameter(1,2.659706);
  TF1* fnct_H_w3 = new TF1("fnct_H_w3","sqrt([0]*[0]+[1]*[1]/x)", 5, 355);
  fnct_H_w3->FixParameter(0,0.984368);
  fnct_H_w3->FixParameter(1,2.765905);
    
  std::vector<ROOT::Math::SVector<double, 3> >consts_EH_vec;
  std::vector<ROOT::Math::SVector<double, 3> >values_EH_vec;
  std::vector<ROOT::Math::SMatrix<double,3, 3, ROOT::Math::MatRepStd<double,3>> > coeff_EH_vec;
  ROOT::Math::SMatrix<double,3, 3, ROOT::Math::MatRepStd<double,3> > coeffs_EH;
  ROOT::Math::SVector<double, 3> consts_EH;
  ROOT::Math::SVector<double, 3> values_EH;

  bool isInverted_EH = false;
  for(int i =0; i<85; i++)
    {
      for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
          coeffs_EH(i,j) = 0.0;
        }
        consts_EH(i) = 0.0;
        values_EH(i) = 0.0;
      }
        consts_EH_vec.push_back(consts_EH);
        values_EH_vec.push_back(values_EH);
        coeff_EH_vec.push_back(coeffs_EH);
    }
  std::vector<ROOT::Math::SVector<double, 2> >consts_H_vec;
  std::vector<ROOT::Math::SVector<double, 2> >values_H_vec;
  std::vector<ROOT::Math::SMatrix<double,2, 2, ROOT::Math::MatRepStd<double,2>> > coeff_H_vec;
  bool isInverted_H = false;
  ROOT::Math::SMatrix<double,2, 2, ROOT::Math::MatRepStd<double,2> > coeffs_H;
  ROOT::Math::SVector<double, 2> consts_H;
  ROOT::Math::SVector<double, 2> values_H;

  for(int i =0; i<85; i++)
    {
      for(int i = 0; i < 2; i++) {
        for(int j = 0; j < 2; j++) {
          coeffs_H(i,j) = 0.0;
        }
        consts_H(i) = 0.0;
        values_H(i) = 0.0;
      }
      consts_H_vec.push_back(consts_H);
      values_H_vec.push_back(values_H);
      coeff_H_vec.push_back(coeffs_H);
    }

  //  cout<<"size"<<"\t"<<consts_H_vec.size()<<endl;
  
  if(DEBUG) cout<<"DEBUG: Configuration = "<<conf_<<endl;
  if(DEBUG) cout<<"DEBUG: TOTAL_ACTIVE_LAYER = "<<TOTAL_ACTIVE_LAYER<<endl;
  if(DEBUG) cout<<"DEBUG: EE_LAYER = "<<EE_LAYER<<endl;
  if(DEBUG) cout<<"DEBUG: FH_LAYER = "<<FH_LAYER<<endl;
  if(DEBUG) cout<<"DEBUG: AH_LAYER = "<<AH_LAYER<<endl;

  if(DEBUG) cout<<"DEBUG: Entering event Loop"<<endl;

  Long64_t jentry;


  float lambda[79];
  for(int i = 0; i < 79; i++) {
    lambda[i] = layer_positions[i+1].at(2); // for nuclear interaction length  & //pi_lambda use -at(3)
    //     cout<<lambda[i]<<endl;
  }

  // float lambda_pi[79];
  // for(int i = 0; i < 79; i++) {
    
  //   lambda_pi[i] = layer_positions[i+1].at(3); // for nuclear interaction length  & //pi_lambda use -at(3)
  //   cout<<lambda_pi[i]<<endl;
  // }
  int ahcal_layer[10]= {43,47,51,55,59,63,67,71,75,79};
  //nentries =1000;
  for (jentry=0; jentry<nentries;jentry++,hgc_jentry++)
  {
    //   for (jentry=0; jentry<10000;jentry++,hgc_jentry++) {
    // ==============print number of events done == == == == == == == =
      double progress = 10.0 * jentry / (1.0 * nentries);
      int k = int (progress);
      if (k > decade)
	cout << 10 * k << " %" << endl;
      decade = k;
    
      // ===============read this entry == == == == == == == == == == ==

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) { break; cout<<"Breaking"<<endl;}
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      //cout<<"insides the event loop: check1"<<endl;

      ////   MESSAGE ////
      // apparently code is not running beyond this point ///
      event_count[0]++;
      
      
      event_count[1]++;
      h_true_beamenergy[4]->Fill(trueBeamEnergy);
      h_particle->Fill(pdgID);
      totalEnergy_inGeV=0;
      EnergySum_SSinEE=0;
      EnergySum_SSinFH=0;
      Esum_rechits_FH=0;
      Esum_rechits_EE=0;
      Esum_rechits_AH=0;
      Esum_rechits_FH_inGeV=0;
      Esum_rechits_EE_inGeV=0;
      Esum_rechits_AH_inGeV=0;
      
      ////////////////////////////////////////////
      //            HGCAL Part                  //
      ////////////////////////////////////////////
      eve_count+=1;
      total_rechits=0;
      rechits_EE=0;
      rechits_FH=0;
      rechits_AH=0;
      // if(event ==3482)
      // 	{	//	continue;
      int nrechit_trimAhcal=0;
      /// Read HGCAL Tree
      if(DEBUG) cout<<"DEBUG: Start Analylizing HGCAL  RecHits!!"<<endl;
      // total_rechits_energy_EE=0;
      // total_rechits_energy_FH=0;
      //      cout<<NRechits<<"\t"<<rechitEn_trimAhcal->size()<<endl;
      //loop over combined rechits for trim AHCAL 
      //float rechitEnergySum_AH=0.0;
      for( int i = 0 ; i < rechitEn_trimAhcal->size(); i++)
      	{
      	  //int layer = rechit_layer->at(i);
      	  float energy=rechitEn_trimAhcal->at(i);
	  //cout<<i<<"\t"<<layer<<"\t"<<energy<<endl;
	  
	  if(comb_rechit_z_trimAhcal->at(i)<54)
	    {
	      Esum_rechits_EE+=(energy/ee_rescaling);
	      rechits_EE++;
	    }
	  else if((comb_rechit_z_trimAhcal->at(i)>54) && (comb_rechit_z_trimAhcal->at(i)<154))
	    {
	      
	      Esum_rechits_FH+=(energy/fh_rescaling);
	    }
	  else if(comb_rechit_z_trimAhcal->at(i)>154)
	    {
	      Esum_rechits_AH+=(energy/ah_rescaling);
	    }
  	  } //Nrechits loop
      //cout<<NRechits<<"\t"<<count<<endl;

      h_true_beamenergy[5]->Fill(trueBeamEnergy);
      nrechit_trimAhcal=0;
      rechits_AH = nAhc;
      total_rechits = nrechits;

      Esum_rechits_EE_inGeV=0.0105*Esum_rechits_EE;
      Esum_rechits_FH_inGeV= 0.0812*Esum_rechits_FH; //updated relative weight
      Esum_rechits_AH_inGeV=0.12504*Esum_rechits_AH;
      double rechitEnergySum_EE = Esum_rechits_EE;
      double rechitEnergySum_FH = Esum_rechits_FH;
      double rechitEnergySum_AH = Esum_rechits_AH;
      double EE_detscale = Esum_rechits_EE_inGeV;//(0.0105rechitEnergySum_EE/EE_scale);                                        
      double FH_detscale = Esum_rechits_FH_inGeV;//(rechitEnergySum_FH/FH_AH_scale);                                                 
      double AH_detscale = Esum_rechits_AH_inGeV;//(alpha_*rechitEnergySum_AH)/FH_AH_scale;                                   
      double full_energy = FH_detscale+AH_detscale;                                                    
      double total_energy= EE_detscale+ FH_detscale+ AH_detscale;
      double ecal_ene_frac = EE_detscale/total_energy;
      
      //file_ecal_ene_frac<< ecal_ene_frac <<"\n";

      ///////////////////////////////////////////////////////////                                                                           
      /////     H hadrons ; Chi2 matrix initialzation     /////                                                                             
      //////////////////////////////////////////////////////////
      //cout<<"alp_check"<<endl;
      double E_beam =0.0;
      double O = 0.0;
      double sigma2_H = (0.09*0.09)*(total_energy*total_energy) + (1.23*1.23)*total_energy;
      double tot_E_gev= total_energy;
      float chi2_funct=0.0;
      for(int i_bin=0;i_bin<85;i_bin++)
	{
	  bool IsIn_en_range = false;
          if(trueBeamEnergy>=Elist[i_bin] && trueBeamEnergy <=Elist[i_bin]+4)
            {
	      //cout<<i_bin<<"\t"<<Elist[i_bin]<<"\t"<<rechit_shower_start_layer<<endl;
	      E_beam = Elist[i_bin]+2;
	      if(rechit_shower_start_layer>28)
		{
		  IsIn_en_range = (total_energy > en_range_FH_xmin_vec[i_bin][Elist[i_bin]+2] && total_energy < en_range_FH_xmax_vec[i_bin][Elist[i_bin]+2]); 
		  if(IsIn_en_range){
		  //cout<<i_bin<<"\t"<<Elist[i_bin]<<endl;

		  ROOT::Math::SMatrix<double,2, 2, ROOT::Math::MatRepStd<double,2> > coeffs_H;
		  ROOT::Math::SVector<double, 2> consts_H;
		  ROOT::Math::SVector<double, 2> values_H;
		  for(int i = 0; i < 2; i++) {
		    for(int j = 0; j < 2; j++) {
		      coeffs_H(i,j) = 0.0;
		    }	
		    consts_H(i) = 0.0;
		    values_H(i) = 0.0;
		  }
		  //cout<<i_bin<<"\t"<<Elist[i_bin]<<endl;

		  
		  if(chi2_method == 0) // method -input to chi2 are in units of MIPs
		    {
		      coeffs_H(0,0) = (rechitEnergySum_FH*rechitEnergySum_FH)/sigma2_H;
		      coeffs_H(0,1) = (rechitEnergySum_AH*rechitEnergySum_FH)/sigma2_H;
		      coeffs_H(1,0) = (rechitEnergySum_FH*rechitEnergySum_AH)/sigma2_H;
		      coeffs_H(1,1) = (rechitEnergySum_AH*rechitEnergySum_AH)/sigma2_H;
		      
		      consts_H(0)   = (E_beam-O)*rechitEnergySum_FH/sigma2_H;
		      consts_H(1)   = (E_beam-O)*rechitEnergySum_AH/sigma2_H;
		      consts_H_vec[i_bin]+= consts_H;
		      coeff_H_vec[i_bin]+=coeffs_H;
		    }
		  else // method -input to chi2 are in units of GeV -> We use this one
		    {
		      coeffs_H(0,0) = (FH_detscale*FH_detscale)/sigma2_H;
		      coeffs_H(0,1) = (AH_detscale*FH_detscale)/sigma2_H;
		      coeffs_H(1,0) = (FH_detscale*AH_detscale)/sigma2_H;
		      coeffs_H(1,1) = (AH_detscale*AH_detscale)/sigma2_H;
		      
		      consts_H(0)   = (E_beam-O)*FH_detscale/sigma2_H;
		      consts_H(1)   = (E_beam-O)*AH_detscale/sigma2_H;
		      //cout<<i_bin<<"\t"<<Elist[i_bin]<<"\t"<<"beforE"<<endl;

		      consts_H_vec[i_bin] += consts_H;
		      coeff_H_vec[i_bin]+=coeffs_H;

		    }
		  }
		}
	      else
		{ // for EH hadrons
		  //cout<<jentry<<"\t"<<"beforE"<<endl;
		    IsIn_en_range = (tot_E_gev > en_range_EE_xmin_vec[i_bin][Elist[i_bin]+2] && tot_E_gev < en_range_EE_xmax_vec[i_bin][Elist[i_bin]+2]);
		    if(IsIn_en_range){
		  double sigma2_EH = (0.084*0.084)*tot_E_gev*tot_E_gev + (1.39*1.39)*tot_E_gev;
		  ROOT::Math::SMatrix<double,3, 3, ROOT::Math::MatRepStd<double,3> > coeffs_EH;
		  ROOT::Math::SVector<double, 3> consts_EH;
		  ROOT::Math::SVector<double, 3> values_EH;

		  for(int i = 0; i < 3; i++) {
		    for(int j = 0; j < 3; j++) {
		      coeffs_EH(i,j) = 0.0;
		    }
		    consts_EH(i) = 0.0;
		    values_EH(i) = 0.0;
		  }
		  

		  if(chi2_method==0)
		    {
		      coeffs_EH(0,0) = (rechitEnergySum_EE*rechitEnergySum_EE)/sigma2_EH;
		      coeffs_EH(0,1) = (rechitEnergySum_EE*rechitEnergySum_FH)/sigma2_EH;
		      coeffs_EH(0,2) = (rechitEnergySum_EE*rechitEnergySum_AH)/sigma2_EH;
		      coeffs_EH(1,0) = (rechitEnergySum_FH*rechitEnergySum_EE)/sigma2_EH;
		      coeffs_EH(1,1) = (rechitEnergySum_FH*rechitEnergySum_FH)/sigma2_EH;
		      coeffs_EH(1,2) = (rechitEnergySum_FH*rechitEnergySum_AH)/sigma2_EH;
		      coeffs_EH(2,0) = (rechitEnergySum_AH*rechitEnergySum_EE)/sigma2_EH;
		      coeffs_EH(2,1) = (rechitEnergySum_AH*rechitEnergySum_FH)/sigma2_EH;
		      coeffs_EH(2,2) = (rechitEnergySum_AH*rechitEnergySum_AH)/sigma2_EH;

		      consts_EH(0)   = (E_beam-O)*rechitEnergySum_EE/sigma2_EH;
		      consts_EH(1)   = (E_beam-O)*rechitEnergySum_FH/sigma2_EH;
		      consts_EH(2)   = (E_beam-O)*rechitEnergySum_AH/sigma2_EH;
		      consts_EH_vec[i_bin]+= consts_EH;
		      coeff_EH_vec[i_bin]+=coeffs_EH;

		    }
		  else
		    {
		      coeffs_EH(0,0) = (EE_detscale*EE_detscale)/sigma2_EH;
		      coeffs_EH(0,1) = (EE_detscale*FH_detscale)/sigma2_EH;
		      coeffs_EH(0,2) = (EE_detscale*AH_detscale)/sigma2_EH;
		      coeffs_EH(1,0) = (FH_detscale*EE_detscale)/sigma2_EH;
		      coeffs_EH(1,1) = (FH_detscale*FH_detscale)/sigma2_EH;
		      coeffs_EH(1,2) = (FH_detscale*AH_detscale)/sigma2_EH;
		      coeffs_EH(2,0) = (AH_detscale*EE_detscale)/sigma2_EH;
		      coeffs_EH(2,1) = (AH_detscale*FH_detscale)/sigma2_EH;
		      coeffs_EH(2,2) = (AH_detscale*AH_detscale)/sigma2_EH;


		      consts_EH(0)   = (E_beam - O)*EE_detscale/sigma2_EH;
		      consts_EH(1)   = (E_beam - O)*FH_detscale/sigma2_EH;
		      consts_EH(2)   = (E_beam - O)*AH_detscale/sigma2_EH;
		      consts_EH_vec[i_bin]+= consts_EH;
		      coeff_EH_vec[i_bin]+=coeffs_EH;
		      
		    }
		  }
		}
	      //applying chi2 weights
	      float trimAhcal_w1, trimAhcal_w2, trimAhcal_w3;
	      if(rechit_shower_start_layer>28)
		{
		  //float w1=0,w2=0,w3=0;
		  if (!strcmp(data,"data"))                    
		    trueBeamEnergy = beamEnergy;

		   trimAhcal_w1 = f_H_trimAhcal_w1->Eval(trueBeamEnergy);
		  trimAhcal_w2 = f_H_trimAhcal_w2->Eval(trueBeamEnergy);
		  trimAhcal_w3 = f_H_trimAhcal_w3->Eval(trueBeamEnergy);
		  chi2_funct = EE_detscale +trimAhcal_w2*FH_detscale + trimAhcal_w3*AH_detscale;
		  hist_resp_total_trimAhcal[i_bin]->Fill(chi2_funct);
		  hist_resp_SS_FH_trimAhcal[i_bin]->Fill(chi2_funct);
		}
	      else
		{
		  if (!strcmp(data,"data"))
                  
		    trueBeamEnergy = beamEnergy;

                  
                  trimAhcal_w1 = f_EH_trimAhcal_w1->Eval(trueBeamEnergy);
                  trimAhcal_w2 = f_EH_trimAhcal_w2->Eval(trueBeamEnergy);
                  trimAhcal_w3 = f_EH_trimAhcal_w3->Eval(trueBeamEnergy);
		  chi2_funct = trimAhcal_w1*EE_detscale +trimAhcal_w2*FH_detscale + trimAhcal_w3*AH_detscale;
                  hist_resp_total_trimAhcal[i_bin]->Fill(chi2_funct);
                  hist_resp_SS_EE_trimAhcal[i_bin]->Fill(chi2_funct);
		  
		}
            file_ecal_ene_frac<< ecal_ene_frac <<"\t"<< trueBeamEnergy <<"\t"<< chi2_funct << "\t" << rechit_shower_start_layer << "\n";// << trimAhcal_w1 <<"\t"<< trimAhcal_w2  <<"\t"<< trimAhcal_w3  <<"\n";
            }

	  
        }//true bin loop   


 
      if(DEBUG) cout<<"DEBUG: End of Entry = "<<jentry<<endl;
      if(DEBUG) cout<<"DEBUG: ****************** "<<endl;
      //      cout<<"\t"<<"beforE"<<endl;

    } // loop over entries
  
  char* name = new char[1000];
  sprintf(name,"Results/chi2_calibFact_H_hadrons_flatEn_%d_v1.txt",chi2_method);
  std::ofstream file_H;
  file_H.open(name,ios::out);//"../txt_maps/chi2_map/chi2_calibFact_H_hadrons_flatEn_%d.txt",chi2_method,ios::out);                        
  sprintf(name,"Results/chi2_calibFact_EH_hadrons_flatEn_%d_v1.txt",chi2_method);
  std::ofstream file_EH;
  file_EH.open(name,ios::out);//"../txt_maps/chi2_map/chi2_calibFact_EH_hadrons_flatEn_%d.txt",chi2_method,ios::out

  /////////////////////////////////////////////////////                                                                                     
  //////                 EH Hadrons            ////////                                                                                     
  /////////////////////////////////////////////////////                                                                                     

   for(int i_en =0; i_en<85; i_en++)
    {
      ROOT::Math::SMatrix<double,3, 3, ROOT::Math::MatRepStd<double,3> > coeffs_EH;
      ROOT::Math::SVector<double, 3> consts_EH;
      ROOT::Math::SVector<double, 3> values_EH;
      ROOT::Math::SMatrix<double,2, 2, ROOT::Math::MatRepStd<double,2> > coeffs_H;
      ROOT::Math::SVector<double, 2> consts_H;
      ROOT::Math::SVector<double, 2> values_H;
       //      bool isInverted_EH = false;                                                                                                  
      for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
          coeffs_EH(i,j) = 0.0;
        }
        consts_EH(i) = 0.0;
        values_EH(i) = 0.0;
      }
       for(int i = 0; i < 2; i++) {
        for(int j = 0; j < 2; j++) {
          coeffs_H(i,j) = 0.0;
        }
        consts_H(i) = 0.0;
        values_H(i) = 0.0;
      }

      consts_EH = consts_EH_vec[i_en];
      values_EH = values_EH_vec[i_en];
      coeffs_EH = coeff_EH_vec[i_en];
      //cout<<"before invert"<<"\t"<<coeffs_EH<<endl;

      isInverted_EH = coeffs_EH.Invert();// cout<<"before invert"<<"\t"<<coeffs_EH<<endl;
      for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
          cout<<coeffs_EH(i, j)<<"\t";
        }
        cout<<endl;
      }

      cout<<" "<<consts_EH(0)<<"\t"<<consts_EH(1)<<"\t"<<consts_EH(2)<<endl;
consts_H = consts_H_vec[i_en];
      values_H = values_H_vec[i_en];
      coeffs_H = coeff_H_vec[i_en];

      isInverted_H = coeffs_H.Invert();
      for(int i = 0; i <2; i++) {
        for(int j = 0; j < 2; j++) {
          cout<<coeffs_H(i, j)<<"\t";
        }
        cout<<endl;
      }
      cout<<" "<<consts_H(0)<<"\t"<<consts_H(1)<<"\t"<<consts_H(2)<<endl;
      if(isInverted_EH) {
        values_EH = coeffs_EH*consts_EH;
	cout<<"EH Hadrons => For E = "<<Elist[i_en]+2<<"GeV, w1 = "<<values_EH(0)<<" ;w2 = "<<values_EH(1)<<" ;w3 = "<<values_EH(2)<<";"<<endl;
	cout<<">>>>>>> Uncertainity  : "<<sqrt(coeffs_EH(0,0))<< " " << sqrt(coeffs_EH(1,1)) << " " << sqrt(coeffs_EH(2,2))<<endl;

	file_EH<<Elist[i_en]+2<<"\t"<<values_EH(0)<<"\t"<<values_EH(1)<<"\t"<<values_EH(2)<<"\t"<<sqrt(coeffs_EH(0,0))<<"\t"<<sqrt(coeffs_EH(1,1))<<"\t"<<sqrt(coeffs_EH(2,2))<<"\n";
        //cout<<"EH Hadrons => For E = "<<Elist[i_en]+2<<"GeV, w1 = "<<values_EH(0)<<" ;w2 = "<<values_EH(1)<<" ;w3 = "<<values_EH(2)<<";"<<endl;
        //file_EH<<Elist[i_en]+2<<"\t"<<values_EH(0)<<"\t"<<values_EH(1)<<"\t"<<values_EH(2)<<"\n";
      }
      else {
        cout<<"Error: Could not invert for EH hadrons..."<<endl;
      }
      //isInverted_H = true;                                                                                                                
      if(isInverted_H) {
        values_H = coeffs_H*consts_H;
        //cout<<"H Hadrons => For E = "<<E_beam<<"GeV, w1 = "<<values_H(0)<<" ;w2 = "<<values_H(1)<<" ;w3 = "<<values_H(2)<<";"<<endl;   

	cout<<"H Hadrons => For E = "<<Elist[i_en]+2<<"GeV, w1 = 0.0 ;w2 = "<<values_H(0)<<" ;w3 = "<<values_H(1)<<";"<<endl;
	 cout<<">>>>>>> Uncertainity  : "<<sqrt(coeffs_H(0,0)) << " " <<sqrt(coeffs_H(1,1))<<endl;

	 file_H<<Elist[i_en]+2<<"\t"<<0<<"\t"<<values_H(0)<<"\t"<<values_H(1)<<"\t"<<0<<"\t"<<sqrt(coeffs_H(0,0))<<"\t"<<sqrt(coeffs_H(1,1))<<"\n";   
	 //cout<<"H Hadrons => For E = "<<Elist[i_en]+2<<"GeV, w1 = 0.0 ;w2 = "<<values_H(0)<<" ;w3 = "<<values_H(1)<<";"<<endl;
	 //file_H<<Elist[i_en]+2<<"\t"<<values_H(0)<<"\t"<<values_H(1)<<"\n";
      }
      else {
        cout<<"Error: Could not invert for H hadrons..."<<endl;
      }
    }


  cout<<endl<<endl;


  ///////////////////////////////////////////////////////////
  ///////  E N D     O F     E N T R Y     L O O P     //////
  ///////////////////////////////////////////////////////////
  //  gSystem->Exit(0);
  cout<<"Got Out "<<jentry<<endl;
  cout<<count_fh<<endl;
  cout<<count_badtrack<<endl;
  cout<<"eve_count="<<eve_count<<endl;
  for(int i =0;i<7;i++)
    {
      cout<<event_count[i]<<endl;
     }
}
