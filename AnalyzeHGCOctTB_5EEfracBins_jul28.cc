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

  char* path_to_save = new char[1000];
  sprintf(path_to_save, "Results/");

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
  file_ecal_ene_frac.open("Results/ecalEneFrac_true_pred_SS.txt",ios::out);  

  int Elist[85] =  {10, 14 , 18 , 22 , 26 , 30,  34 , 38 , 42,  46,  50,  54,  58,  62,  66,  70,  74,  78,
                     82,  86,  90,  94,  98, 102, 106, 110, 114 ,118, 122, 126, 130, 134, 138, 142 ,146, 150,
                     154 ,158 ,162 ,166 ,170 ,174 ,178, 182, 186, 190 ,194 ,198, 202, 206, 210, 214 ,218, 222,
                  226, 230, 234 ,238, 242 ,246 ,250, 254, 258 ,262 ,266, 270, 274, 278, 282, 286, 290, 294,
                     298, 302 ,306, 310, 314, 318 ,322, 326, 330, 334 ,338 ,342 ,346};
  
  float ecal_ene_frac_bins[6] = {0.0, 0.06, 0.2, 0.4, 0.6, 1.01};

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



  //fit function

  TF1* f_w1_bin1 = new TF1("f_w1_bin1","sqrt([0]*[0]+[1]/x)+[2]/x+[3]",5,360);
  f_w1_bin1->FixParameter(0,0);
  f_w1_bin1->FixParameter(1,0);
  f_w1_bin1->FixParameter(2,0);
  f_w1_bin1->FixParameter(3,0);

  TF1* f_w2_bin1 = new TF1("f_w2_bin1","sqrt([0]*[0]+[1]/x)+[2]/x+[3]",5,360);
  f_w2_bin1->FixParameter(0, 0.1211832);
  f_w2_bin1->FixParameter(1, 3.3331938);
  f_w2_bin1->FixParameter(2, -2.275214);
  f_w2_bin1->FixParameter(3, 0.7325118);

  TF1* f_w3_bin1 = new TF1("f_w3_bin1","sqrt([0]*[0]+[1]/x)+[2]/x+[3]",5,360);
  f_w3_bin1->FixParameter(0, 1.0794749);
  f_w3_bin1->FixParameter(1, 6.2016315);
  f_w3_bin1->FixParameter(2, -0.0002443);
  f_w3_bin1->FixParameter(3, 0.01856074);


  TF1* f_w1_bin2 = new TF1("f_w1_bin2","sqrt([0]*[0]+[1]/x)+[2]",5,360);
  f_w1_bin2->FixParameter(0, 0.2372573);
  f_w1_bin2->FixParameter(1, 13.9916573);
  f_w1_bin2->FixParameter(2, 1.1076621);

  TF1* f_w2_bin2 = new TF1("f_w2_bin2","sqrt([0]*[0]+[1]/(x*x))+[2]*x*x+[3]/x",5,360);
  f_w2_bin2->FixParameter(0, 0.94995797);
  f_w2_bin2->FixParameter(1, 57.88050471);
  f_w2_bin2->FixParameter(2, -4.00533886e-07);
  f_w2_bin2->FixParameter(3, 3.849738014);

  TF1* f_w3_bin2 = new TF1("f_w3_bin2","sqrt([0]*[0]+[1]/x)+[2]",5,360);
  f_w3_bin2->FixParameter(0, 1.1707796);
  f_w3_bin2->FixParameter(1, 5.08021159);
  f_w3_bin2->FixParameter(2, -0.0062870);


  TF1* f_w1_bin3 = new TF1("f_w1_bin3","sqrt([0]*[0]+[1]/x)+[2]",5,360);
  f_w1_bin3->FixParameter(0, 1.897244578e-06);
  f_w1_bin3->FixParameter(1, 9.783924263);
  f_w1_bin3->FixParameter(2, 1.162540391);

  TF1* f_w2_bin3 = new TF1("f_w2_bin3","sqrt([0]*[0]+[1]/x)+[2]/x+[3]",5,360);
  f_w2_bin3->FixParameter(0, 0.069849295);
  f_w2_bin3->FixParameter(1, 3.899134583);
  f_w2_bin3->FixParameter(2, -3.07192980);
  f_w2_bin3->FixParameter(3, 0.842306012);


  TF1* f_w3_bin3 = new TF1("f_w3_bin3","sqrt([0]*[0]+[1]/x)+[2]/x+[3]",5,360);
  f_w3_bin3->FixParameter(0, 1.108526756);
  f_w3_bin3->FixParameter(1, 28.71782046);
  f_w3_bin3->FixParameter(2, -6.576580516);
  f_w3_bin3->FixParameter(3, -0.015224433);


  TF1* f_w1_bin4 = new TF1("f_w1_bin4","sqrt([0]*[0]+[1]/x)+[2]/x+[3]",5,360);
  f_w1_bin4->FixParameter(0, 0.066766587);
  f_w1_bin4->FixParameter(1, 26.56138793);
  f_w1_bin4->FixParameter(2, -6.799449467);
  f_w1_bin4->FixParameter(3, 0.9479815981);

  TF1* f_w2_bin4 = new TF1("f_w2_bin4","sqrt([0]*[0]+[1]/x)+[2]/x+[3]",5,360);
  f_w2_bin4->FixParameter(0, 1.293058311e-05);
  f_w2_bin4->FixParameter(1, 0.2152833864);
  f_w2_bin4->FixParameter(2, 1.4448171546);
  f_w2_bin4->FixParameter(3, 1.0239101593);

  TF1* f_w3_bin4 = new TF1("f_w3_bin4","sqrt([0]*[0]+[1]/x)+[2]/x+[3]",5,360);
  f_w3_bin4->FixParameter(0, 1.10204086);
  f_w3_bin4->FixParameter(1, 0.01203362);
  f_w3_bin4->FixParameter(2, 4.66728660);
  f_w3_bin4->FixParameter(3, -0.01402613);


  TF1* f_w1_bin5 = new TF1("f_w1_bin5","sqrt([0]*[0]+[1]/x)+[2]/x+[3]",5,360);
  f_w1_bin5->FixParameter(0, 0.154847123);
  f_w1_bin5->FixParameter(1, 5.133857756);
  f_w1_bin5->FixParameter(2, 2.081803157);
  f_w1_bin5->FixParameter(3, 0.849460533);

  TF1* f_w2_bin5 = new TF1("f_w2_bin5","sqrt([0]*[0]+[1]/x)+[2]/x+[3]",5,360);
  f_w2_bin5->FixParameter(0, 8.791864739e-06);
  f_w2_bin5->FixParameter(1, 72.41649617);
  f_w2_bin5->FixParameter(2, -19.96473945);
  f_w2_bin5->FixParameter(3, 0.93810884218);

  TF1* f_w3_bin5 = new TF1("f_w3_bin5","sqrt([0]*[0]+[1]/x)+[2]/x+[3]",5,360);
  f_w3_bin5->FixParameter(0, 0.087574234);
  f_w3_bin5->FixParameter(1, 39.35367556);
  f_w3_bin5->FixParameter(2, -14.35196768);
  f_w3_bin5->FixParameter(3, 0.9863167649);





 
  std::vector<ROOT::Math::SVector<double, 3> >consts_high_ecal_vec;
  std::vector<ROOT::Math::SVector<double, 3> >values_high_ecal_vec;
  std::vector<ROOT::Math::SMatrix<double,3, 3, ROOT::Math::MatRepStd<double,3>> > coeff_high_ecal_vec;
  ROOT::Math::SMatrix<double,3, 3, ROOT::Math::MatRepStd<double,3> > coeffs_high_ecal;
  ROOT::Math::SVector<double, 3> consts_high_ecal;
  ROOT::Math::SVector<double, 3> values_high_ecal;

  std::vector<std::vector<ROOT::Math::SVector<double, 3> > > consts_high_ecal_vec_diff_bins;
  std::vector<std::vector<ROOT::Math::SVector<double, 3> > > values_high_ecal_vec_diff_bins;
  std::vector<std::vector<ROOT::Math::SMatrix<double,3, 3, ROOT::Math::MatRepStd<double,3>> > >coeff_high_ecal_vec_diff_bins;

  
  for (int ecal_bin=0; ecal_bin<4; ecal_bin++){
    for(int i =0; i<85; i++)
        {
        for(int i = 0; i < 3; i++) {
            for(int j = 0; j < 3; j++) {
            coeffs_high_ecal(i,j) = 0.0;
            }
            consts_high_ecal(i) = 0.0;
            values_high_ecal(i) = 0.0;
        }
            consts_high_ecal_vec.push_back(consts_high_ecal);
            values_high_ecal_vec.push_back(values_high_ecal);
            coeff_high_ecal_vec.push_back(coeffs_high_ecal);
        }

        consts_high_ecal_vec_diff_bins.push_back(consts_high_ecal_vec);
        values_high_ecal_vec_diff_bins.push_back(values_high_ecal_vec);
        coeff_high_ecal_vec_diff_bins.push_back(coeff_high_ecal_vec);
//        cout<<ecal_ene_frac_bins[ecal_bin]<<endl;
    }
  std::vector<ROOT::Math::SVector<double, 2> >consts_low_ecal_vec;
  std::vector<ROOT::Math::SVector<double, 2> >values_low_ecal_vec;
  std::vector<ROOT::Math::SMatrix<double,2, 2, ROOT::Math::MatRepStd<double,2>> > coeff_low_ecal_vec;
  bool isInverted = false;
  ROOT::Math::SMatrix<double,2, 2, ROOT::Math::MatRepStd<double,2> > coeffs_low_ecal;
  ROOT::Math::SVector<double, 2> consts_low_ecal;
  ROOT::Math::SVector<double, 2> values_low_ecal;

  for(int i =0; i<85; i++)
    {
      for(int i = 0; i < 2; i++) {
        for(int j = 0; j < 2; j++) {
          coeffs_low_ecal(i,j) = 0.0;
        }
        consts_low_ecal(i) = 0.0;
        values_low_ecal(i) = 0.0;
      }
      consts_low_ecal_vec.push_back(consts_low_ecal);
      values_low_ecal_vec.push_back(values_low_ecal);
      coeff_low_ecal_vec.push_back(coeffs_low_ecal);
    }

  //  cout<<"size"<<"\t"<<consts_low_ecal_vec.size()<<endl;
  
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


  int ahcal_layer[10]= {43,47,51,55,59,63,67,71,75,79};
  //nentries =4187091;
  nentries =100090;
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
      /////     Matrix initialzation     /////                                                                             
      //////////////////////////////////////////////////////////
      //cout<<"alp_check"<<endl;
      double E_beam =0.0;
      double O = 0.0;
      double sigma2;// = (0.121*0.121)*(total_energy*total_energy) + (1.441*1.441)*total_energy;
      double tot_E_gev= total_energy;
      float chi2_funct=0.0;
      //double ecal_ene_frac_bins[6] = {0.0, 0.06, 0.2, 0.4, 0.6, 1.01};
      for (int ecal_bins=0; ecal_bins<5; ecal_bins++)
      {//cout<<"ecal_bins="<<ecal_bins<<endl; 
        if(ecal_ene_frac>=ecal_ene_frac_bins[ecal_bins] && ecal_ene_frac<ecal_ene_frac_bins[ecal_bins+1])
        {//cout<<"ecal_ene_frac="<<ecal_ene_frac<<endl;
            for(int i_bin=0;i_bin<85;i_bin++)
            {
                bool IsIn_en_range = true;
                if(trueBeamEnergy>=Elist[i_bin] && trueBeamEnergy <=Elist[i_bin]+4)
                {
                    //cout<<"Elist[i_bin]="<<"\t"<<Elist[i_bin]<<"\t"<<"trueBeamEnergy"<<trueBeamEnergy<<endl;

                    if (rechit_shower_start_layer>28)
                    {//cout<<"rechit_shower_start_layer="<<rechit_shower_start_layer<<endl;
                    IsIn_en_range = (total_energy > en_range_FH_xmin_vec[i_bin][Elist[i_bin]+2] && total_energy < en_range_FH_xmax_vec[i_bin][Elist[i_bin]+2]);
                   /* double*/ sigma2 = (0.09*0.09)*(total_energy*total_energy) + (1.23*1.23)*total_energy;
                    //cout<<"sigma2"<<sigma2<<endl; 
                    }
                    else
                    {//cout<<"rechit_shower_start_layer="<<rechit_shower_start_layer<<endl;
                    IsIn_en_range = (tot_E_gev > en_range_EE_xmin_vec[i_bin][Elist[i_bin]+2] && tot_E_gev < en_range_EE_xmax_vec[i_bin][Elist[i_bin]+2]);
                   /* double*/ sigma2 = (0.084*0.084)*(total_energy*total_energy) + (1.39*1.39)*total_energy;
                    //cout<<"sigma2"<<sigma2<<endl;
                    }
                    if(IsIn_en_range){
                        E_beam = Elist[i_bin]+2;
                        if(ecal_bins<1)
                        {
                        
                            
                            //cout<<i_bin<<"\t"<<Elist[i_bin]<<endl;

                            ROOT::Math::SMatrix<double,2, 2, ROOT::Math::MatRepStd<double,2> > coeffs_low_ecal;
                            ROOT::Math::SVector<double, 2> consts_low_ecal;
                            ROOT::Math::SVector<double, 2> values_low_ecal;
                            for(int i = 0; i < 2; i++) {
                                for(int j = 0; j < 2; j++) {
                                coeffs_low_ecal(i,j) = 0.0;
                                }	
                                consts_low_ecal(i) = 0.0;
                                values_low_ecal(i) = 0.0;
                            }
                            //cout<<i_bin<<"\t"<<Elist[i_bin]<<endl;

                            
                            if(chi2_method == 0) // method -input to chi2 are in units of MIPs
                                {
                                coeffs_low_ecal(0,0) = (rechitEnergySum_FH*rechitEnergySum_FH)/sigma2;
                                coeffs_low_ecal(0,1) = (rechitEnergySum_AH*rechitEnergySum_FH)/sigma2;
                                coeffs_low_ecal(1,0) = (rechitEnergySum_FH*rechitEnergySum_AH)/sigma2;
                                coeffs_low_ecal(1,1) = (rechitEnergySum_AH*rechitEnergySum_AH)/sigma2;
                                
                                consts_low_ecal(0)   = (E_beam-O)*rechitEnergySum_FH/sigma2;
                                consts_low_ecal(1)   = (E_beam-O)*rechitEnergySum_AH/sigma2;
                                consts_low_ecal_vec[i_bin]+= consts_low_ecal;
                                coeff_low_ecal_vec[i_bin]+=coeffs_low_ecal;

                                }
                            else // method -input to chi2 are in units of GeV -> We use this one
                                {
                                coeffs_low_ecal(0,0) = (FH_detscale*FH_detscale)/sigma2;
                                coeffs_low_ecal(0,1) = (AH_detscale*FH_detscale)/sigma2;
                                coeffs_low_ecal(1,0) = (FH_detscale*AH_detscale)/sigma2;
                                coeffs_low_ecal(1,1) = (AH_detscale*AH_detscale)/sigma2;
                                
                                consts_low_ecal(0)   = (E_beam-O)*FH_detscale/sigma2;
                                consts_low_ecal(1)   = (E_beam-O)*AH_detscale/sigma2;
                                //cout<<i_bin<<"\t"<<Elist[i_bin]<<"\t"<<"beforE"<<endl;

                                consts_low_ecal_vec[i_bin] += consts_low_ecal;
                                coeff_low_ecal_vec[i_bin]+=coeffs_low_ecal;

                                }
                        }
                        
                        else
                        { // for EH hadrons
                            //cout<<jentry<<"\t"<<"beforE"<<endl;
                            
                                //if(IsIn_en_range){
                            //double sigma2 = (0.084*0.084)*tot_E_gev*tot_E_gev + (1.39*1.39)*tot_E_gev;
                            ROOT::Math::SMatrix<double,3, 3, ROOT::Math::MatRepStd<double,3> > coeffs_high_ecal;
                            ROOT::Math::SVector<double, 3> consts_high_ecal;
                            ROOT::Math::SVector<double, 3> values_high_ecal;

                            for(int i = 0; i < 3; i++) {
                                for(int j = 0; j < 3; j++) {
                                coeffs_high_ecal(i,j) = 0.0;
                                }
                                consts_high_ecal(i) = 0.0;
                                values_high_ecal(i) = 0.0;
                            }
                            

                            if(chi2_method==0)
                                {
                                coeffs_high_ecal(0,0) = (rechitEnergySum_EE*rechitEnergySum_EE)/sigma2;
                                coeffs_high_ecal(0,1) = (rechitEnergySum_EE*rechitEnergySum_FH)/sigma2;
                                coeffs_high_ecal(0,2) = (rechitEnergySum_EE*rechitEnergySum_AH)/sigma2;
                                coeffs_high_ecal(1,0) = (rechitEnergySum_FH*rechitEnergySum_EE)/sigma2;
                                coeffs_high_ecal(1,1) = (rechitEnergySum_FH*rechitEnergySum_FH)/sigma2;
                                coeffs_high_ecal(1,2) = (rechitEnergySum_FH*rechitEnergySum_AH)/sigma2;
                                coeffs_high_ecal(2,0) = (rechitEnergySum_AH*rechitEnergySum_EE)/sigma2;
                                coeffs_high_ecal(2,1) = (rechitEnergySum_AH*rechitEnergySum_FH)/sigma2;
                                coeffs_high_ecal(2,2) = (rechitEnergySum_AH*rechitEnergySum_AH)/sigma2;

                                consts_high_ecal(0)   = (E_beam-O)*rechitEnergySum_EE/sigma2;
                                consts_high_ecal(1)   = (E_beam-O)*rechitEnergySum_FH/sigma2;
                                consts_high_ecal(2)   = (E_beam-O)*rechitEnergySum_AH/sigma2;
                                //consts_high_ecal_vec[i_bin]+= consts_high_ecal;
                                //coeff_high_ecal_vec[i_bin]+=coeffs_high_ecal;
                                consts_high_ecal_vec_diff_bins[ecal_bins-1][i_bin] += consts_high_ecal;
                                coeff_high_ecal_vec_diff_bins[ecal_bins-1][i_bin] += coeffs_high_ecal; 

                                }
                            else
                                {
                                coeffs_high_ecal(0,0) = (EE_detscale*EE_detscale)/sigma2;
                                coeffs_high_ecal(0,1) = (EE_detscale*FH_detscale)/sigma2;
                                coeffs_high_ecal(0,2) = (EE_detscale*AH_detscale)/sigma2;
                                coeffs_high_ecal(1,0) = (FH_detscale*EE_detscale)/sigma2;
                                coeffs_high_ecal(1,1) = (FH_detscale*FH_detscale)/sigma2;
                                coeffs_high_ecal(1,2) = (FH_detscale*AH_detscale)/sigma2;
                                coeffs_high_ecal(2,0) = (AH_detscale*EE_detscale)/sigma2;
                                coeffs_high_ecal(2,1) = (AH_detscale*FH_detscale)/sigma2;
                                coeffs_high_ecal(2,2) = (AH_detscale*AH_detscale)/sigma2;


                                consts_high_ecal(0)   = (E_beam - O)*EE_detscale/sigma2;
                                consts_high_ecal(1)   = (E_beam - O)*FH_detscale/sigma2;
                                consts_high_ecal(2)   = (E_beam - O)*AH_detscale/sigma2;
                                //consts_high_ecal_vec[i_bin]+= consts_high_ecal;
                                //coeff_high_ecal_vec[i_bin]+= coeffs_high_ecal;
                                consts_high_ecal_vec_diff_bins[ecal_bins-1][i_bin] += consts_high_ecal;
                                coeff_high_ecal_vec_diff_bins[ecal_bins-1][i_bin] += coeffs_high_ecal;
                                }

                            

                        }
                    }
                }

            
            }//true bin loop   
        
        float w1=0,w2=0,w3=0;
        float pred = 0.0;

        if(ecal_bins==0)
          {
           w1 = 0.0;
           w2 = f_w2_bin1->Eval(trueBeamEnergy);
           w3 = f_w3_bin1->Eval(trueBeamEnergy);

           pred = EE_detscale + w2*FH_detscale + w3*AH_detscale;
           //pred = w2*FH_detscale + w3*AH_detscale;
          }

        else if(ecal_bins==1)
          {
           w1 = f_w1_bin2->Eval(trueBeamEnergy);;
           w2 = f_w2_bin2->Eval(trueBeamEnergy);
           w3 = f_w3_bin2->Eval(trueBeamEnergy);

           pred = w1*EE_detscale + w2*FH_detscale + w3*AH_detscale;
          }

        else if(ecal_bins==2)
          {
           w1 = f_w1_bin3->Eval(trueBeamEnergy); 
           w2 = f_w2_bin3->Eval(trueBeamEnergy);
           w3 = f_w3_bin3->Eval(trueBeamEnergy);

           pred = w1*EE_detscale + w2*FH_detscale + w3*AH_detscale;
          }

        else if(ecal_bins==3)
          {
           w1 = f_w1_bin4->Eval(trueBeamEnergy);
           w2 = f_w2_bin4->Eval(trueBeamEnergy);
           w3 = f_w3_bin4->Eval(trueBeamEnergy);

           pred = w1*EE_detscale + w2*FH_detscale + w3*AH_detscale;
          }

        
        else
          {
           w1 = f_w1_bin5->Eval(trueBeamEnergy);
           w2 = f_w2_bin5->Eval(trueBeamEnergy);
           w3 = f_w3_bin5->Eval(trueBeamEnergy);

           pred = w1*EE_detscale + w2*FH_detscale + w3*AH_detscale;
          }

        file_ecal_ene_frac<< ecal_ene_frac <<"\t"<< trueBeamEnergy <<"\t"<< pred << "\t"<< ecal_bins+1<< "\t" << rechit_shower_start_layer <<"\n";

        }
        }//end of ecal_bins loop
      if(DEBUG) cout<<"DEBUG: End of Entry = "<<jentry<<endl;
      if(DEBUG) cout<<"DEBUG: ****************** "<<endl;
      //      cout<<"\t"<<"beforE"<<endl;

    } // loop over entries
  
  char* name = new char[1000];


  ////////////////////////////////////////////////////                                                                                     
  //////  Matrix Multiplication and Saving    ////////                                                                                     
  ////////////////////////////////////////////////////                                                                                
  for (int ecal_bins=0; ecal_bins<5; ecal_bins++)
  { 
    sprintf(name,"Results/chi2_opt_ecal_ene_fracFact_%f_to_%f_flatEn.txt",ecal_ene_frac_bins[ecal_bins],ecal_ene_frac_bins[ecal_bins+1]);
    std::ofstream file_bin1;
    file_bin1.open(name,ios::out);
    

    for(int i_en =0; i_en<85; i_en++)
        {
        ROOT::Math::SMatrix<double,3, 3, ROOT::Math::MatRepStd<double,3> > coeffs_high_ecal;
        ROOT::Math::SVector<double, 3> consts_high_ecal;
        ROOT::Math::SVector<double, 3> values_high_ecal;
        ROOT::Math::SMatrix<double,2, 2, ROOT::Math::MatRepStd<double,2> > coeffs_low_ecal;
        ROOT::Math::SVector<double, 2> consts_low_ecal;
        ROOT::Math::SVector<double, 2> values_low_ecal;
        bool isInverted = false;

        if(ecal_bins<1)
        {    
            for(int i = 0; i < 2; i++) {
                for(int j = 0; j < 2; j++) {
                coeffs_low_ecal(i,j) = 0.0;
                }
                consts_low_ecal(i) = 0.0;
                values_low_ecal(i) = 0.0;
            }
            consts_low_ecal = consts_low_ecal_vec[i_en];
            values_low_ecal = values_low_ecal_vec[i_en];
            coeffs_low_ecal = coeff_low_ecal_vec[i_en];

            isInverted = coeffs_low_ecal.Invert();
            for(int i = 0; i <2; i++) {
                for(int j = 0; j < 2; j++) {
                cout<<coeffs_low_ecal(i, j)<<"\t";
                }
                cout<<endl;
            }
            cout<<" "<<consts_low_ecal(0)<<"\t"<<consts_low_ecal(1)<<"\t"<<consts_low_ecal(2)<<endl;
                                                                                                                            
            if(isInverted) {
                values_low_ecal = coeffs_low_ecal*consts_low_ecal;
                //cout<<"H Hadrons => For E = "<<E_beam<<"GeV, w1 = "<<values_low_ecal(0)<<" ;w2 = "<<values_low_ecal(1)<<" ;w3 = "<<values_low_ecal(2)<<";"<<endl;   

                cout<<"weights => For E = "<<Elist[i_en]+2<<"GeV, w1 = 0.0 ;w2 = "<<values_low_ecal(0)<<" ;w3 = "<<values_low_ecal(1)<<";"<<endl;
                cout<<">>>>>>> Uncertainity  : "<<sqrt(coeffs_low_ecal(0,0)) << " " <<sqrt(coeffs_low_ecal(1,1))<<endl;

                file_bin1<<Elist[i_en]+2<<"\t"<<0<<"\t"<<values_low_ecal(0)<<"\t"<<values_low_ecal(1)<<"\t"<<0<<"\t"<<sqrt(coeffs_low_ecal(0,0))<<"\t"<<sqrt(coeffs_low_ecal(1,1))<<"\n";   
            
            }
            else {
                cout<<"Error: Could not invert for low_ecal_ene_frac hadrons..."<<endl;
            }
        }

        else
        {
            
            for(int i = 0; i < 3; i++) {
                for(int j = 0; j < 3; j++) {
                coeffs_high_ecal(i,j) = 0.0;
                }
                consts_high_ecal(i) = 0.0;
                values_high_ecal(i) = 0.0;
            }
            consts_high_ecal = consts_high_ecal_vec_diff_bins[ecal_bins-1][i_en];
            values_high_ecal = values_high_ecal_vec_diff_bins[ecal_bins-1][i_en];
            coeffs_high_ecal = coeff_high_ecal_vec_diff_bins[ecal_bins-1][i_en];
            //cout<<"before invert"<<"\t"<<coeffs_high_ecal<<endl;

            isInverted = coeffs_high_ecal.Invert();// cout<<"before invert"<<"\t"<<coeffs_high_ecal<<endl;
            for(int i = 0; i < 3; i++) {
                for(int j = 0; j < 3; j++) {
                cout<<coeffs_high_ecal(i, j)<<"\t";
                }
                cout<<endl;
            }
            cout<<" "<<consts_high_ecal(0)<<"\t"<<consts_high_ecal(1)<<"\t"<<consts_high_ecal(2)<<endl;
            if(isInverted) {
                values_high_ecal = coeffs_high_ecal*consts_high_ecal;
                cout<<"weights => For E = "<<Elist[i_en]+2<<"GeV, w1 = "<<values_high_ecal(0)<<" ;w2 = "<<values_high_ecal(1)<<" ;w3 = "<<values_high_ecal(2)<<";"<<endl;
                cout<<">>>>>>> Uncertainity  : "<<sqrt(coeffs_high_ecal(0,0))<< " " << sqrt(coeffs_high_ecal(1,1)) << " " << sqrt(coeffs_high_ecal(2,2))<<endl;

                file_bin1<<Elist[i_en]+2<<"\t"<<values_high_ecal(0)<<"\t"<<values_high_ecal(1)<<"\t"<<values_high_ecal(2)<<"\t"<<sqrt(coeffs_high_ecal(0,0))<<"\t"<<sqrt(coeffs_high_ecal(1,1))<<"\t"<<sqrt(coeffs_high_ecal(2,2))<<"\n";
                
            }
            else {
         //       cout<<"Error: Could not invert for high_ecal_ene_frac hadrons..."<<endl;
            }
        }
        }


    file_bin1.close();    
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
