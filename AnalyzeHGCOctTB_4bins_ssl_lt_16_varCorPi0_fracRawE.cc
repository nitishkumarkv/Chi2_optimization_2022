#define AnalyzeHGCOctTB_cxx
#include <iostream>
#include <vector>
#include <cstring>
#include "AnalyzeHGCOctTB.h"
#include "Math/SMatrix.h"
#include "Math/SVector.h"
#include <math.h>
#include<TF1.h>

#include <vector>
#include <algorithm>

using namespace std;



// chip 3022,44,3028

///////////////////////////////

float ene_uptonX0(float n, int ssl, vector<double> r_len, vector<double> ene_lay){
    
    int ssl_idx = ssl-1;
    
    vector<float> dummy(50, 0.0);
    for(int i=0; i<50; i++){
       dummy[i] = abs(r_len[i] - r_len[ssl_idx] - n);
    }
    
    auto it = minmax_element(dummy.begin(), dummy.end());
    int max_lay_idx = distance(dummy.begin(), it.first);
    
    float ene = 0.0;
    for(int lay_=ssl_idx; lay_<=max_lay_idx; lay_++){
        ene += ene_lay[lay_];
    }

    return ene;
}


//////////////////////////////////////

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
  sprintf(path_to_save, "Results_new_chi2_4bins_ssl_lt_16_varCorPi0_fracRawE/");

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
  file_ecal_ene_frac.open("Results_new_chi2_4bins_ssl_lt_16_varCorPi0_fracRawE/ecal_hcal_Eneupto13X0_true_pred_SS.txt",ios::out);  

  int Elist[85] =  {10, 14 , 18 , 22 , 26 , 30,  34 , 38 , 42,  46,  50,  54,  58,  62,  66,  70,  74,  78,
                     82,  86,  90,  94,  98, 102, 106, 110, 114 ,118, 122, 126, 130, 134, 138, 142 ,146, 150,
                     154 ,158 ,162 ,166 ,170 ,174 ,178, 182, 186, 190 ,194 ,198, 202, 206, 210, 214 ,218, 222,
                  226, 230, 234 ,238, 242 ,246 ,250, 254, 258 ,262 ,266, 270, 274, 278, 282, 286, 290, 294,
                     298, 302 ,306, 310, 314, 318 ,322, 326, 330, 334 ,338 ,342 ,346};
  
  float ss_ecal_ene_frac_bins[11] = {0.0, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8, 999999999.0};
  float mips_hcal_ene_frac_bins[6] = {0.0, 0.4, 0.8, 1.2, 999999999.0};
  float var_pi0_frac_bins[5] = {0.0, 0.2, 0.4, 0.6, 999999999.0};

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


  std::vector<double>param_vector;

  sprintf(name1,"Results_new_chi2_4bins_ssl_lt_16_varCorPi0_fracRawE/fitted_parameters_4bins_ssl_lt_16_varCorPi0_fracRawE.txt");
  std::fstream params;
  params.open(name1,ios::in);
 
  if(!params.is_open())
    {
       std::cout << " file not opened" << std::endl;
    }
  else
    {
      double p;
      while (params >> p)
      {
         param_vector.push_back(p);

      }
 
    }

  //fit function MipsInEE

//  TF1* f_w1_bin1_mips = new TF1("f_w1_bin1_mips","sqrt([0]*[0]+[1]/x)+[2]/x+[3]",5,360);
//  f_w1_bin1_mips->FixParameter(0,param_vector[0]);
//  f_w1_bin1_mips->FixParameter(1,param_vector[1]);
//  f_w1_bin1_mips->FixParameter(2,param_vector[2]);
//  f_w1_bin1_mips->FixParameter(3,param_vector[3]);
//
//  TF1* f_w2_bin1_mips = new TF1("f_w2_bin1_mips","sqrt([0]*[0]+[1]/x)+[2]*x+[3]",5,360);
//  f_w2_bin1_mips->FixParameter(0, param_vector[4]);
//  f_w2_bin1_mips->FixParameter(1, param_vector[5]);
//  f_w2_bin1_mips->FixParameter(2, param_vector[6]);
//  f_w2_bin1_mips->FixParameter(3, param_vector[7]);
//
//  TF1* f_w3_bin1_mips = new TF1("f_w3_bin1_mips","sqrt([0]*[0]+[1]/x)+[2]/x+[3]",5,360);
//  f_w3_bin1_mips->FixParameter(0, param_vector[8]);
//  f_w3_bin1_mips->FixParameter(1, param_vector[9]);
//  f_w3_bin1_mips->FixParameter(2, param_vector[10]);
//  f_w3_bin1_mips->FixParameter(3, param_vector[11]);
//
//  TF1* f_w2_bin2_mips = new TF1("f_w2_bin2_mips","sqrt([0]*[0]+[1]/x)+[2]*x+[3]",5,360);
//  f_w2_bin2_mips->FixParameter(0, param_vector[16]);
//  f_w2_bin2_mips->FixParameter(1, param_vector[17]);
//  f_w2_bin2_mips->FixParameter(2, param_vector[18]);
//  f_w2_bin2_mips->FixParameter(3, param_vector[19]);
//
//  TF1* f_w3_bin2_mips = new TF1("f_w3_bin2_mips","sqrt([0]*[0]+[1]*x)+[2]/x+[3]",5,360);
//  f_w3_bin2_mips->FixParameter(0, param_vector[20]);
//  f_w3_bin2_mips->FixParameter(1, param_vector[21]);
//  f_w3_bin2_mips->FixParameter(2, param_vector[22]);
//  f_w3_bin2_mips->FixParameter(3, param_vector[23]);
//
//  TF1* f_w2_bin3_mips = new TF1("f_w2_bin3_mips","sqrt([0]*[0]+[1]/x)+[2]/x+[3]*x",5,360);
//  f_w2_bin3_mips->FixParameter(0, param_vector[28]);
//  f_w2_bin3_mips->FixParameter(1, param_vector[29]);
//  f_w2_bin3_mips->FixParameter(2, param_vector[30]);
//  f_w2_bin3_mips->FixParameter(3, param_vector[31]);
//
//  TF1* f_w3_bin3_mips = new TF1("f_w3_bin3_mips","sqrt([0]*[0]+[1]*x)+[2]*x+[3]",5,360);
//  f_w3_bin3_mips->FixParameter(0, param_vector[32]);
//  f_w3_bin3_mips->FixParameter(1, param_vector[33]);
//  f_w3_bin3_mips->FixParameter(2, param_vector[34]);
//  f_w3_bin3_mips->FixParameter(3, param_vector[35]);
//
//  TF1* f_w2_bin4_mips = new TF1("f_w2_bin4_mips","sqrt([0]*[0]+[1]/x)+[2]/x+[3]",5,360);
//  f_w2_bin4_mips->FixParameter(0, param_vector[40]);
//  f_w2_bin4_mips->FixParameter(1, param_vector[41]);
//  f_w2_bin4_mips->FixParameter(2, param_vector[42]);
//  f_w2_bin4_mips->FixParameter(3, param_vector[43]);
//
//  TF1* f_w3_bin4_mips = new TF1("f_w3_bin4_mips","sqrt([0]*[0]+[1]*x)+[2]*x+[3]",5,360);
//  f_w3_bin4_mips->FixParameter(0, param_vector[44]);
//  f_w3_bin4_mips->FixParameter(1, param_vector[45]);
//  f_w3_bin4_mips->FixParameter(2, param_vector[46]);
//  f_w3_bin4_mips->FixParameter(3, param_vector[47]);

//  TF1* f_w2_bin5_mips = new TF1("f_w2_bin5_mips","sqrt([0]*[0]+[1]/x)+[2]/x+[3]",5,360);
//  f_w2_bin5_mips->FixParameter(0, param_vector[52]);
//  f_w2_bin5_mips->FixParameter(1, param_vector[53]);
//  f_w2_bin5_mips->FixParameter(2, param_vector[54]);
//  f_w2_bin5_mips->FixParameter(3, param_vector[55]);
//
//  TF1* f_w3_bin5_mips = new TF1("f_w3_bin5_mips","sqrt([0]*[0]+[1]/x)+[2]/x+[3]",5,360);
//  f_w3_bin5_mips->FixParameter(0, param_vector[56]);
//  f_w3_bin5_mips->FixParameter(1, param_vector[57]);
//  f_w3_bin5_mips->FixParameter(2, param_vector[58]);
//  f_w3_bin5_mips->FixParameter(3, param_vector[59]);

//fit function SSInEE

  TF1* f_w1_bin1_ss = new TF1("f_w1_bin1_ss","sqrt([0]*[0]+[1]/x)+[2]/x+[3]",5,360);
  f_w1_bin1_ss->FixParameter(0,param_vector[60]);
  f_w1_bin1_ss->FixParameter(1,param_vector[61]);
  f_w1_bin1_ss->FixParameter(2,param_vector[62]);
  f_w1_bin1_ss->FixParameter(3,param_vector[63]);

  TF1* f_w2_bin1_ss = new TF1("f_w2_bin1_ss","sqrt([0]*[0]+[1]/x)+[2]/x+[3]",5,360);
  f_w2_bin1_ss->FixParameter(0, param_vector[64]);
  f_w2_bin1_ss->FixParameter(1, param_vector[65]);
  f_w2_bin1_ss->FixParameter(2, param_vector[66]);
  f_w2_bin1_ss->FixParameter(3, param_vector[67]);

  TF1* f_w3_bin1_ss = new TF1("f_w3_bin1_ss","sqrt([0]*[0]+[1]*x)+[2]*x+[3]",5,360);
  f_w3_bin1_ss->FixParameter(0, param_vector[68]);
  f_w3_bin1_ss->FixParameter(1, param_vector[69]);
  f_w3_bin1_ss->FixParameter(2, param_vector[70]);
  f_w3_bin1_ss->FixParameter(3, param_vector[71]);


  TF1* f_w1_bin2_ss = new TF1("f_w1_bin2_ss","sqrt([0]*[0]+[1]/x)+[2]*x+[3]",5,360);
  f_w1_bin2_ss->FixParameter(0, param_vector[72]);
  f_w1_bin2_ss->FixParameter(1, param_vector[73]);
  f_w1_bin2_ss->FixParameter(2, param_vector[74]);
  f_w1_bin2_ss->FixParameter(3, param_vector[75]);

  TF1* f_w2_bin2_ss = new TF1("f_w2_bin2_ss","sqrt([0]*[0]+[1]/x)+[2]*x+[3]",5,360);
  f_w2_bin2_ss->FixParameter(0, param_vector[76]);
  f_w2_bin2_ss->FixParameter(1, param_vector[77]);
  f_w2_bin2_ss->FixParameter(2, param_vector[78]);
  f_w2_bin2_ss->FixParameter(3, param_vector[79]);

  TF1* f_w3_bin2_ss = new TF1("f_w3_bin2_ss","sqrt([0]*[0]+[1]/x)+[2]/x+[3]",5,360);
  f_w3_bin2_ss->FixParameter(0, param_vector[80]);
  f_w3_bin2_ss->FixParameter(1, param_vector[81]);
  f_w3_bin2_ss->FixParameter(2, param_vector[82]);
  f_w3_bin2_ss->FixParameter(3, param_vector[83]);


  TF1* f_w1_bin3_ss = new TF1("f_w1_bin3_ss","sqrt([0]*[0]+[1]/x)+[2]*x+[3]",5,360);
  f_w1_bin3_ss->FixParameter(0, param_vector[84]);
  f_w1_bin3_ss->FixParameter(1, param_vector[85]);
  f_w1_bin3_ss->FixParameter(2, param_vector[86]);
  f_w1_bin3_ss->FixParameter(3, param_vector[87]);

  TF1* f_w2_bin3_ss = new TF1("f_w2_bin3_ss","sqrt([0]*[0]+[1]/x)+[2]*x+[3]",5,360);
  f_w2_bin3_ss->FixParameter(0, param_vector[88]);
  f_w2_bin3_ss->FixParameter(1, param_vector[89]);
  f_w2_bin3_ss->FixParameter(2, param_vector[90]);
  f_w2_bin3_ss->FixParameter(3, param_vector[91]);

  TF1* f_w3_bin3_ss = new TF1("f_w3_bin3_ss","sqrt([0]*[0]+[1]/x)+[2]*x+[3]",5,360);
  f_w3_bin3_ss->FixParameter(0, param_vector[92]);
  f_w3_bin3_ss->FixParameter(1, param_vector[93]);
  f_w3_bin3_ss->FixParameter(2, param_vector[94]);
  f_w3_bin3_ss->FixParameter(3, param_vector[95]);


  TF1* f_w1_bin4_ss = new TF1("f_w1_bin4_ss","sqrt([0]*[0]+[1]/x)+[2]/x+[3]",5,360);
  f_w1_bin4_ss->FixParameter(0, param_vector[96]);
  f_w1_bin4_ss->FixParameter(1, param_vector[97]);
  f_w1_bin4_ss->FixParameter(2, param_vector[98]);
  f_w1_bin4_ss->FixParameter(3, param_vector[99]);

  TF1* f_w2_bin4_ss = new TF1("f_w2_bin4_ss","sqrt([0]*[0]+[1]/x)+[2]*x+[3]",5,360);
  f_w2_bin4_ss->FixParameter(0, param_vector[100]);
  f_w2_bin4_ss->FixParameter(1, param_vector[101]);
  f_w2_bin4_ss->FixParameter(2, param_vector[102]);
  f_w2_bin4_ss->FixParameter(3, param_vector[103]);

  TF1* f_w3_bin4_ss = new TF1("f_w3_bin4_ss","sqrt([0]*[0]+[1]/x)+[2]/x+[3]",5,360);
  f_w3_bin4_ss->FixParameter(0, param_vector[104]);
  f_w3_bin4_ss->FixParameter(1, param_vector[105]);
  f_w3_bin4_ss->FixParameter(2, param_vector[106]);
  f_w3_bin4_ss->FixParameter(3, param_vector[107]);

//  TF1* f_w1_bin5_ss = new TF1("f_w1_bin5_ss","sqrt([0]*[0]+[1]/x)+[2]/x+[3]*x",5,360);
//  f_w1_bin5_ss->FixParameter(0, param_vector[108]);
//  f_w1_bin5_ss->FixParameter(1, param_vector[109]);
//  f_w1_bin5_ss->FixParameter(2, param_vector[110]);
//  f_w1_bin5_ss->FixParameter(3, param_vector[111]);
//
//  TF1* f_w2_bin5_ss = new TF1("f_w2_bin5_ss","sqrt([0]*[0]+[1]*x)+[2]/x+[3]",5,360);
//  f_w2_bin5_ss->FixParameter(0, param_vector[112]);
//  f_w2_bin5_ss->FixParameter(1, param_vector[113]);
//  f_w2_bin5_ss->FixParameter(2, param_vector[114]);
//  f_w2_bin5_ss->FixParameter(3, param_vector[115]);
//
//  TF1* f_w3_bin5_ss = new TF1("f_w3_bin5_ss","sqrt([0]*[0]+[1]*x)+[2]/x+[3]",5,360);
//  f_w3_bin5_ss->FixParameter(0, param_vector[116]);
//  f_w3_bin5_ss->FixParameter(1, param_vector[117]);
//  f_w3_bin5_ss->FixParameter(2, param_vector[118]);
//  f_w3_bin5_ss->FixParameter(3, param_vector[119]);
//
//
//
// 
//  TF1* f_w1_bin6_ss = new TF1("f_w1_bin6_ss","sqrt([0]*[0]+[1]/x)+[2]*x+[3]",5,360);
//  f_w1_bin6_ss->FixParameter(0,param_vector[120]);
//  f_w1_bin6_ss->FixParameter(1,param_vector[121]);
//  f_w1_bin6_ss->FixParameter(2,param_vector[122]);
//  f_w1_bin6_ss->FixParameter(3,param_vector[123]);
//
//  TF1* f_w2_bin6_ss = new TF1("f_w2_bin6_ss","sqrt([0]*[0]+[1]*x)+[2]/x+[3]",5,360);
//  f_w2_bin6_ss->FixParameter(0, param_vector[124]);
//  f_w2_bin6_ss->FixParameter(1, param_vector[125]);
//  f_w2_bin6_ss->FixParameter(2, param_vector[126]);
//  f_w2_bin6_ss->FixParameter(3, param_vector[127]);
//
//  TF1* f_w3_bin6_ss = new TF1("f_w3_bin6_ss","sqrt([0]*[0]+[1]*x)+[2]/x+[3]",5,360);
//  f_w3_bin6_ss->FixParameter(0, param_vector[128]);
//  f_w3_bin6_ss->FixParameter(1, param_vector[129]);
//  f_w3_bin6_ss->FixParameter(2, param_vector[130]);
//  f_w3_bin6_ss->FixParameter(3, param_vector[131]);
//
//  TF1* f_w1_bin7_ss = new TF1("f_w1_bin7_ss","sqrt([0]*[0]+[1]/x)+[2]*x+[3]",5,360);
//  f_w1_bin7_ss->FixParameter(0, param_vector[132]);
//  f_w1_bin7_ss->FixParameter(1, param_vector[133]);
//  f_w1_bin7_ss->FixParameter(2, param_vector[134]);
//  f_w1_bin7_ss->FixParameter(3, param_vector[135]);
//
//  TF1* f_w2_bin7_ss = new TF1("f_w2_bin7_ss","sqrt([0]*[0]+[1]*x)+[2]/x+[3]",5,360);
//  f_w2_bin7_ss->FixParameter(0, param_vector[136]);
//  f_w2_bin7_ss->FixParameter(1, param_vector[137]);
//  f_w2_bin7_ss->FixParameter(2, param_vector[138]);
//  f_w2_bin7_ss->FixParameter(3, param_vector[139]);
//
//  TF1* f_w3_bin7_ss = new TF1("f_w3_bin7_ss","sqrt([0]*[0]+[1]*x)+[2]/x+[3]",5,360);
//  f_w3_bin7_ss->FixParameter(0, param_vector[140]);
//  f_w3_bin7_ss->FixParameter(1, param_vector[141]);
//  f_w3_bin7_ss->FixParameter(2, param_vector[142]);
//  f_w3_bin7_ss->FixParameter(3, param_vector[143]);
//
//  TF1* f_w1_bin8_ss = new TF1("f_w1_bin8_ss","sqrt([0]*[0]+[1]/x)+[2]*x+[3]",5,360);
//  f_w1_bin8_ss->FixParameter(0, param_vector[144]);
//  f_w1_bin8_ss->FixParameter(1, param_vector[145]);
//  f_w1_bin8_ss->FixParameter(2, param_vector[146]);
//  f_w1_bin8_ss->FixParameter(3, param_vector[147]);
//
//  TF1* f_w2_bin8_ss = new TF1("f_w2_bin8_ss","sqrt([0]*[0]+[1]*x)+[2]/x+[3]",5,360);
//  f_w2_bin8_ss->FixParameter(0, param_vector[148]);
//  f_w2_bin8_ss->FixParameter(1, param_vector[149]);
//  f_w2_bin8_ss->FixParameter(2, param_vector[150]);
//  f_w2_bin8_ss->FixParameter(3, param_vector[151]);  
// 
//  TF1* f_w3_bin8_ss = new TF1("f_w3_bin8_ss","sqrt([0]*[0]+[1]*x)+[2]/x+[3]",5,360);
//  f_w3_bin8_ss->FixParameter(0, param_vector[152]);
//  f_w3_bin8_ss->FixParameter(1, param_vector[153]);
//  f_w3_bin8_ss->FixParameter(2, param_vector[154]);
//  f_w3_bin8_ss->FixParameter(3, param_vector[155]);
//
//  TF1* f_w1_bin9_ss = new TF1("f_w1_bin9_ss","sqrt([0]*[0]+[1]/x)+[2]/x+[3]*x",5,360);
//  f_w1_bin9_ss->FixParameter(0, param_vector[156]);
//  f_w1_bin9_ss->FixParameter(1, param_vector[157]);
//  f_w1_bin9_ss->FixParameter(2, param_vector[158]);
//  f_w1_bin9_ss->FixParameter(3, param_vector[159]);
//
//  TF1* f_w2_bin9_ss = new TF1("f_w2_bin9_ss","sqrt([0]*[0]+[1]/x)+[2]/x+[3]",5,360);
//  f_w2_bin9_ss->FixParameter(0, param_vector[160]);
//  f_w2_bin9_ss->FixParameter(1, param_vector[161]);
//  f_w2_bin9_ss->FixParameter(2, param_vector[162]);
//  f_w2_bin9_ss->FixParameter(3, param_vector[163]);
//
//  TF1* f_w3_bin9_ss = new TF1("f_w3_bin9_ss","sqrt([0]*[0]+[1]*x)+[2]/x+[3]",5,360);
//  f_w3_bin9_ss->FixParameter(0, param_vector[164]);
//  f_w3_bin9_ss->FixParameter(1, param_vector[165]);
//  f_w3_bin9_ss->FixParameter(2, param_vector[166]);
//  f_w3_bin9_ss->FixParameter(3, param_vector[167]);

//  TF1* f_w1_bin10_ss = new TF1("f_w1_bin10_ss","sqrt([0]*[0]+[1]/x)+[2]/x+[3]",5,360);
//  f_w1_bin10_ss->FixParameter(0, param_vector[168]);
//  f_w1_bin10_ss->FixParameter(1, param_vector[169]);
//  f_w1_bin10_ss->FixParameter(2, param_vector[170]);
//  f_w1_bin10_ss->FixParameter(3, param_vector[171]);
//
//  TF1* f_w2_bin10_ss = new TF1("f_w2_bin10_ss","sqrt([0]*[0]+[1]/x)+[2]/x+[3]",5,360);
//  f_w2_bin10_ss->FixParameter(0, param_vector[172]);
//  f_w2_bin10_ss->FixParameter(1, param_vector[173]);
//  f_w2_bin10_ss->FixParameter(2, param_vector[174]);
//  f_w2_bin10_ss->FixParameter(3, param_vector[175]);
//
//  TF1* f_w3_bin10_ss = new TF1("f_w3_bin10_ss","sqrt([0]*[0]+[1]/x)+[2]/x+[3]",5,360);
//  f_w3_bin10_ss->FixParameter(0, param_vector[176]);
//  f_w3_bin10_ss->FixParameter(1, param_vector[177]);
//  f_w3_bin10_ss->FixParameter(2, param_vector[178]);
//  f_w3_bin10_ss->FixParameter(3, param_vector[179]); 
 

  
   std::vector<ROOT::Math::SVector<double, 3> >consts_ssInEE_vec;
   std::vector<ROOT::Math::SVector<double, 3> >values_ssInEE_vec;
   std::vector<ROOT::Math::SMatrix<double,3, 3, ROOT::Math::MatRepStd<double,3>> > coeff_ssInEE_vec;
   ROOT::Math::SMatrix<double,3, 3, ROOT::Math::MatRepStd<double,3> > coeffs_ssInEE;
   ROOT::Math::SVector<double, 3> consts_ssInEE;
   ROOT::Math::SVector<double, 3> values_ssInEE;
 
   std::vector<std::vector<ROOT::Math::SVector<double, 3> > > consts_ssInEE_vec_diff_bins;
   std::vector<std::vector<ROOT::Math::SVector<double, 3> > > values_ssInEE_vec_diff_bins;
   std::vector<std::vector<ROOT::Math::SMatrix<double,3, 3, ROOT::Math::MatRepStd<double,3>> > >coeff_ssInEE_vec_diff_bins;

   vector <double> ene_lay_ecal(40, 0.0);
   
   for (int pi0Frac_bin=0; pi0Frac_bin<4; pi0Frac_bin++){
     for(int i =0; i<85; i++)
         {
         for(int i = 0; i < 3; i++) {
             for(int j = 0; j < 3; j++) {
             coeffs_ssInEE(i,j) = 0.0;
             }
             consts_ssInEE(i) = 0.0;
             values_ssInEE(i) = 0.0;
         }
             consts_ssInEE_vec.push_back(consts_ssInEE);
             values_ssInEE_vec.push_back(values_ssInEE);
             coeff_ssInEE_vec.push_back(coeffs_ssInEE);
         }
 
         consts_ssInEE_vec_diff_bins.push_back(consts_ssInEE_vec);
         values_ssInEE_vec_diff_bins.push_back(values_ssInEE_vec);
         coeff_ssInEE_vec_diff_bins.push_back(coeff_ssInEE_vec);
     }
 
   //  cout<<"size"<<"\t"<<consts_mipsInEE_vec.size()<<endl;
   
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

  vector <double> rad_len = {1.00114,  1.98059,  2.91411,  3.89356,  4.82708,  5.80653,
        6.74005,  7.7195 ,  8.65302,  9.63247, 10.566  , 11.5454 ,
       12.479  , 13.4584 , 14.3919 , 15.3714 , 16.3049 , 17.2844 ,
       18.2179 , 19.1973 , 20.1308 , 21.2774 , 22.211  , 23.3575 ,
       24.2911 , 25.5403 , 26.4738 , 27.7602 , 30.597  , 33.4725 ,
       36.348  , 39.2235 , 42.099  , 44.8892 , 48.7909 , 51.6664 ,
       54.5419 , 57.6129 , 60.6003 , 63.3922, 67.4867,  71.6301,  75.7736,  79.917 ,  84.0604,  88.2038,
        92.3473,  96.4907, 100.634 , 106.735};


  int ahcal_layer[10]= {43,47,51,55,59,63,67,71,75,79};
  //nentries =4187091;
  //nentries =1000000;
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


      for (int lay_=0; lay_<40; lay_++){
          ene_lay_ecal[lay_] = 0.0;
      }

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
    	//double ecal_ene_frac = EE_detscale/total_energy;
      	//double hcal_ene_frac = FH_detscale/total_energy;
      	double ecal_ene_frac = EE_detscale/trueBeamEnergy;
      	double hcal_ene_frac = FH_detscale/trueBeamEnergy;


      if (rechit_shower_start_layer<=15)
      {
     
          for( int i = 0 ; i < rechit_energy->size(); i++)
          {
             int rechit_lay_idx = rechit_layer->at(i)-1;
             float energy = rechit_energy->at(i);
             ene_lay_ecal[rechit_lay_idx] += (energy*0.0105); 
          }
      }


      double ene_upto_13X0_frac = ene_uptonX0(13, rechit_shower_start_layer, rad_len, ene_lay_ecal)/total_energy; 
      //cout<<"ene_upto_13X0_frac = "<<ene_upto_13X0_frac<<"  ssl = "<<rechit_shower_start_layer<<endl;

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

      
      if (rechit_shower_start_layer<=15) 
      {
        for (int ecal_bins=0; ecal_bins<4; ecal_bins++)
        {//cout<<"ecal_bins="<<ecal_bins<<endl; 
          if(ene_upto_13X0_frac>=var_pi0_frac_bins[ecal_bins] && ene_upto_13X0_frac<var_pi0_frac_bins[ecal_bins+1])
          {//cout<<"ecal_ene_frac="<<ecal_ene_frac<<endl;
              for(int i_bin=0;i_bin<85;i_bin++)
              {
                  bool IsIn_en_range = true;
                  if(trueBeamEnergy>=Elist[i_bin] && trueBeamEnergy <=Elist[i_bin]+4)
                  {
                      //cout<<"Elist[i_bin]="<<"\t"<<Elist[i_bin]<<"\t"<<"trueBeamEnergy"<<trueBeamEnergy<<endl;

                    //cout<<"rechit_shower_start_layer="<<rechit_shower_start_layer<<endl;
                      IsIn_en_range = (tot_E_gev > en_range_EE_xmin_vec[i_bin][Elist[i_bin]+2] && tot_E_gev < en_range_EE_xmax_vec[i_bin][Elist[i_bin]+2]);
                      sigma2 = (0.084*0.084)*(total_energy*total_energy) + (1.39*1.39)*total_energy;
                      //cout<<"sigma2"<<sigma2<<endl;
                      
                      if(IsIn_en_range){
                          E_beam = Elist[i_bin]+2;


                        // for EH hadrons
                              //cout<<jentry<<"\t"<<"beforE"<<endl;
                              
                                  //if(IsIn_en_range){
                              //double sigma2 = (0.084*0.084)*tot_E_gev*tot_E_gev + (1.39*1.39)*tot_E_gev;
                              ROOT::Math::SMatrix<double,3, 3, ROOT::Math::MatRepStd<double,3> > coeffs_ssInEE;
                              ROOT::Math::SVector<double, 3> consts_ssInEE;
                              ROOT::Math::SVector<double, 3> values_ssInEE;

                              for(int i = 0; i < 3; i++) {
                                  for(int j = 0; j < 3; j++) {
                                  coeffs_ssInEE(i,j) = 0.0;
                                  }
                                  consts_ssInEE(i) = 0.0;
                                  values_ssInEE(i) = 0.0;
                              }
                              

                              if(chi2_method==0)
                                  {
                                  coeffs_ssInEE(0,0) = (rechitEnergySum_EE*rechitEnergySum_EE)/sigma2;
                                  coeffs_ssInEE(0,1) = (rechitEnergySum_EE*rechitEnergySum_FH)/sigma2;
                                  coeffs_ssInEE(0,2) = (rechitEnergySum_EE*rechitEnergySum_AH)/sigma2;
                                  coeffs_ssInEE(1,0) = (rechitEnergySum_FH*rechitEnergySum_EE)/sigma2;
                                  coeffs_ssInEE(1,1) = (rechitEnergySum_FH*rechitEnergySum_FH)/sigma2;
                                  coeffs_ssInEE(1,2) = (rechitEnergySum_FH*rechitEnergySum_AH)/sigma2;
                                  coeffs_ssInEE(2,0) = (rechitEnergySum_AH*rechitEnergySum_EE)/sigma2;
                                  coeffs_ssInEE(2,1) = (rechitEnergySum_AH*rechitEnergySum_FH)/sigma2;
                                  coeffs_ssInEE(2,2) = (rechitEnergySum_AH*rechitEnergySum_AH)/sigma2;

                                  consts_ssInEE(0)   = (E_beam-O)*rechitEnergySum_EE/sigma2;
                                  consts_ssInEE(1)   = (E_beam-O)*rechitEnergySum_FH/sigma2;
                                  consts_ssInEE(2)   = (E_beam-O)*rechitEnergySum_AH/sigma2;
                                  //consts_ssInEE_vec[i_bin]+= consts_ssInEE;
                                  //coeff_ssInEE_vec[i_bin]+=coeffs_ssInEE;
                                  consts_ssInEE_vec_diff_bins[ecal_bins][i_bin] += consts_ssInEE;
                                  coeff_ssInEE_vec_diff_bins[ecal_bins][i_bin] += coeffs_ssInEE; 

                                  }
                              else
                                  {
                                  coeffs_ssInEE(0,0) = (EE_detscale*EE_detscale)/sigma2;
                                  coeffs_ssInEE(0,1) = (EE_detscale*FH_detscale)/sigma2;
                                  coeffs_ssInEE(0,2) = (EE_detscale*AH_detscale)/sigma2;
                                  coeffs_ssInEE(1,0) = (FH_detscale*EE_detscale)/sigma2;
                                  coeffs_ssInEE(1,1) = (FH_detscale*FH_detscale)/sigma2;
                                  coeffs_ssInEE(1,2) = (FH_detscale*AH_detscale)/sigma2;
                                  coeffs_ssInEE(2,0) = (AH_detscale*EE_detscale)/sigma2;
                                  coeffs_ssInEE(2,1) = (AH_detscale*FH_detscale)/sigma2;
                                  coeffs_ssInEE(2,2) = (AH_detscale*AH_detscale)/sigma2;


                                  consts_ssInEE(0)   = (E_beam - O)*EE_detscale/sigma2;
                                  consts_ssInEE(1)   = (E_beam - O)*FH_detscale/sigma2;
                                  consts_ssInEE(2)   = (E_beam - O)*AH_detscale/sigma2;
                                  //consts_ssInEE_vec[i_bin]+= consts_ssInEE;
                                  //coeff_ssInEE_vec[i_bin]+= coeffs_ssInEE;
                                  consts_ssInEE_vec_diff_bins[ecal_bins][i_bin] += consts_ssInEE;
                                  coeff_ssInEE_vec_diff_bins[ecal_bins][i_bin] += coeffs_ssInEE;
                                  }
                          
                      }
                  }

              
              }//true bin loop   
          
          float w1=0,w2=0,w3=0;
          float pred = 0.0;

          if(ecal_bins==0)
            {
            w1 = f_w1_bin1_ss->Eval(trueBeamEnergy);
            w2 = f_w2_bin1_ss->Eval(trueBeamEnergy);
            w3 = f_w3_bin1_ss->Eval(trueBeamEnergy);

            if(w1<0.0) w1=0.0;
            if(w2<0.0) w2=0.0;
            if(w3<0.0) w3=0.0;

            pred = w1*EE_detscale + w2*FH_detscale + w3*AH_detscale;
            }

          else if(ecal_bins==1)
            {
            w1 = f_w1_bin2_ss->Eval(trueBeamEnergy);;
            w2 = f_w2_bin2_ss->Eval(trueBeamEnergy);
            w3 = f_w3_bin2_ss->Eval(trueBeamEnergy);

            if(w1<0.0) w1=0.0;
            if(w2<0.0) w2=0.0;
            if(w3<0.0) w3=0.0;

            pred = w1*EE_detscale + w2*FH_detscale + w3*AH_detscale;
            }

          else if(ecal_bins==2)
            {
            w1 = f_w1_bin3_ss->Eval(trueBeamEnergy); 
            w2 = f_w2_bin3_ss->Eval(trueBeamEnergy);
            w3 = f_w3_bin3_ss->Eval(trueBeamEnergy);

            if(w1<0.0) w1=0.0;
            if(w2<0.0) w2=0.0;
            if(w3<0.0) w3=0.0;

            pred = w1*EE_detscale + w2*FH_detscale + w3*AH_detscale;
            }

          else
            {
            w1 = f_w1_bin4_ss->Eval(trueBeamEnergy);
            w2 = f_w2_bin4_ss->Eval(trueBeamEnergy);
            w3 = f_w3_bin4_ss->Eval(trueBeamEnergy);

            if(w1<0.0) w1=0.0;
            if(w2<0.0) w2=0.0;
            if(w3<0.0) w3=0.0;

            pred = w1*EE_detscale + w2*FH_detscale + w3*AH_detscale;
            }

          
       //   else if(ecal_bins==4)
       //     {
       //     w1 = f_w1_bin5_ss->Eval(trueBeamEnergy);
       //     w2 = f_w2_bin5_ss->Eval(trueBeamEnergy);
       //     w3 = f_w3_bin5_ss->Eval(trueBeamEnergy);

       //     if(w1<0.0) w1=0.0;
       //     if(w2<0.0) w2=0.0;
       //     if(w3<0.0) w3=0.0;

       //     pred = w1*EE_detscale + w2*FH_detscale + w3*AH_detscale;
       //     }

       //   else if(ecal_bins==5)
       //     {
       //     w1 = f_w1_bin6_ss->Eval(trueBeamEnergy);
       //     w2 = f_w2_bin6_ss->Eval(trueBeamEnergy);
       //     w3 = f_w3_bin6_ss->Eval(trueBeamEnergy);

       //     if(w1<0.0) w1=0.0;
       //     if(w2<0.0) w2=0.0;
       //     if(w3<0.0) w3=0.0;

       //     pred = w1*EE_detscale + w2*FH_detscale + w3*AH_detscale;
       //     }

       //   else if(ecal_bins==6)
       //     {
       //     w1 = f_w1_bin7_ss->Eval(trueBeamEnergy);;
       //     w2 = f_w2_bin7_ss->Eval(trueBeamEnergy);
       //     w3 = f_w3_bin7_ss->Eval(trueBeamEnergy);

       //     if(w1<0.0) w1=0.0;
       //     if(w2<0.0) w2=0.0;
       //     if(w3<0.0) w3=0.0;

       //     pred = w1*EE_detscale + w2*FH_detscale + w3*AH_detscale;
       //     }

       //   else if(ecal_bins==7)
       //     {
       //     w1 = f_w1_bin8_ss->Eval(trueBeamEnergy); 
       //     w2 = f_w2_bin8_ss->Eval(trueBeamEnergy);
       //     w3 = f_w3_bin8_ss->Eval(trueBeamEnergy);

       //     if(w1<0.0) w1=0.0;
       //     if(w2<0.0) w2=0.0;
       //     if(w3<0.0) w3=0.0;

       //     pred = w1*EE_detscale + w2*FH_detscale + w3*AH_detscale;
       //     }

       //   else 
       //     {
       //     w1 = f_w1_bin9_ss->Eval(trueBeamEnergy);
       //     w2 = f_w2_bin9_ss->Eval(trueBeamEnergy);
       //     w3 = f_w3_bin9_ss->Eval(trueBeamEnergy);

       //     if(w1<0.0) w1=0.0;
       //     if(w2<0.0) w2=0.0;
       //     if(w3<0.0) w3=0.0;

       //     pred = w1*EE_detscale + w2*FH_detscale + w3*AH_detscale;
       //     }

          //else
          //  {
          //  w1 = f_w1_bin10_ss->Eval(trueBeamEnergy);
          //  w2 = f_w2_bin10_ss->Eval(trueBeamEnergy);
          //  w3 = f_w3_bin10_ss->Eval(trueBeamEnergy);

          //  pred = w1*EE_detscale + w2*FH_detscale + w3*AH_detscale;
          //  }


          file_ecal_ene_frac<< ecal_ene_frac <<"\t"<< hcal_ene_frac <<"\t"<< ene_upto_13X0_frac <<"\t"<< trueBeamEnergy <<"\t"<< pred << "\t"<< ecal_bins+1<< "\t" << rechit_shower_start_layer<<"\n";

          }
          }//end of ecal_bins loop

      }

      if(DEBUG) cout<<"DEBUG: End of Entry = "<<jentry<<endl;
      if(DEBUG) cout<<"DEBUG: ****************** "<<endl;
      //      cout<<"\t"<<"beforE"<<endl;

    } // loop over entries
  
  char* name = new char[1000];


  ////////////////////////////////////////////////////                                                                                     
  //////  Matrix Multiplication and Saving    ////////                                                                                     
  ////////////////////////////////////////////////////    
                                                                            

   	   for (int ecal_bins=0; ecal_bins<4; ecal_bins++)
   	   { 
   	   sprintf(name,"Results_new_chi2_4bins_ssl_lt_16_varCorPi0_fracRawE/chi2_opt_ssInEE_ecal_ene_fracFact_%f_to_%f_flatEn.txt",var_pi0_frac_bins[ecal_bins],var_pi0_frac_bins[ecal_bins+1]);
   	   std::ofstream file_bin1;
   	   file_bin1.open(name,ios::out);

   	   for(int i_en =0; i_en<85; i_en++)
   	       {
   	       ROOT::Math::SMatrix<double,3, 3, ROOT::Math::MatRepStd<double,3> > coeffs_ssInEE;
   	       ROOT::Math::SVector<double, 3> consts_ssInEE;
   	       ROOT::Math::SVector<double, 3> values_ssInEE;
   	       bool isInverted = false;


   	           for(int i = 0; i < 3; i++) {
   	               for(int j = 0; j < 3; j++) {
   	               coeffs_ssInEE(i,j) = 0.0;
   	               }
   	               consts_ssInEE(i) = 0.0;
   	               values_ssInEE(i) = 0.0;
   	           }
   	           consts_ssInEE = consts_ssInEE_vec_diff_bins[ecal_bins][i_en];
   	           values_ssInEE = values_ssInEE_vec_diff_bins[ecal_bins][i_en];
   	           coeffs_ssInEE = coeff_ssInEE_vec_diff_bins[ecal_bins][i_en];
   	           //cout<<"before invert"<<"\t"<<coeffs_ssInEE<<endl;

   	           isInverted = coeffs_ssInEE.Invert();
   	           // cout<<"before invert"<<"\t"<<coeffs_ssInEE<<endl;
   	           for(int i = 0; i < 3; i++) {
   	               for(int j = 0; j < 3; j++) {
   	              // cout<<coeffs_ssInEE(i, j)<<"\t";
   	               }
   	              // cout<<endl;
   	           }
   	           //cout<<" "<<consts_ssInEE(0)<<"\t"<<consts_ssInEE(1)<<"\t"<<consts_ssInEE(2)<<endl;
   	           if(isInverted) {
   	               values_ssInEE = coeffs_ssInEE*consts_ssInEE;
   	              // cout<<"EH Hadrons => For E = "<<Elist[i_en]+2<<"GeV, w1 = "<<values_ssInEE(0)<<" ;w2 = "<<values_ssInEE(1)<<" ;w3 = "<<values_ssInEE(2)<<"\t"<<ecal_bins<<"\t"<<1<<endl;
   	              // cout<<">>>>>>> Uncertainity  : "<<sqrt(coeffs_ssInEE(0,0))<< " " << sqrt(coeffs_ssInEE(1,1)) << " " << sqrt(coeffs_ssInEE(2,2))<<endl;

   	               file_bin1<<Elist[i_en]+2<<"\t"<<values_ssInEE(0)<<"\t"<<values_ssInEE(1)<<"\t"<<values_ssInEE(2)<<"\t"<<sqrt(coeffs_ssInEE(0,0))<<"\t"<<sqrt(coeffs_ssInEE(1,1))<<"\t"<<sqrt(coeffs_ssInEE(2,2))<<"\n";
   	               
   	           }
   	           else {
   	             // cout<<"Error: Could not invert for ssInEE_ene_frac hadrons..."<<endl;
   	             // cout<<"EH Hadrons => For E = "<<Elist[i_en]+2<<"GeV, w1 = "<<values_ssInEE(0)<<" ;w2 = "<<values_ssInEE(1)<<" ;w3 = "<<values_ssInEE(2)<<"\t"<<ecal_bins<<"\t"<<1<<endl;
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
