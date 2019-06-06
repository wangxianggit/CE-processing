#include "mapping.h"
#include <TH1.h>
#include <TF1.h>
#include <TH2D.h>
#include <TFitResultPtr.h>
#include <TStyle.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TRandom3.h>
#include <fstream>
#include <map>
#include <iostream>
#include <algorithm>
#include <unistd.h>
#include <fstream>
#include <sstream>
#include <vector>
#include "TVirtualFFT.h"
using namespace std;

 Double_t xpar[128][2];
 Double_t ypar[48][2];
 
 Double_t xpar_xia[128][2];
 Double_t ypar_xia[48][2];
 Double_t sl=2.56038,cept=-33.8167;//linear
 //Double_t sl=2.79903,cept=-75.0199;//pol2

void mapping::Readpar(){
  //ifstream in("fbnorm.txt");
  ifstream in("EF_equal.dat");
  if(!in.is_open()){
	 cout<<"cannot open file !"<<endl;
	 for(Int_t i=0;i<128;i++){
		xpar[i][0]=1;
		xpar[i][0]=0;
	}
	 for(Int_t i=0;i<48;i++){
		ypar[i][0]=1;
		ypar[i][0]=0;
	}
	 return;
  }
  Double_t xp,yp;
  Int_t ch;
  while(1){
    in>>ch>>yp>>xp;
	if(!in.good()) break;
	if(ch<128) {
	  xpar[ch][0]=xp;
	  xpar[ch][1]=yp;
	  //cout<<ch<<"  "<<xp<<"  "<<yp<<"  "<<zp<<endl;
	} else{
	  ypar[ch-128][0]=xp;
	  ypar[ch-128][1]=yp;
	  //cout<<ch<<"  "<<xp<<"  "<<yp<<"  "<<zp<<endl;
	}
  }
}

 TH1D *hx[128],*hnx,*htx,*hehx;
 TH1D *hy[48],*hny,*hty,*hehy;
 Int_t nx[128],ny[48];

void mapping::SetHist()
{
  TString sx,sy;
  Int_t nx[128],ny[48];
  for(Int_t i=0;i<128;i++) {
    sx.Form("hx%03i",i);
    hx[i]=new TH1D(sx.Data(),sx.Data(),740,0,740);
    nx[i]=0;
    if(i<48) {
      sy.Form("hy%03i",i);
      hy[i]=new TH1D(sy.Data(),sy.Data(),740,0,740);
      ny[i]=0;
    } 
  }
    hnx=new TH1D("hnx","hnx",128,0,128);
    hny=new TH1D("hny","hny",48,0,48);
    htx=new TH1D("htx","htx",128,0,128);
    hty=new TH1D("hty","hty",48,0,48);
    hehx=new TH1D("hehx","hehx",128,0,128);
    hehy=new TH1D("hehy","hehy",48,0,48);
}

TF1 *f1;
TF1 *f2;
TF1 *f3;
TGraph *gr = new TGraph;
TCanvas *c1 = new TCanvas;
TGraph *gra = new TGraph;
void mapping::FillHist()
{
  if(det!=2) return;
  if(tsp[0]<125||tsp[0]>137) return;
  if(np!=1) return;
  if(gfl==1) return;
/*  delete c1;
  c1=new TCanvas;
  delete gra;
  gra=new TGraph;
       for(int i=0;i<750;i++)
        gra->SetPoint(i,i,dwave[i]);
        
        gra->Draw();
        c1->SaveAs("c1.eps");
        sleep(5);*/

  if(fr==2) {
    Int_t x=str;
    if(nx[x]<1000 && e>4000 && e<10000) {//choose a peak.5200-5140//4800-5800
    //align for tsp[0]==132, superpulse range from 0-740;
     Int_t tx0 = tsp[0];
     for(Int_t i=0;i<50;i++)
       hx[x]->Fill(i,0);
     
     for(Int_t i=132;i>0;i--){
       hx[x]->Fill(i,dwave[tx0]*2000./er); //since the tx0
       tx0--;
       if(tx0<1) break;
       }
     
     tx0 = tsp[0]+1;
     for(Int_t i=133;i<740;i++){
       hx[x]->Fill(i,dwave[tx0]*2000./er);
       tx0++;
       if(tx0>749) break;
       }
       nx[x]++;
   }
  }
  
  if(fr==1) {
    Int_t y=str;
    if(ny[y]<1000 && e>4000 && e<10000) {
     Int_t ty0 = tsp[0];
     for(Int_t i=0;i<50;i++)
       hy[y]->Fill(i,0);
     
     for(Int_t i=132;i>0;i--){
       hy[y]->Fill(i,dwave[ty0]*2000./er);
       ty0--;
       if(ty0<1)break;
       }
     
     ty0 = tsp[0]+1;//repeat fill increase the bin value
     for(Int_t i=133;i<740;i++) {
       hy[y]->Fill(i,dwave[ty0]*2000./er);
       ty0++;
       if(ty0>749) break;
       }
       ny[y]++;
    }
  }
}

void mapping::WriteHist(Int_t runnum)
{
  char opfh[124];
  number=runnum;
  sprintf(opfh,"./hist%05d.root",number);
  TFile *opf1=new TFile(opfh,"RECREATE");
  for(Int_t i=0;i<128;i++) {
    if(nx[i]<1)continue;
    hx[i]->Scale(1./nx[i]);
    hx[i]->Write();
    hnx->Fill(i,nx[i]);
  //  cout<<"i,nx[i]:"<<i<<" "<<nx[i]<<endl;
    for(Int_t j=0;j<750;j++) dwave[j]=hx[i]->GetBinContent(j+1);
    ResetOpt();
    det=2;
    fr=2;
    str=i;
    PuDet();
    htx->Fill(i,tsp[0]);
    hehx->Fill(i,eh[0]);
  }

  for(Int_t i=0;i<48;i++) {
   if(ny[i]<1)continue;
   hy[i]->Scale(1./ny[i]);
   hy[i]->Write();
   hny->Fill(i,ny[i]);
    for(Int_t j=0;j<750;j++) dwave[j]=hy[i]->GetBinContent(j+1);
    ResetOpt();
    det=2;
    fr=1;
    str=i;
    PuDet();
    hty->Fill(i,tsp[0]);
    hehy->Fill(i,eh[0]);
  //  if(np==1) cout<<fr<<" y "<<str<<" "<<np<<" "<<tsp[0]<<" "<<eh[0]<<endl;
 //   if(np==2) cout<<fr<<" y "<<str<<" "<<np<<" "<<tsp[0]<<" "<<eh[0]<<" "<<tsp[1]<<" "<<eh[1]<<endl;  
  // cout<<"i,ny[i]:"<<i<<" "<<ny[i]<<endl;
   }
  hnx->Write();
  hny->Write();
  htx->Write();
  hty->Write();
  hehx->Write();
  hehy->Write();
  opf1->Close();
}

  TRandom3 *ru;
void mapping::Loop(TTree *ipt,Int_t bdn,TTree *opt)
{
   Init(ipt);
   boardn=bdn;
   ru=new TRandom3();

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   
   Long64_t nbytes = 0, nb = 0;
   Int_t be = 2705;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      //if(Data.TriggerTimeTag!=3981594868210)continue;
      // if (Cut(ientry) < 0) continue;
      if(det==4)continue;
	  ResetOpt();
	  if(Mapping()){
	    opt->Fill();
            //if(jentry*1000/nentries%50==0) cout<<" processing "<<Int_t(jentry*1000./nentries)/10.<<"%"<<endl;
#ifdef HIST
            FillHist();//fill histogram
#endif
	    nevt++;
	   // if(jentry%1000==0)cout<<"Processing..."<<"  "<<jentry<<"  Total events:"<<nentries<<endl;
	   }
   }
}

mapping::mapping(TTree *opt) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  if(opt==0) {
    cout<<"no output tree!"<<endl;
	exit(1);
  }
  BranchOpt(opt);
  nevt=0;
}

mapping::~mapping()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t mapping::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t mapping::LoadTree(Long64_t entry)
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

void mapping::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // poInt_ters of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch poInt_ters
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("BankNo", &BankNo, &b_BankNo);
   fChain->SetBranchAddress("ChannelNo", &ChannelNo, &b_ChannelNo);
   fChain->SetBranchAddress("SysTime", &SysTime, &b_SysTime);
   fChain->SetBranchAddress("Data", &Data, &b_Data);
   fChain->SetBranchAddress("Sample", Sample, &b_Sample);
   Notify();
}

Bool_t mapping::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void mapping::Show(Long64_t entry)
{
// PrInt_t contents of entry.
// If entry is not specified, prInt_t current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t mapping::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void mapping::ResetOpt(){
  np=0;
  trig0=125;
  ene_thr=1;//1MeV
  det=0;
  fr=0;
  str=-1;
  er=0;
  flu=0;
  e=0;
  st=false;
  thr1=0;
  thr2=0;
  min=0;
  memset(ratio,0,sizeof(ratio));//trapezoid calibrate
  memset(trap,0,sizeof(trap));
  memset(trap1,0,sizeof(trap1));
  memset(msd,0,sizeof(msd));
  memset(tsp,0,sizeof(tsp));
  memset(msd1,0,sizeof(msd1));
  memset(tsp1,0,sizeof(tsp1));

}

Bool_t mapping::BranchOpt(TTree *opt){
  if(opt==0) return 0;
  opt->Branch("det",&det,"det/I");
  opt->Branch("fr",&fr,"fr/I");
  opt->Branch("str",&str,"str/I");
  opt->Branch("ts",&ts,"ts/l");
  opt->Branch("er",&er,"er/D");
  opt->Branch("e",&e,"e/D");
  opt->Branch("np",&np,"np/I");
  opt->Branch("tsp",&tsp,"tsp[np]/l");
  opt->Branch("eh",&eh,"eh[np]/D");
  opt->Branch("msd",&msd,"msd[np]/D");
  opt->Branch("flu",&flu,"flu/D");
#ifdef TEST
  opt->Branch("a2",&a2,"a2/D");
  opt->Branch("off2",&off2,"off2/I");
#endif
  char varp[120];
  opt->Branch("dwave",&dwave,"dwave[750]/D");
  sprintf(varp,"rwave[%d]/S",NSAMPLE_V1724);
  opt->Branch("rwave",&rwave,varp);
  opt->Branch("trap",&trap,"trap[750]/D");
  opt->Branch("trap1",&trap1,"trap1[750]/D");
  opt->Branch("pu",&npu,"npu/I");
  opt->Branch("sample",sample,"sample[750]/S");
  
  opt->Branch("bd",&boardn,"bd/I");
  opt->Branch("ch",&ChannelNo,"ch/b");
  opt->Branch("nevt",&nevt,"nevt/I");
  opt->Branch("st",&st,"st/B");
  opt->Branch("gfl",&gfl,"gfl/I");
  return 1;
}

Bool_t mapping::Mapping(){
 
  if(boardn>=0&&boardn<16){
    det=2;//dssd
    fr=2;
    str=boardn*8+ChannelNo;
  }
  else if(boardn>=16&&boardn<22) {
    det=2;
    fr=1;
    str=(boardn-16)*8+ChannelNo;
  }
  else if(boardn==22&&ChannelNo==5){
    det=1;//mwpc
  }
  else if(boardn==24&&ChannelNo>=0&&ChannelNo<5){
    //return 0; 
    det=4;//gamma detector
    str=ChannelNo;
  }
  else if(boardn==22&&ChannelNo<3&&ChannelNo>=0){
    det=3;//veto
    str=ChannelNo;
  }
  else if(boardn==23&&ChannelNo==0){
    det=2;
    fr=2;
    str=45;
  }
  else {
    return 0;
  }
  
  //if(ChannelNo!=0)return 0;
  ts=Data.TriggerTimeTag;
  er=Data.Energy;
  
  //if(er<=0) return 0;
  if(det==2){
  if(fr==2) e=xpar[str][0]+xpar[str][1]*er;
  if(fr==1) e=ypar[str][0]+ypar[str][1]*er;
  
  e=e*sl+cept;
  }
  else e=er;
  //if(det!=2&&er>0)cout<<e<<"  "<<det<<endl;
  memcpy(rwave,Data.Probe2,sizeof(Short_t)*NSAMPLE_V1724);
  memcpy(wave,Data.Probe2,sizeof(Short_t)*NSAMPLE_V1724);
  memcpy(dp1,Data.Dprobe1,sizeof(Short_t)*NSAMPLE_V1724);
  memcpy(dp2,Data.Dprobe2,sizeof(Short_t)*NSAMPLE_V1724);
  
  for(Int_t i=0;i<750;i++){//for even sample
  sample[i] = i;
  dwave[i]=wave[i*2];
  }

  if(det==2 ) {
    if(!BaseLineCorrection()) return 0;
    if(!PuDet()) return 0;
    np=2;
 }
  return 1;
}

bool judge(const pair<Double_t,int> a, const pair<Double_t,int> b)
{
  return a.first>b.first;
}

bool judge_t(const pair<Double_t,int> a, const pair<Double_t,int> b)
{
  return a.second<b.second;
}

bool compare(int a, int b)
{
  return a<b;
}

Bool_t mapping::PuDet()
{     
    
      rt=AmpRatio();
      for(int i=0;i<750;i++) twave[i]=rt*dwave[i];//amplitude normalization
   
      //fast filter for test
      Int_t n=10, gg=1;
      Int_t ll=n;
      gfl=0;
	  int G=gg, L=ll;
	  int i1, i2, j1, j2;
      
      Double_t s0=0,s1=0,tmax=0;
      for(int i=0; i<ll+1; i++) s0 += twave[i];
      for(int i=ll+gg; i<2*ll+gg+1; i++) s1 +=twave[i];
      Int_t gf=0;
      for(int i=0;i<750;i++) {
        if(i>650&&i<700&&dwave[i]<40) gf++;

        //trapz shape
	  if(i>=2*ll+gg+1 && i<750) {
	  s0 += twave[i-ll-gg]-twave[i-2*ll-gg-1];
	  s1 += twave[i]-twave[i-ll-1];
	  int k=i-ll-1-gg/2;
	  trap[k]=(s1-s0)/(ll+1);
	  if(trap[k]<0)trap[k]=0;
	  if(tmax<trap[k]) tmax=trap[k];
	}
      }
      if(gf>10) gfl=1;
          
     //iterator element
      vector<pair<Double_t,int> > vpeak;
      for(int i=n+1;i<750-n-1;i++) {
    if(trap[i]<=0) continue;
    int bp=0;
    for(int j=1;j<3;j++)
      if(trap[i]>=trap[i+j] && trap[i]>trap[i-j])
        bp++;
    
    if(bp==2) vpeak.push_back(make_pair(trap[i],i));
      }
      
      // another parameters
      n=10, gg=1;
      ll=n;
      // initialize trap
      vector<pair<Double_t,int> > vpeak1;
      tmax=0;

      G=0;
      L=11;//L=ll+1 for triangle shape
      for(int i=0;i<750;i++){
        s0=0,s1=0;
        i1 = i-2*L-G+1;
        j1 = i-L+1;
        i2 = i-L-G;
        j2 = i;
        int k=i-L-G/2;
        if(i1>=0&&j1>=0&&i2<750&&j2<750){
          for(int a = i1;a<=i2;a++) s0+=twave[a];
          for(int b = j1;b<=j2;b++) s1+=twave[b];
          trap1[k]=(s1-s0)/L;
        }
        else trap1[k]=0;
      }
      
      //iterator element
      for(int i=n+1;i<750-n-1;i++) {
	if(trap1[i]<=0) continue;
        int bp=0;
	for(int j=1;j<3;j++)
	  if(trap1[i]>=trap1[i+j] && trap1[i]>trap1[i-j])
	    bp++;
	if(bp==2) vpeak1.push_back(make_pair(trap1[i],i));
      }
      
       //select largest amplitude of trap
      sort(vpeak1.begin(),vpeak1.end(),judge);
      
      if(vpeak1.size()>3)
        thr2=vpeak1[3].first;
        else thr2=3.;
      
      if(vpeak1.size()>2)
        vpeak1.erase(vpeak1.begin()+2,vpeak1.end());//remove redundant element 

      sort(vpeak1.begin(),vpeak1.end(),judge_t);//timestamp sort
      Int_t np=0;
      for(int i=0;i<int(vpeak1.size());i++) {
        msd[i]=vpeak1[i].first;
        tsp[i]=vpeak1[i].second;
        np++;
      }
      // cout<<np1<<endl;
      vpeak1.clear();
      
        //fluctuation estimation
        Int_t ta=tsp[np-1]+1,tb=740;
        flu=10.;
        
        if(ta<730)
        if(np>1){
        Int_t mm=2;
        Int_t dt=tb-ta;
        for(int i=1;i<np;i++) 
         if(dt<tsp[i]-tsp[i-1]) {
           ta=tsp[i-1]+mm;
           tb=ta+5;
           dt=tsp[i]-tsp[i-1];
          }
         }
          else if(np==1){
              ta=tsp[0]+2;
              tb=ta+5;
              }
              else{ 
         	ta=150;
         	tb=ta+5;
         	}
           else{
             ta=734;
             tb=739;
           }
         	
        for(int i=ta;i<tb;i++)
          flu+=abs(dwave[i+1]-dwave[i])/(tb-ta);
    //amplitude estimation and filter
     double h0[10],h1[10],hend=0;
     for(int i=745;i<750;i++) hend+=dwave[i]/5.;
     h0[np]=hend;
     Int_t anp=np;
     np=0;
    for(int i=0;i<anp;i++) {
      h0[i]=0;
      h1[i]=0;
      double af,bf;
      af=0.;
      bf=0.;
      int r=tsp[i]+5;
      if(tsp[i]>742) r=749;
      for(int j=tsp[i]-4;j<tsp[i];j++)
        af+=dwave[j]/4.;
      for(int j=tsp[i]+1;j<r;j++)
        bf+=dwave[j]/(r-tsp[i]-1);
        
      int n0 = tsp[i]-2;
      int n1 = tsp[i]+4;
      
      if (i>0&&np>0) n0 = (tsp[np-1]+2 > tsp[i]-2) ? (tsp[np-1]+1) : (tsp[i]-3);
      if (i>0 && i+1<anp) n1 = (tsp[i+1]-1 < tsp[i]+4) ? (tsp[i+1]-1) : (tsp[i]+4);
      for(int j=n0;j<tsp[i]-1;j++) h0[i] += dwave[j] / (tsp[i]-1-n0);
      for(int j=tsp[i]+2;j<n1;j++) h1[i] += dwave[j] / (n1-tsp[i]-2);
      
      if(i==0&&anp>1&&tsp[anp-1]-tsp[i]==3) h1[i]=dwave[tsp[i]+1];
      
      if(h0[i]<0) h0[i]=0;
       
      eh[i]=h1[i]-h0[i];

      tsp[np]=tsp[i];
      eh[np]=eh[i];
      msd[np]=msd[i];
      np++;
    }
       if(np>=0&&np<5)
         return 1;
       else
         return 0;
}

Double_t mapping::AmpRatio()
{  
  for(Int_t i=120;i<140;i++)
    if(dwave[i]<min) min=dwave[i];
  
  if(det==2&&fr==2) {
    if(str>15&&str<32) return 55./14.5;
    if(str>111&&str<128) return 55./17.;
    return 55./19.5;
  }

  if(det==2&&fr==1) {
    if(str==0 || (str>31&&str<48)) return 55./13.;
    return 55./17.5;
  }
  //average amplitude normalie to 3000.
  /*
  Double_t ave=0;
  for(Int_t i=150;i<750;i++)
    ave+=dwave[i]/600.;*/
  
  //return 3000./ave;
    return 1.;
}
  TGraph *gra1 = new TGraph;

//perform the baseline correction and determined trigger positon from dp1
bool mapping::BaseLineCorrection()
{
     // for(Int_t i=0;i<750;i++) cout<<i*2<<" "<<dwave[i*2]<<endl;
     Int_t i0=trig0-10;
     Double_t ave=0;

     for(Int_t j=0;j<i0;j++) ave+= dwave[j];
     ave /=i0;
     npu=0;
     for(Int_t i=0;i<750;i++) {
       dwave[i]-=ave;
       if(dp1[i]>0) npu++;
    //   if(i>trig0+50 && dwave[i]<thr()) return 0; //abnormal fluctuation due to noise.
     }
     for(Int_t i=0;i<trig0-10;i++)
     dwave[i]=0;//for superpulse fitting
     
     return 1;
}