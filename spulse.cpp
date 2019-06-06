#define spulse_cxx
#include "spulse.h"
#include <TPad.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH1.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TString.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
using namespace std;

 Double_t xpar[128][2];
 Double_t ypar[48][2];
 Double_t sl=2.56038,cept=-33.8167;

void spulse::Readpar(){
  ifstream in("EF_equal.dat");
  if(!in.is_open()){
         cout<<"cannot open file !"<<endl;
         for(int i=0;i<128;i++){
                xpar[i][0]=1;
                xpar[i][0]=0;
        }
         for(int i=0;i<48;i++){
                ypar[i][0]=1;
                ypar[i][0]=0;
        }
         return;
  }

  Double_t xp,yp;
  Int_t ch;
  while(1){
    in>>ch>>xp>>yp;
        if(!in.good()) break;
        if(ch<128) {
          xpar[ch][0]=xp;
          xpar[ch][1]=yp;
        } else{
          ypar[ch-128][0]=xp;
          ypar[ch-128][1]=yp;
        }
  }
}

spulse::spulse(TTree *tree, Int_t ifn) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("map00315.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("map00315.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
   fn=ifn;
}

spulse::~spulse()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t spulse::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t spulse::LoadTree(Long64_t entry)
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

void spulse::Init(TTree *tree)
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

   fChain->SetBranchAddress("det", &det, &b_det);
   fChain->SetBranchAddress("fr", &fr, &b_fr);
   fChain->SetBranchAddress("str", &str, &b_str);
   fChain->SetBranchAddress("ts", &ts, &b_ts);
   fChain->SetBranchAddress("er", &er, &b_er);
   fChain->SetBranchAddress("e", &e, &b_e);
   
   fChain->SetBranchAddress("np", &np, &b_np);
   fChain->SetBranchAddress("tsp", tsp, &b_tsp);
   fChain->SetBranchAddress("eh", eh, &b_eh);
   fChain->SetBranchAddress("msd", msd, &b_msd);
   fChain->SetBranchAddress("dwave", dwave, &b_dwave);
   fChain->SetBranchAddress("rwave", rwave, &b_rwave);
   //fChain->SetBranchAddress("cr", cr, &b_cr);
   fChain->SetBranchAddress("sample", sample, &b_sample);
   fChain->SetBranchAddress("bd", &bd, &b_bd);
   fChain->SetBranchAddress("ch", &ch, &b_ch);
   fChain->SetBranchAddress("nevt", &nevt, &b_nevt);
   fChain->SetBranchAddress("flu",&flu,&b_flu);
   fChain->SetBranchAddress("st",&st,&b_st);
   fChain->SetBranchAddress("flu",&flu,&b_flu);
   fChain->SetBranchAddress("trap",&trap,&b_trap);
   fChain->SetBranchAddress("trap1",&trap1,&b_trap1);
   fChain->SetBranchAddress("gfl", &gfl, &b_gfl);
   //  fChain->SetBranchAddress("off2", &off2, &b_off2);
   Notify();
}

Bool_t spulse::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is stad when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void spulse::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

Int_t spulse::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}


Double_t esx[128],esy[48];
Double_t tsx[128],tsy[48];
Double_t sx[128][750],sy[48][750];
Double_t sp[750];


void spulse::ReadSuperPulse(TString spulse)
{
 TH1D *hx[128],*htx,*hehx;
 TH1D *hy[48],*hty,*hehy;
  TFile *falpha= new TFile(spulse);
  TString sh;
  //  int t0=261;
  double e,t;
  htx=(TH1D*)gDirectory->Get("htx");
  hty=(TH1D*)gDirectory->Get("hty");
  hehx=(TH1D*)gDirectory->Get("hehx");
  hehy=(TH1D*)gDirectory->Get("hehy");
  for(int i=0;i<128;i++) {
    sh.Form("hx%03d",i);
    //    cout<<sh<<endl;
    hx[i]=(TH1D*)gDirectory->Get(sh);
    for(int j=0;j<740;j++) sx[i][j]=hx[i]->GetBinContent(j+1);
    tsx[i]=htx->GetBinContent(i+1);
    esx[i]=hehx->GetBinContent(i+1);
	//cout<<i<<" "<<esx[i]<<" "<<tsx[i]<<endl;
    
    if(i<48) {
      sh.Form("hy%03d",i);
      hy[i]=(TH1D*)gDirectory->Get(sh);
      for(int j=0;j<740;j++) sy[i][j]=hy[i]->GetBinContent(j+1);
      tsy[i]=hty->GetBinContent(i+1);
      esy[i]=hehy->GetBinContent(i+1);
	//cout<<i<<" "<<esy[i]<<" "<<tsy[i]<<endl;
    }
  }
}


TH1D *h;//=new TH1D("h","h",1500,0,1500);
Double_t func(Double_t *x, Double_t *par)
{
  Double_t yy=0;
  Double_t ss=0;
  Int_t j;

  for(int i=0;i<3;i++) {
    j=Int_t(x[0])-Int_t(par[2*i+1]);// j=Int_t(x[0]-par[2*i+1]), it will cause some problems;
    if(j<0 || j>739) continue;
    else ss=sp[j];
    yy+=par[i*2]*ss;    
  }
  return yy;
}

//fit with baseline
Double_t func1(Double_t *x, Double_t *par)
{
  Double_t yy=0;
  Double_t ss=0;
  Int_t j;

  for(int i=0;i<3;i++) {
    j=Int_t(x[0])-Int_t(par[i*2+4]);// j=Int_t(x[0]-par[2*i+1]), it will cause some problems;
    if(j<0 || j>739) continue;
    else ss=sp[j];
    yy+=par[i*2+3]*ss;    
  }
  return par[0]+par[1]*exp(par[2]*x[0])+yy;
}

Double_t func2(Double_t *x, Double_t *par){
  return par[0]+par[1]*exp(par[2]*x[0]);
}

void spulse::BranchOpt(TTree *opt) {

  opt->Branch("det",&det,"det/I");
  opt->Branch("fr",&fr,"fr/I");
  opt->Branch("str",&str,"str/I");
  opt->Branch("ts",&ts,"ts/l");
  opt->Branch("er",&er,"er/D");
  opt->Branch("e",&e,"e/D");
  opt->Branch("fe",&fe,"fe/D"); 
  opt->Branch("ip",&ip,"ip/I");
  opt->Branch("np",&np,"np/I");
  //opt->Branch("ew",ew,"ew[np]/D");
  opt->Branch("ep",ep,"ep[np]/D");
  opt->Branch("fep",fep,"fep[np]/D");
  opt->Branch("epw",epw,"epw[np]/D");
  opt->Branch("tsp",tsp,"tsp[np]/l");
  opt->Branch("tsw",&tsw,"tsw/I");
  opt->Branch("dwave",dwave,"dwave[750]/D");
  opt->Branch("rwave",rwave,"rwave[750]/S");
  opt->Branch("sample",sample,"sample[750]/S");
  opt->Branch("chi2",chi2,"chi2[np]/D");
  opt->Branch("trapezoidE",&trapezoidE,"trapezoidE/D");
  opt->Branch("chi",&chi,"chi/D");
  opt->Branch("para",para,"para[10]/D");
  opt->Branch("ndf",&ndf,"ndf/I");
  opt->Branch("flu",&flu,"flu/D");
  opt->Branch("chi2h",&chi2h,"chi2h/D");
  opt->Branch("flu",&flu,"flu/D");
  opt->Branch("fwave",fwave,"fwave[730]/D");
  
  opt->Branch("trap",trap,"trap[750]/D");
  opt->Branch("trap1",trap1,"trap1[750]/D");
}

void spulse::Reset()
{
  memset(para,0,sizeof(para));
  memset(ew,0,sizeof(ew));
  memset(tsp,0,sizeof(tsp));
  memset(tsp1,0,sizeof(tsp1));
  memset(chi2,0,sizeof(chi2));
  memset(eph,0,sizeof(eph));
  memset(fslope,0,sizeof(fslope));
  np=0;
  ip=-1;
  tsw=0;
}

Bool_t reject;
Double_t fitf(Double_t *x,Double_t *par) {
  if (reject && (x[0]<Int_t(par[4]) || (x[0]>745 && x[0]<Int_t(par[2])) || x[0]>Int_t(par[3]))) {
    TF1::RejectPoint();
    return 0;
  }
  Double_t fitval = par[0]*exp(-par[1]*(x[0]-Int_t(par[4])));
  return fitval;
}


TCanvas *c1=new TCanvas("c1","c1",500,500);
TF1 *f;

TGraph *gr=new TGraph;
TH1D *h1 = new TH1D("h1","h1",750,0,750);
TH1D *h2 = new TH1D("h2","h2",750,0,750);  
TH1D *gr1 = new TH1D("h",",h",750,0,750);
void spulse::Loop(TTree *opt)
{
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast()*0.1;
   cout<<nentries<<endl;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      Reset();
      fChain->GetEntry(jentry);   
      if(det!=2 ) {
        opt->Fill();
        continue;
      }
      if(np==0||gfl==1) continue;
      bool t=false;      

      Double_t hsmax=0;
      int tpulse=0;
      //x,128
      memset(sp,0,sizeof(sp));
      if(fr==1) {
	memcpy(sp,sy[str],sizeof(Double_t)*740);
        hsmax=esy[str];
        tpulse=tsy[str];
        }
      
      if(fr==2) {
	memcpy(sp,sx[str],sizeof(Double_t)*740);
        hsmax=esx[str];
        tpulse=tsx[str];
        }
      Double_t a[3],off[3]; 
      Double_t a1[3],off1[3];
      Double_t tpar[5];
      Double_t fpar[9];
      Double_t tspar[4];
      Double_t cpar[3];
      ULong64_t ftsp[3];
      memset(cpar,0,sizeof(cpar));
      memset(fpar,0,sizeof(fpar));
      for(int i=0;i<3;i++) {
        a[i]=0;
        off[i]=0;
        }
      
      for(int i=0;i<np;i++) {
        a[i]=eh[i]/hsmax;
        off[i]=double(tsp[i])-tpulse;
        ftsp[i]=tsp[i];
        }
      
      /***************************** two step fit ***********************************/
      c1->Clear();
      gr->Clear();
      for(int i=0;i<750;i++) gr->SetPoint(i,i,dwave[i]);
      
      //Two step fit, Get the first parameters for np==2
      if(np>1){
      f= new TF1("f",func,tsp[0]+1,tsp[1]-2,2);
      f->SetParameters(a[0],off[0]);
      f->SetParLimits(0,a[0]*0.20,a[0]*5.0);
      gr->Fit("f","RQ");
      f->GetParameters(tpar);
      a[0]=tpar[0];
      off[0]=tpar[1];
      chi2[0]=f->GetChisquare()/f->GetNDF();
      }
      
      Int_t lr,hr;
      lr=120;

      if(tsp[np-1]>700&&tsp[np-1]<735) hr=tsp[np-1]+5;
      else if(tsp[np-1]>=735) hr=739;
      else hr=700;//735,55
      
      f = new TF1("f",func,lr,hr,6);
      f->SetParameters(a[0],off[0],a[1],off[1],a[2],off[2]);
      if(np==1) {
        f->SetParLimits(0,0.2*a[0],5.0*a[0]);
        f->FixParameter(1,off[0]);
        f->FixParameter(2,0);
        f->FixParameter(3,0);
        f->FixParameter(4,0);
        f->FixParameter(5,0);
        }
      
      if(np==2){
        f->FixParameter(1,off[0]);
        f->SetParLimits(2,0.2*a[1],5.0*a[1]);
        f->FixParameter(3,off[1]); 
        f->FixParameter(4,0);
        f->FixParameter(5,0); 
        }
     
      TFitResultPtr r = gr->Fit(f, "RQS");
      for(Int_t i=0;i<730;i++)fwave[i]=f->Eval(i);

      ndf=f->GetNDF();
      
      f->GetParameters(&para[0]);
      
      if(np==1)chi2[0]=f->GetChisquare()/f->GetNDF();
      if(np==2)chi2[1]=f->GetChisquare()/f->GetNDF();
      if(np==3)chi2[2]=f->GetChisquare()/f->GetNDF();
      delete f;
      
       for(int i=0;i<np;i++){
       ew[i]=para[2*i]*2000.;
       epw[i]=parah[2*i]*2000.;
       tsp1[0]=ts;
       if(i>0)tsp1[i]=para[2*i+1]*20.-para[1]*20.+ts;
        
       if(fr==1) ew[i]=ew[i]*ypar[str][0]+ypar[str][1];
       if(fr==2) ew[i]=ew[i]*xpar[str][0]+xpar[str][1];
       ew[i]=ew[i]*sl+cept;
       
       if(fr==1) epw[i]=epw[i]*ypar[str][0]+ypar[str][1];
       if(fr==2) epw[i]=epw[i]*xpar[str][0]+xpar[str][1];
       epw[i]=epw[i]*sl+cept;
       }
       
       for(int j=0;j<np;j++) {
          tsp[j]=tsp1[j];
          ep[j]=ew[j];
       }
      
      er=e;       
      for(int i=0;i<np;i++){
          ts=tsp1[i];
          e=ep[i];
          fe=fep[i];
          ip=i;
          chi=chi2[i];
          if(i==0)trapezoidE=er;
            else trapezoidE=0;
            
       opt->Fill();
     }
      if(jentry%1000==0) cout<<int(jentry*100./nentries)<<" %  ";
   }
}

