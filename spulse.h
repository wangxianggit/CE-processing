//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jul  1 15:00:03 2016 by ROOT version 5.34/36
// from TTree tree/tree
// found on file: map00315.root
//////////////////////////////////////////////////////////

#ifndef spulse_h
#define spulse_h
#include <TString.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class spulse {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           det;
   Int_t           fr;
   Int_t           str;
   ULong64_t       ts;
   Double_t        er;
   Double_t        e;
   Double_t        fe; //free parameters fit
   Int_t           np;
   Double_t        ep[10];   //[np]
   Double_t        fep[10];  //free parameters fit
   Double_t        eph[10];
   ULong64_t       tsp[10];   //[np]
   ULong64_t       tsp1[10];
   Double_t        eh[10];   //[np]
   Double_t        msd[10];
   Double_t        sd[750];
   Double_t        fd[750];
   Short_t         cr[750];
   Double_t        dwave[750];
   Double_t        fwave[730];//fit range
   UShort_t        sample[750];
   UShort_t        rwave[1500];
   Double_t        bwave[750];
   Double_t        d1wave[750];
   Int_t           npt;
   Int_t           bd;
   UChar_t         ch;
   Int_t           nevt;
   Double_t        chi2[10];
   Double_t        chi2h;//for np==1;
   Int_t           ip;
   Int_t           ndf;
   Double_t        flu;
   Bool_t          st;
   Double_t        eph1[10];
   Int_t           tsw;
   Int_t           gfl;
   Int_t           fn;
   
   //trapezoid...
   Double_t        trapezoidE;
   Double_t        chi;
   Double_t        trap[750];
   Double_t        trap1[750];
   
   //new branch
   Double_t        ew[10];
   Double_t        para[10];
   Double_t        parah[10];
   Double_t        pa[10];
   Double_t        fslope[750];
   Double_t        dif;
   Double_t        epw[10];
   //   Double_t        a2;
   //Int_t           off2;

   // List of branches
   TBranch        *b_det;   //!
   TBranch        *b_fr;   //!
   TBranch        *b_str;   //!
   TBranch        *b_ts;   //!
   TBranch        *b_er;   //!
   TBranch        *b_e;   //!
   TBranch        *b_te;
   TBranch        *b_rte;
   TBranch        *b_cate;
   TBranch        *b_slowfilter;
   TBranch        *b_np;   //!
   TBranch        *b_tsp;   //!
   TBranch        *b_eh;   //!
   TBranch        *b_msd;   //!
   TBranch        *b_dwave;   //!
   TBranch        *b_rwave;
   TBranch        *b_sd;   //!
   TBranch        *b_fd;
   TBranch        *b_flu;
   TBranch        *b_trap;
   TBranch        *b_trap1;
   //   TBranch        *b_cr;   //!
   TBranch        *b_sample;   //!
   TBranch        *b_bd;   //!
   TBranch        *b_ch;   //!
   TBranch        *b_nevt;   //!
   TBranch        *b_slope;
   TBranch        *b_st;
   TBranch        *b_gfl;   //!
   //  TBranch        *b_off2;   //!

   spulse(TTree *tree=0, Int_t ifn=0);
   virtual ~spulse();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Readpar();
   virtual void     Loop(TTree *tree);
   virtual void     BranchOpt(TTree *tree);
   virtual void     ReadSuperPulse(TString spulse);
   virtual Bool_t   Notify();
   virtual void     Reset();
   virtual void     Show(Long64_t entry = -1);
};

#endif

