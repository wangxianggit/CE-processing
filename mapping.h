//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Apr 24 12:04:05 2016 by ROOT version 5.34/21
// from TTree DGTZ00/A Tree for Board 00
// found on file: ../../run00307_wave.root
//////////////////////////////////////////////////////////

#ifndef mapping_h
#define mapping_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TCanvas.h>
#include "setup.h"

//#define TEST
//#define HIST
// Header file for the classes stored in the TTree if any.
 struct DPP_PHA_Event_V1724{
	       UShort_t TT;
	       UShort_t PU;
	       UShort_t Extras;
	       UShort_t Energy;
	       ULong64_t TriggerTimeTag;
	       UShort_t Dprobe1[NSAMPLE_V1724];
	       UShort_t Dprobe2[NSAMPLE_V1724];
	       Short_t Probe1[NSAMPLE_V1724];
	       Short_t Probe2[NSAMPLE_V1724];
  };
  
// Fixed size dimensions of array or collections stored in the TTree if any.

class mapping {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          BankNo;
   UChar_t         ChannelNo;
   UInt_t          SysTime;
   DPP_PHA_Event_V1724  Data;
   UShort_t        Sample[NSAMPLE_V1724];

   // List of branches
   TBranch        *b_BankNo;   //!
   TBranch        *b_ChannelNo;   //!
   TBranch        *b_SysTime;   //!
   TBranch        *b_Data;   //!
   TBranch        *b_Sample;   //!
  
   
   double ratio[2];//trapezoid cali
  
   // output variables
   Int_t		det; // 1 mwpc 2 dssd 3 veto 4 ge 
   Int_t		fr; // 1 front 2 back
   Int_t		str;
   ULong64_t	        ts;
   Double_t		e;//calibrated
   Double_t             er;//raw
   Short_t		wave[NSAMPLE_V1724];
   Short_t              rwave[NSAMPLE_V1724];
   Short_t              dp1[NSAMPLE_V1724];
   Short_t              dp2[NSAMPLE_V1724];
   Double_t             dwave[750];
   Double_t             twave[750];
   Double_t             trap[750];
   Double_t             trap1[750];
   Double_t             rt;
   //graph
   
   Int_t		boardn;
   Int_t                nevt;
   Int_t                np;
   Int_t                np1;
   Int_t                npu;
   Short_t              pu;
   Double_t             flu;
   Short_t              sample[750];
   ULong64_t            tsp[10];
   ULong64_t            tsp1[10];
   Double_t             thr1;
   Double_t             thr2;
   Double_t             min;
   
   Double_t             eh[10];
   Double_t             msd[10];
   Double_t             msd1[10];
   Float_t              ene_thr;
   Int_t                trig0;
   Double_t             a2;
   Int_t                off2;
   Double_t             bslope;
   Double_t             slope;
   Double_t             slopechi2;
   Bool_t               st;//saturation
   // Double_t             fluc;
   
   TBranch        *b_er;
   TBranch        *b_rte;
   TBranch        *b_str;
   TBranch        *b_np;
   TBranch        *b_det;
   TBranch        *b_fr;
   Int_t          number;
   Int_t          gfl;//Flag for noise

   mapping(TTree *tree=0);
   virtual ~mapping();
   virtual void     Readpar();
   virtual void     SetHist();
   virtual void     FillHist();
   virtual void     WriteHist(int runnum);
   virtual bool     BaseLineCorrection();
   virtual Bool_t     PuDet();
   virtual Double_t AmpRatio();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TTree *ipt=0,Int_t bdn=0,TTree *opt=0);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
//   virtual void     TrapezoidStd();
   Bool_t 
   BranchOpt(TTree *opt);
   void ResetOpt();
   Bool_t Mapping();
};

#endif
