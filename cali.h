
Int_t xhit,yhit;
Int_t x[10],y[10];
Double_t xe[10],ye[10],mwe,ve[3];
ULong64_t xt[10],yt[0],mwt,vt[3];

void ResetOpt(){
  xhit=0;
  yhit=0;
  memset(x,0,sizeof(x));
  memset(y,0,sizeof(y));
  memset(xe,0,sizeof(xe));
  memset(ye,0,sizeof(ye));
  memset(xt,0,sizeof(xt));
  memset(yt,0,sizeof(yt));
  memset(ve,0,sizeof(ve));
  memset(vt,0,sizeof(vt));
  mwe=0;
  mwt=0;
}
