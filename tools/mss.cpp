/*
  MAGIC Source Simulator. The macro takes a given spectral shape, 
  and  using MAGIC performance in the range 40GeV-16TeV from 
  Aleksic, J., et al., 2016, Astroparticle Physics, 72, 76 
  (please cite this paper if you want to use this macro in a publication)
  estimates what kind of signal can be observed by MAGIC.
  It computes significances of each spectral point according to
  Li & Ma 1983, ApJ, 272, 317, Eq. 17 definition.

  Terminal Output:
  For each estimated energy range one gets the number of excess events,
  signal-to-background ratio, significance, and information if a given
  bin satisfies the conditions for the detection. 

  Plot Output:
  Spectral points following the assumed spectrum are shown, their 
  error bars reflect the performance of the instrument for such a source.
  Significance for each bin is given (bins without detection
  have gray numbers. In the top of the plot a simple number using 
  information from all the shown bins is given to evaluate 
  if the source can be detected (roughly if sum of sigma_i / sqrt(N) >~5) 

  requirements:
  - ROOT (https://root.cern.ch), the macro was tested on 5.34/34 version
    and on the newest ROOT 6.28/00

  caveats:
  - the macro is operating on estimated energy doing simple
  comparisons with differential in estimated energy rates seen for
  Crab Nebula, for soft sources the differences in energy migration 
  will result in different performance than calculated here
  - the threatment of extended sources is very approximate, only
  increase in background is taken into account (without energy
  dependence of PSF) for extended sources for extension >~0.4 deg the
  dependence on the offset from the centre of the camera will further
  worsen the performance w.r.t. the one calculated here
  - significances are given for each differential energy bin
  separatelly, but to detect a source one normally applies a cut that
  keeps a broad range of energies inside resulting in better integral
  sensitivity than differential one. Also, optimization of cuts for 
  a broad energy range usually results in somewhat better sensitivity
  than what one can get by simply integrating the used here signal in 
  differential energy bins. As a very crude approximation for detection
  capability we calculate here also a sum of significances of all the
  points in SED divided by the sqrt of the number of those points

  usage:
  - modify the "basic settings" below and run in terminal:
  root mss.cpp++

  for questions ask Julian Sitarek (jsitarek [at] uni.lodz.pl)

  version history:
  1.9 (2024/10/25): Added high zenith performance for MAGIC+LST1 and MAGIC, fixed the combination of sensitivities
  1.8 (2024/01/29): Changed the calculation of combined sensitivity
  1.7 (2023/11/08): updated the LST1+MAGIC mid zenith sensitivity to the one in the paper
  1.6 (2023/11/04): fix for ROOT6
  1.5 (2022/10/24): again updated MAGIC+LST1 performance
  1.4 (2022/09/19): updated MAGIC+LST1 performance
  1.3 (2021/10/13): added also MAGIC+LST1 performance from MC study presented in ICRC 2019
  1.2 (2020/07/24): a number of additional features:
      - EBL attenuation with Dominguez+11 model
      - SUMT performance (only for low zenith)
      - pulsar mode
  1.1 (2017/10/16): fixed a bug which made the errorbars too small
 */

#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TGraph2D.h>
#include <TH1.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMath.h>

#include <iostream>
#include <fstream>

using namespace std;

// basic settings
const Double_t timeh=20.; //[h], time of observations
const Double_t extension=0.; // [deg] extension radius of the source

// Source settings
// [cm^-2 s^-1 TeV^-1] assumed (intrinsic) source spectrum (dN/dS dt dE) as a function of energy in GeV
// the string understands many mathematical expressions, some common examples below  
// const TString assumedstr="2.7e-9 * pow(x/100, -3.18)"; // simple power-law
const TString assumedstr="2.e-11 * pow(x/1000, -1.99) * exp(-x/100)"; // power-law with an exponential cut-off
// const TString assumedstr="1.e-11 * pow(x/1000, -1.9) * exp(-pow(x/300, 0.5))"; // power-law with a non-exponential cut-off
// const TString assumedstr="0.007* 3.39e-11 * pow(x/1000, -2.51-0.21*log10(x/1000))"; //Crab Nebula [Aleksic et al. 2016]
// const TString assumedstr="0.1* exp(-x/200)* 3.39e-11 * pow(x/1000, -2.51-0.21*log10(x/1000))"; //10% Crab with 200 GeV cut-off
const Double_t redshift=-1; // redshift of the source (for the EBL absorption), if -1 then no absorption
const TString srcname="PWL+cutoff"; // name of the source

// const Double_t redshift=0.94; 
// const TString assumedstr="2.7e-9 * pow(x/100, -3.18)"; // simple power-law
// const TString srcname="PKS1441"; // name of the source

// const TString assumedstr="2.28e-9 * pow(x/32.1, -5.62)"; // simple power-law
// const TString srcname="Geminga"; // name of the source

// const TString assumedstr="4.e-9 * pow(x/100., -2.6)"; // simple power-law
// const TString srcname="PG 1553"; // name of the source

// const TString assumedstr="1.4e-12 * pow(x/1310., -2.74)"; // simple power-law
// const TString srcname="J1835-069"; // name of the source

// const TString assumedstr="7.9e-11 * pow(x/175, -3.26)"; // simple power-law
// const Double_t redshift=0.360; // redshift of the source (for the EBL absorption), if -1 then no absorption
// const TString srcname="PKS1510low"; // name of the source

enum Zeniths {lowzd, midzd, hizd};
// const Bool_t ismidzd=kFALSE; // false: performance for 0-30deg zenith angles, true: for 30-45deg zenith
Zeniths zenith=lowzd; // performance for 0-30deg zenith angles
// Zeniths zenith=midzd; // performance for 30-45deg zenith angles
// Zeniths zenith=hizd; // performance for ~60deg zenith angles

// you can check visibility of your source e.g. here: http://www.magic.iac.es/scheduler/
const Bool_t isSUMT=kFALSE;

// advanced settings
const Int_t numoff=3; // number of background estimation regions
const Double_t minev=10; // minimum number of events
const Double_t minSBR=0.05; // minimum ratio of excess to background
const Double_t PSF=0.1; //[deg] PSF (for worsening the performance for extended sources)

const Double_t offsetdegrad=1.; // degradation factor (applied to gamma and bgd. rates) if observations are taken at higher offset then 0.4 deg from camera center
// for best sensitivity region (~300 GeV at low zenith) and MAGIC-only or MAGIC+LST1 observations it can be estimated as 1.1*exp(-0.8 * offset^2)
// for proposals focussed on multi-TeV emission a slower drop can be assumed: ~exp(-0.3*offset^2)

const TString crabstr="3.39e-11 * pow(x/1000, -2.51-0.21*log10(x/1000))"; //Crab Nebula [Aleksic et al. 2016]
const Double_t eplotmin=30; //[GeV] x plot range
const Double_t eplotmax=20.e3; //[GeV] x plot range
const Double_t yplotmin=1.e-14; // [TeV cm^-2 s^-1] y plot range
const Double_t yplotmax=1.e-9; //[TeV cm^-2 s^-1] y plot range
const Double_t minerror=2; // showing only points with value > minerror * error
const Bool_t drawsigma=kTRUE; // whether to draw also sigmas on the plot

// pulsar mode settings 
const Bool_t pulsarmode=kFALSE; // if true the background is reduced to on phase (see below) and SBR cut is ignored
const Double_t pulsarOnRange=0.092;// range of ON phases used for pulsar mode
const Double_t pulsarOffRange=0.25;// range of OFF phases used for pulsar mode

const Bool_t isLSTmode=kFALSE; // uses MAGIC+LST1 (software trigger) performance (data in mid zenith, MC in low zenith)


// global variables (DO NOT MODIFY)
const Int_t npoints = 13;
Double_t crabrate[npoints];// [min^-1]
Double_t bgdrate[npoints];// [min^-1]
Double_t enbins[npoints+1]; // [GeV]
TGraph *gEBL=0; // tau vs log E for interpolation
TString pathebl="./dominguez_ebl_tau.txt"; // path with EBL model of Dominguez+11
TF1 *ffluxint=0; // intrinsic source flux
const TString version="1.9";
Bool_t ismc=kFALSE; // warning message if the performance curve is MC-based

const char *const conError  = "\033[31m";
const char *const conWarn   = "\033[33m\033[1m";
const char *const conReset  = "\033[0m";
//=======================================================================================
Double_t fluxobs(Double_t *xx, Double_t *par) // observed flux dN/dE
{
  Double_t en = xx[0]; // GeV
  Double_t tau = gEBL->Eval(log10(en)); // simple linear interpolation in log E
  Double_t atten = exp(- tau);
  return ffluxint->Eval(en) * atten;
}

//=======================================================================================
Double_t sedobs(Double_t *xx, Double_t *par) // observed sed E^2 dN/dE
{
  Double_t en = xx[0]*1.e-3; // GeV ==>TeV
  Double_t f = fluxobs(xx, par); // TeV^-1 cm^-2 s^-1
  return f * en *en; // TeV cm^-2 s^-1
}

//=======================================================================================
Bool_t LoadEBL(Double_t z)
{
  gEBL = new TGraph();
  ifstream tabfile (pathebl.Data());
  Int_t nz;
  const Int_t maxnbinsz=40;
  Double_t zz[maxnbinsz], en, tau;
  tabfile>>nz;
  for (int i=0; i<nz; i++)
    tabfile>>zz[i];
  // cout<<nz<<" values of z: ";
  // for (int i=0; i<nz; i++)
  //   cout<<zz[i]<<", ";
  // cout<<endl;
  
  TGraph2D *grebl = new TGraph2D();

  while (!tabfile.eof())
    {
      tabfile>>en;
      en*=1.e3; // GeV
      if (tabfile.eof())
	break;

      // cout<<"new energy: "<<en<<" TeV"<<endl;
      for (int i=0; i<nz; i++)
	{
	  tabfile>>tau;
	  grebl->SetPoint(grebl->GetN(), log10(en), zz[i], tau);
	}
    }

  tabfile.close();
  cout<<nz<<" redshifts readout in EBL attenuation "<<grebl->GetN()<<" total points"<<endl; 
  // TCanvas *c3d = new TCanvas ("c3c", "", 640, 480);
  // c3d->SetLogz();
  // grebl->GetXaxis()->SetTitle("log10(E/TeV)");
  // grebl->GetYaxis()->SetTitle("z");
  // grebl->GetZaxis()->SetTitle("tau");
  // grebl->Draw("lego");

    
  Int_t npoint=100; 
  Double_t enmin=10., enmax=29.99e3;
  for (int i=0; i<=npoint; i++)
    {
      en=enmin* pow(enmax/enmin, 1.*i/npoint);
      tau=grebl->Interpolate(log10(en), z);
      gEBL->SetPoint(gEBL->GetN(), log10(en), tau);
    }
  delete grebl;
  return kTRUE;
}

//=======================================================================================
void Prepare()
{
  if (isLSTmode)
    {
      for (int i=0; i<=npoints; i++)
	enbins[i]=100* pow(10, (i-2)*0.2);
      if (zenith==hizd) // MC values values computed at 59 zd angle, might by ~20% too optimistic >~1 TeV
      {
	crabrate[0] = 0; // 50GeV, below threshold
	crabrate[1] = 0. ; // 79 GeV, below threshol
	crabrate[2] = 0. ; // 126 GeV below threshold
	crabrate[3] = 0. ; // 200 GeV, below threshold
	crabrate[4] = 0.6929; // 316 GeV
	crabrate[5] = 2.0531;
	crabrate[6] = 1.5028;
	crabrate[7] = 0.8963;
	crabrate[8] = 0.5915;
	crabrate[9] = 0.3513;
	crabrate[10]= 0.2314;
	crabrate[11]= 0.111; // 7 TeV
	crabrate[12]= 0.0488; // 12.6 TeV

	bgdrate[0] = 0; // below threshold
	bgdrate[1] = 0.;
	bgdrate[2] = 0.;
	bgdrate[3] = 0.; // 200 GeV below threshold
	bgdrate[4] = 0.349;
	bgdrate[5] = 0.2163;
	bgdrate[6] = 0.0393;
	bgdrate[7] = 0.0126;
	bgdrate[8] = 0.0037;
	bgdrate[9] = 0.0013;
	bgdrate[10]= 0.0018;
	bgdrate[11]= 0.0005;
	bgdrate[12]= 0.0006; // 12.6 TeV 
	ismc=kTRUE;

	// the MC sensitivites are often too optimistic. While at low zenith MAGIC+LST1 sensitivity more or less matches
	// between the data and MCs, at mid zenith there is some ~20% difference >~1 TeV. Very likely the same happens
	// at high zenith, so this is why we artificially increase background by 40% to emulate this effect
	// and make more fair MAGIC-only vs MAGIC+LST1 comparisons at high zenith
	for (int i=0; i<npoints; i++)
	  if (sqrt(enbins[i]*enbins[i+1])>1000)
	    bgdrate[i]*=1.4;
      }
      if (zenith==midzd) // values from https://doi.org/10.1051/0004-6361/202346927 
      {
	crabrate[0] = 0; // 50GeV, below threshold
	crabrate[1] = 0.59 ;
	crabrate[2] = 2.43 ;
	crabrate[3] = 2.65 ;
	crabrate[4] = 2.03 ;
	crabrate[5] = 1.171;
	crabrate[6] = 0.899;
	crabrate[7] = 0.806;
	crabrate[8] = 0.319;
	crabrate[9] = 0.185;
	crabrate[10]= 0.113;
	crabrate[11]= 0.084;
	crabrate[12]=0 ; // missing data

	bgdrate[0] = 0; // below threshold
	bgdrate[1] = 0.715;
	bgdrate[2] = 1.42;
	bgdrate[3] = 0.307;
	bgdrate[4] = 0.093;
	bgdrate[5] = 0.0229;
	bgdrate[6] = 0.0093;
	bgdrate[7] = 0.0096;
	bgdrate[8] = 0.00264;
	bgdrate[9] = 0.0007;
	bgdrate[10]= 0.00138;
	bgdrate[11]= 0.00148;
	bgdrate[12]=0 ; // missing data
      }
      else if (zenith==lowzd)
      {
	// those numbers are taken from MC, but they agree
	// with a small bunch of data available 
	crabrate[0] = 0.3795; // 50GeV
	crabrate[1] = 2.5808;  // 80GeV
	crabrate[2] = 2.8257;
	crabrate[3] = 1.6654;
	crabrate[4] = 1.3289;
	crabrate[5] = 1.2871;
	crabrate[6] = 0.8834;
	crabrate[7] = 0.7045;
	crabrate[8] = 0.3908;
	crabrate[9] = 0.1963;
	crabrate[10]= 0.089;
	crabrate[11]= 0.0368;
	crabrate[12]=0 ; // missing data

	bgdrate[0] = 8.8870e-01; 
	bgdrate[1] = 1.6505e+00;
	bgdrate[2] = 5.7688e-01;
	bgdrate[3] = 5.0804e-02;
	bgdrate[4] = 2.0521e-02;
	bgdrate[5] = 2.1718e-02;
	bgdrate[6] = 5.6106e-03;
	bgdrate[7] = 3.9491e-03;
	bgdrate[8] = 3.2053e-03;
	bgdrate[9] = 1.8074e-03;
	bgdrate[10]= 7.3253e-04;
	bgdrate[11]= 1.1352e-05;
	bgdrate[12]=0 ; // missing data
	ismc=kTRUE;
      }
      
      return;
    }
  
  for (int i=0; i<=npoints; i++)
    enbins[i]=100* pow(10, (i-2)*0.2);

  if (isSUMT)
    {
      if (zenith==lowzd)
	{
	  // 2018/19 Crab data analyzed with generic ST0307 SUMT MCs
	  crabrate[0] =1.39684; // ~50GeV 
	  crabrate[1] =3.12657; 
	  crabrate[2] =3.09145; 
	  crabrate[3] =2.40268; 
	  crabrate[4] =1.32915; 
	  crabrate[5] =0.86180;
	  crabrate[6] =0.51666; 
	  crabrate[7] =0.31533; 
	  crabrate[8] =0.16207; 
	  crabrate[9] =0.09279; 
	  crabrate[10]=0.04624; 
	  crabrate[11]=0.02345; 
	  crabrate[12]=0.00874; // ~12TeV

	  bgdrate[0] =3.33321; 
	  bgdrate[1] =3.24046; 
	  bgdrate[2] =1.32361; 
	  bgdrate[3] =0.406588; 
	  bgdrate[4] =0.091944; 
	  bgdrate[5] =0.032226; 
	  bgdrate[6] =0.007277; 
	  bgdrate[7] =0.003123; 
	  bgdrate[8] =0.001487; 
	  bgdrate[9] =0.001464; 
	  bgdrate[10]=0.001231; 
	  bgdrate[11]=0.001152; 
	  bgdrate[12]=0.000957;
	}
      else
	{
	  // to be implemented
	}
    }
  else
    {
      if (zenith==hizd) // Crab results from 2.5 hrs of 55-62 zd data from 2016-2018. Dataset from Juliane van Scherpenberg.
      {
	crabrate[0] = 0; // 50GeV, below threshold
	crabrate[1] = 0. ; // 79 GeV, below threshol
	crabrate[2] = 0. ; // 126 GeV below threshold
	crabrate[3] = 0. ; // 200 GeV, below threshold
	crabrate[4] = 0.503462; // 316 GeV
	crabrate[5] = 1.60232 ;
	crabrate[6] = 2.26558 ;
	crabrate[7] = 0.928094;
	crabrate[8] = 0.698335;
	crabrate[9] = 0.305662;
	crabrate[10]= 0.173859;
	crabrate[11]= 0.083892; // 7 TeV
	crabrate[12]= 0.069938; // 12.6 TeV

	bgdrate[0] = 0; // below threshold
	bgdrate[1] = 0.;
	bgdrate[2] = 0.;
	bgdrate[3] = 0.; // 200 GeV below threshold
	bgdrate[4] = 0.750815;
	bgdrate[5] = 0.597588;
	bgdrate[6] = 0.564753;
	bgdrate[7] = 0.089775;
	bgdrate[8] = 0.0568584;
	bgdrate[9] = 0.00954936;
	bgdrate[10]= 0.00344762;
	bgdrate[11]= 0.00147755;
	bgdrate[12]= 0.00229841; // 12.6 TeV

	ismc=kFALSE;
      }
      else if (zenith==midzd) //Aleksic et al 2016
	{
	  crabrate[0]=0; // below threshold !!!
	  crabrate[1]=0.404836; 
	  crabrate[2]=3.17608; 
	  crabrate[3]=2.67108; 
	  crabrate[4]=2.86307; 
	  crabrate[5]=1.76124;
	  crabrate[6]=1.43988; 
	  crabrate[7]=0.944385; 
	  crabrate[8]=0.673335; 
	  crabrate[9]=0.316263; 
	  crabrate[10]=0.200331; 
	  crabrate[11]=0.0991222; 
	  crabrate[12]=0.0289831;

	  bgdrate[0]=1.67777; // below threashold
	  bgdrate[1]=2.91732; 
	  bgdrate[2]=2.89228; 
	  bgdrate[3]=0.542563; 
	  bgdrate[4]=0.30467; 
	  bgdrate[5]=0.0876449; 
	  bgdrate[6]=0.0375621; 
	  bgdrate[7]=0.0197085; 
	  bgdrate[8]=0.0111295; 
	  bgdrate[9]=0.00927459; 
	  bgdrate[10]=0.00417356; 
	  bgdrate[11]=0.00521696; 
	  bgdrate[12]=0.000231865;      
	}
      else if (zenith==lowzd)//Aleksic et al 2016 
	{
	  // 0-30 deg values
	  crabrate[0]=0.818446;
	  crabrate[1]=3.01248;
	  crabrate[2]=4.29046;
	  crabrate[3]=3.3699;
	  crabrate[4]=1.36207;
	  crabrate[5]=1.21791;
	  crabrate[6]=0.880268;
	  crabrate[7]=0.579754;
	  crabrate[8]=0.299179;
	  crabrate[9]=0.166192;
	  crabrate[10]=0.0931911;
	  crabrate[11]=0.059986;
	  crabrate[12]=0.017854;
      
	  bgdrate[0]=3.66424;
	  bgdrate[1]=4.05919;
	  bgdrate[2]=2.41479;
	  bgdrate[3]=0.543629;
	  bgdrate[4]=0.0660764;
	  bgdrate[5]=0.0270313;
	  bgdrate[6]=0.0132653;
	  bgdrate[7]=0.00592351;
	  bgdrate[8]=0.00266975;
	  bgdrate[9]=0.00200231;
	  bgdrate[10]=0.00141831;
	  bgdrate[11]=0.00458864;
	  bgdrate[12]=0.0016686;
	}
    }

  for (int i=0; i<npoints; i++)
    {
      crabrate[i]*=offsetdegrad;
      bgdrate[i]*=offsetdegrad;
    }
}

// Li & Ma eq 17. significance, function abridged from MARS 
Double_t SignificanceLiMa(Double_t s, Double_t b, Double_t alpha)
{
    if( b==0 && s==0)
      return 0;

    b = b==0 ? FLT_MIN:b; // Guarantee that a value is output even in the case b or s == 0
    s = s==0 ? FLT_MIN:s; //   (which is not less allowed, possible, or meaningful than
                          //    doing it with b,s = 0.001)
    const Double_t sum = s+b;

    if (sum<0 || alpha<=0)
        return -1;
    const Double_t l = s*TMath::Log(s/sum*(alpha+1)/alpha);
    const Double_t m = b*TMath::Log(b/sum*(alpha+1)      );
    return l+m<0 ? -1 : TMath::Sqrt((l+m)*2);
}

Bool_t Checks()
{
  if (extension>1)
    {
      cout<<conError<<"Extension comparable to the size of the MAGIC camera cannot be simulated"<<conReset<<endl;
      return kFALSE;
    }
  if((numoff<=0)||(numoff>7))
    {
      cout<<conError<<"Number of OFF estimation regions must be in the range 1-7"<<conReset<<endl;
      return kFALSE;
    }
  if ((extension >0.5)&&(numoff>1))
    {
      cout<<conError<<"For large source extensions 1 OFF estimation region (numoff) should be used"<<conReset<<endl;
      return kFALSE;
    }

  if (isSUMT && (zenith!=lowzd))
    {
      cout<<conError<<"Only low zenith sensitivity is available for SUMT :-("<<conReset<<endl;
      return kFALSE;
    }

  if (isLSTmode && (isSUMT) )
    {
      cout<<conError<<"LST mode is not compatible with SUMT"<<conReset<<endl;
      return kFALSE;
    }

  if (offsetdegrad>1.00001)
    {
      cout<<conError<<"No cheating! the performance degradation ("<<offsetdegrad<<") should not be larger then 1"<<conReset<<endl;
      return kFALSE;
    }

  if (pulsarmode) 
    {
      if (pulsarOnRange<=0 ||pulsarOnRange>=1)
	{
	  cout<<conError<<"Pulsar mode ON phase range is "<<pulsarOnRange<<", and it should be in range (0,1)"<<conReset<<endl;
	  return kFALSE;
	}
      if (pulsarOffRange<=0 ||pulsarOffRange>=1)
	{
	  cout<<conError<<"Pulsar mode OFF phase range is "<<pulsarOffRange<<", and it should be in range (0,1)"<<conReset<<endl;
	  return kFALSE;
	}
	
      if (redshift >0)
	cout<<conWarn<<"Do you really want to observe a pulsar at redshift of "<<redshift<<" ??"<<conReset<<endl;
    }

  return kTRUE;
}


void mss()
{
  if (!Checks())
    {
      cout<<"exiting"<<endl;
      return;
    }

  if (redshift>0)
    LoadEBL(redshift);

  Prepare();
  // preparing reference SED graph;
  TF1 *ffluxcrab = new TF1 ("ffluxcrab", crabstr, eplotmin, eplotmax);
  TF1 *fsedcrab = new TF1 ("fsedcrab", Form("1.e-6*x*x*%s", crabstr.Data()), eplotmin, eplotmax);
  fsedcrab->SetLineColor(kGray+2);

  // TF1 *ffluxass = new TF1 ("ffluxass", assumedstr, eplotmin, eplotmax);
  ffluxint = new TF1 ("ffluxint", assumedstr, eplotmin, eplotmax);

  TF1 *ffluxass=0, *fsedass=0;
  if (redshift>0)
    {
      ffluxass = new TF1("ffluxass", fluxobs, eplotmin, eplotmax, 0);
      fsedass = new TF1("fsedass", sedobs, eplotmin, eplotmax, 0);
    }
  else
    {
      ffluxass=ffluxint;
      fsedass = new TF1 ("fsedass", Form("1.e-6*x*x*%s", assumedstr.Data()), eplotmin, eplotmax);
    }
  fsedass->SetLineColor(kGreen);

  TGraphErrors *gsed = new TGraphErrors();


  TCanvas *can = new TCanvas ("can", "", 640, 480);
  can->SetLogx(); can->SetLogy();
  can->SetGrid();
  TH1F *hframe = can->DrawFrame(eplotmin, yplotmin, eplotmax, yplotmax);
  hframe->GetXaxis()->SetTitle("E [GeV]");
  hframe->GetYaxis()->SetTitle("E^{2} dN/dE [TeV cm^{-2} s^{-1}]");
  hframe->GetYaxis()->SetTitleOffset(1.3);
  fsedcrab->Draw("same");
  fsedass->Draw("same");

  Double_t nexc_all=0, noff_all=0;
  Double_t best_int_e=-1, best_int_sigma=-1;

  Double_t pulsarOnOffRatio = pulsarOnRange/pulsarOffRange;
  for (int i=npoints-1; i>=0; i--)
    {
      #if ROOT_VERSION_CODE < ROOT_VERSION(6, 00, 00)
      Double_t dummy[3];
      Double_t intcrab = ffluxcrab->Integral(enbins[i], enbins[i + 1], dummy, 1.e-15);
      Double_t intass = ffluxass->Integral(enbins[i], enbins[i + 1], dummy, 1.e-15);
      #else
      Double_t intcrab = ffluxcrab->Integral(enbins[i], enbins[i + 1], 1.e-15);
      Double_t intass = ffluxass->Integral(enbins[i], enbins[i + 1], 1.e-15);
      #endif
      Double_t noff = bgdrate[i]*timeh*60; // number of off events
      noff*=(PSF*PSF+extension*extension)/(PSF*PSF); // larger integration cut due to extension

      Double_t dnoff=sqrt(noff/numoff); // error on the number of off events (computed from numoff regions)

      Double_t nexc = crabrate[i]*timeh*60 * intass/intcrab; 
      Double_t dexc = sqrt(nexc + noff + dnoff*dnoff); // error on the number of excess events 

      Double_t noffon=0;
      if (pulsarmode)
	{
	  // cout<<nexc<<" "<<noff<<" ";
	  noffon=noff *pulsarOnRange; // number of bgd events in ON phase
	  noff*=pulsarOffRange; // number of bgd events in OFF phase
	  dnoff=sqrt(noff)*pulsarOnOffRatio; //ignoring numoff for pulsars and scaling for the phase difference
	  dexc = sqrt(nexc + noffon + dnoff*dnoff); // error on the number of excess events 
	  // cout<<dnoff<<endl;
	}
      nexc_all+=nexc;
      noff_all+=noff;
      
      Double_t sigma = 0; 
      // for tiny excesses (1.e-7 events) the function below can have numerical problems, and either way sigma should be 0 then
      if (nexc>0.01)
	{ 
	  if (pulsarmode)
	    {
	      // cout<<nexc+noffon<<" "<<noff<<" "<<pulsarOnRange/pulsarOffRange<<endl;
	      sigma = SignificanceLiMa(nexc+noffon, noff, pulsarOnOffRatio);
	      noff=noffon; // needed later for SBR 
	    }
	  else
	    sigma = SignificanceLiMa(nexc+noff, noff*numoff, 1./numoff);
	}
      Bool_t detect=kFALSE; 
      if (pulsarmode)
	{
	  if ((sigma>=5.0)&&(nexc>minev))
	    detect=kTRUE;
	}
      else
	if ((sigma>=5.0)&&(nexc/noff>minSBR)&&(nexc>minev))
	  detect=kTRUE;
      cout<<Form("%.1f-%.1f GeV: exc. = %.1f+-%.1f ev., SBR=%.2f%%, sigma = %.1f", 
		 enbins[i], enbins[i+1], nexc, dexc, 100.*nexc/noff, sigma)<<
	(detect?" DETECTION":"")<<endl;

      
      Double_t sigma_all=-1;
      if (nexc_all>0.01)
	{
	  if (pulsarmode)
	    sigma_all = SignificanceLiMa(nexc_all+noff_all*pulsarOnOffRatio, noff_all, pulsarOnOffRatio);
	  else
	    sigma_all = SignificanceLiMa(nexc_all+noff_all, noff_all*numoff, 1./numoff);
	  cout<<Form("Integral significance > %.1f GeV = %.2f with %.1f excess",enbins[i],sigma_all,nexc_all)<<
	    (pulsarmode?"":Form(", SBR=%.2f",nexc_all/noff_all))<<endl;
	  if ((nexc_all>minev) && (pulsarmode || (nexc_all>minSBR*noff_all) ))
	    if (sigma_all>best_int_sigma)
	      {
		best_int_sigma=sigma_all;
		best_int_e = enbins[i];
	      }
	}
      
      if (nexc>minerror*dexc)
	{
	  Double_t en = sqrt(enbins[i]*enbins[i+1]);
	  Double_t sed = fsedass->Eval(en);
	  Double_t dsed = sed * dexc/nexc;
	  gsed->SetPoint(gsed->GetN(), en, sed);
	  gsed->SetPointError(gsed->GetN()-1, 0, dsed);
	  if (drawsigma)
	    {
	      TLatex *text = new TLatex (en, 2*sed, Form("%.1f#sigma", sigma));
	      text->SetTextColor(detect?kBlack:kGray);
	      text->SetTextSizePixels(14);
	      text->SetTextAngle(90);
	      text->Draw();
	    }
	}
    }
  TString inst=TString(isLSTmode?"MAGIC+LST1":"MAGIC");
  gsed->SetMarkerStyle(20);
  
  TLegend *leg = new TLegend (0.5, 0.75, 0.99, 0.99);
  TString zdtag="lowZd";
  if (zenith==midzd) zdtag="midZd";
  if (zenith==hizd) zdtag="highZd";
  
  TString header= inst + " " + zdtag + TString (isSUMT?", SUMT":"");
  leg->SetHeader(header);
  if (redshift>0)
    leg->AddEntry(fsedass, Form("%s (Assumed, z=%.3f)", srcname.Data(), redshift), "L");
  else
    leg->AddEntry(fsedass, Form("%s (Assumed)", srcname.Data()), "L");
  leg->AddEntry(gsed, Form("expected SED (T_{obs} = %.1f h)", timeh), "PE");
  leg->AddEntry(fsedcrab, "Crab (Aleksic et al 2016)", "L");

  gsed->Draw("PE");
  leg->Draw();
  if (drawsigma && (best_int_sigma>0))
    {
      TLatex *text = new TLatex();
      text->DrawLatexNDC(0.03, 0.93,  Form("#sigma (>%.1f GeV) = %.2f", best_int_e, best_int_sigma));
    }
  TLatex *lversion = new TLatex();
  lversion->SetTextSize(0.04);
  lversion->DrawLatexNDC(0.12, 0.12,
			 "mss v"+version
			 +(ismc?", (MC based)":"")
			 +(offsetdegrad<0.99?Form(", offset degr.=%.2f",offsetdegrad):"")
			 );
}
