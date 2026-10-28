# ## MAGIC Source Simulator. 
# 
#    MAGIC Source Simulator. The macro takes a given spectral shape, 
#   and  using MAGIC performance in the range 40GeV-16TeV from 
#   Aleksic, J., et al., 2016, Astroparticle Physics, 72, 76 
#   (please cite this paper if you want to use this macro in a publication)
#   estimates what kind of signal can be observed by MAGIC.
#   It computes significances of each spectral point according to
#   Li & Ma 1983, ApJ, 272, 317, Eq. 17 definition.
# 
# ####   Terminal Output:
#   For each estimated energy range one gets the number of excess events,
#   signal-to-background ratio, significance, and information if a given
#   bin satisfies the conditions for the detection. 
# 
# ####   Plot Output:
#   Spectral points following the assumed spectrum are shown, their 
#   error bars reflect the performance of the instrument for such a source.
#   Significance for each bin is given (bins without detection
#   have gray numbers. In the top of the plot a simple number using 
#   information from all the shown bins is given to evaluate 
#   if the source can be detected (roughly if sum of sigma_i / sqrt(N) >~5) 
# 
# ####   requirements:
# - python version: Python 3, uses scipy
# 
# #### caveats:
# - the macro is operating on estimated energy doing simple
#   comparisons with differential in estimated energy rates seen for
#   Crab Nebula, for soft sources the differences in energy migration 
#   will result in different performance than calculated here
# - the threatment of extended sources is very approximate, only
#   increase in background is taken into account (without energy
#   dependence of PSF) for extended sources for extension >~0.4 deg the
#   dependence on the offset from the centre of the camera will further
#   worsen the performance w.r.t. the one calculated here
# - significances are given for each differential energy bin
#   separatelly, but to detect a source one normally applies a cut that
#   keeps a broad range of energies inside resulting in better integral
#   sensitivity than differential one. Also, optimization of cuts for 
#   a broad energy range usually results in somewhat better sensitivity
#   than what one can get by simply integrating the used here signal in 
#   differential energy bins. As a very crude approximation for detection
#   capability we calculate here also a sum of significances of all the
#   points in SED divided by the sqrt of the number of those points
# 
# #### usage:
# - modify the Assumed spectrum and "basic settings" below and run in terminal:
#       python3 mss.py 
# - for questions ask Julian Sitarek (jsitarek [at] uni.lodz.pl)
# - add ana.babic@fer.hr for python version
# 
# 
# #### version history of mss.cpp:
# - 1.9 (2024/10/25): Added high zenith performance for MAGIC+LST1 and MAGIC, fixed the combination of sensitivities
# - 1.8 (2024/01/29): Changed the calculation of combined sensitivity
# - 1.7 (2023/11/08): updated the LST1+MAGIC mid zenith sensitivity to the one in the paper
# - 1.5 (2022/10/24): again updated MAGIC+LST1 performance
# - 1.4 (2022/09/19): updated MAGIC+LST1 performance
# - 1.3 (2021/10/13): added also MAGIC+LST1 performance from MC study presented in ICRC 2019
# - 1.2 (2020/07/24): a number of additional features:
#       - EBL attenuation with Dominguez+11 model
#       - SUMT performance (only for low zenith)
#       - pulsar mode
# - 1.1 (2017/10/16): fixed a bug which made the errorbars too small
# mss.py was tested for compatibility with the ROOT (5.34/34) version mss.cpp


#get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from scipy.special import erfc,erfcinv
from scipy import interpolate
from enum import Enum
import matplotlib.lines as mlines

# Source settings
# [cm^-2 s^-1 TeV^-1] assumed (intrinsic) source spectrum (dN/dS dt dE) as a function of energy in GeV
# you can create your own, some common examples below  
def Assumed(x):
    #return 2.7e-9*pow(x/100., -3.18)                                              # simple power-law 
    return 2.0e-11*pow(x/1000., -1.99)*np.exp(-x/100)                               # power-law with an exponential cut-off
    #return 1.0e-11*pow(x/1000., -1.9)*np.exp(-pow(x/300, 0.5))                     # power-law with a non-exponential cut-off
    #return 0.007*3.39e-11*pow(x/1000., -2.51-0.21*np.log10(x/1000.))               # Crab Nebula [Aleksic et al. 2016]
    #return 0.1*np.exp(-x/200)*3.39e-11*pow(x/1000., -2.51-0.21*np.log10(x/1000.))  # 10% Crab with 200 GeV cut-off


# More examples of assumed spectrum:
# 
#     const Double_t redshift=0.94; 
#     const TString assumedstr="2.7e-9 * pow(x/100, -3.18)"; // simple power-law
#     const TString srcname="PKS1441"; // name of the source
# 
#     const TString assumedstr="2.28e-9 * pow(x/32.1, -5.62)"; // simple power-law
#     const TString srcname="Geminga"; // name of the source
# 
#     const TString assumedstr="4.e-9 * pow(x/100., -2.6)"; // simple power-law
#     const TString srcname="PG 1553"; // name of the source
# 
#     const TString assumedstr="1.4e-12 * pow(x/1310., -2.74)"; // simple power-law
#     const TString srcname="J1835-069"; // name of the source
# 
#     const TString assumedstr="7.9e-11 * pow(x/175, -3.26)"; // simple power-law
#     const Double_t redshift=0.360; // redshift of the source (for the EBL absorption), if -1 then no absorption
#     const TString srcname="PKS1510low"; // name of the source


# basic settings
timeh=20               # [h], time of observations
extension=0.0          # [deg] extension radius of the source

# source settings
redshift = -1   # redshift of the source (for the EBL absorption), if -1 then no absorption
srcname ="PWL+cutoff" # name of the source

Zeniths = Enum('Zeniths', {'lowzd':'lowZd', 'midzd':'midZd', 'hizd':'highZd'})

zenith = Zeniths.lowzd     # performance for 0-30deg zenith angles
# zenith = Zeniths.midzd   # performance for 30-45deg zenith angles
# zenith = Zeniths.hizd    # performance for ~60deg zenith angles

# you can check visibility of your source e.g. here: http://www.magic.iac.es/scheduler/
isSUMT = False

# advanced settings
numoff = 3        # number of background estimation regions
minev = 10.0      # minimum number of events
minSBR = 0.05     # minimum ratio of excess to background
PSF = 0.1         # [deg] PSF (for worsening the performance for extended sources)
offsetdegrad = 1.0 # degradation factor (applied to gamma and bgd. rates) if observations are taken at higher offset then 0.4 deg from camera center
# for best sensitivity region (~300 GeV at low zenith) and MAGIC-only or MAGIC+LST1 observations it can be estimated as 1.1*exp(-0.8 * offset^2)
# for proposals focussed on multi-TeV emission a slower drop can be assumed: ~exp(-0.3*offset^2)
eplotmin = 30     # [GeV] x plot range
eplotmax = 20.e3  # [GeV] x plot range
yplotmin = 1.e-14 # [TeV cm^-2 s^-1] y plot range
yplotmax = 1.e-9  # [TeV cm^-2 s^-1] y plot range
minerror = 2      # showing only points with value > minerror * error
drawsigma = True  # whether to draw also sigmas on the plot


# pulsar mode settings 
pulsarmode = False     # if true the background is reduced to on phase (see below) and SBR cut is ignored
pulsarOnRange = 0.092  # range of ON phases used for pulsar mode
pulsarOffRange = 0.25  # range of OFF phases used for pulsar mode

isLSTmode = False # uses MAGIC+LST1 (software trigger) performance (data in mid zenith, MC in low zenith)

# global variables (DO NOT MODIFY)
npoints = 13
crabrate = np.zeros(npoints) # [min^-1]
bgdrate  = np.zeros(npoints) # [min^-1]
enbins   = np.zeros(npoints+1) # [GeV]
pathebl = "dominguez_ebl_tau.txt"  # path with EBL model of Dominguez+11

version = "1.9"
ismc = False

#Crab Nebula [Aleksic et al. 2016]
def Crab(x):
    return 3.39e-11*pow(x/1000., -2.51-0.21*np.log10(x/1000.))

# Calculate SED value from energy and flux, for display
def CalcSED(x,flux):
    return 1.e-6*x*x*flux

def LoadEBL():
    f1 = open(pathebl, 'r')
    line = f1.readline()
    firstlineEBL = list(map(float, line.split()))
    nz = firstlineEBL[0]
    zz = firstlineEBL[1:]
    #print(nz, zz)
    f1.close()
    eblbody = np.loadtxt(pathebl,delimiter=' ', skiprows=1)
    energies = eblbody[:,0]*1.e3 #in GeV  #en*=1.e3; // GeV
    taus = eblbody[:,1:]
    if len(taus>0): 
        return zz, energies, taus

def FluxObs(z, xx, fluxint):
    from scipy.interpolate import RectBivariateSpline
    en = xx
    EBLz, EBLen, EBLtaus = LoadEBL()
    ftau = RectBivariateSpline(EBLz, EBLen, EBLtaus.T,kx=3, ky=3, s=0) # cubic, no smoothing
    tau =  ftau(z, en)[0] 
    atten = np.exp(-tau)
    return fluxint*atten

def Prepare():
    global crabrate, bgdrate, enbins, ismc
    npoints_array = np.arange(0, npoints+1)
    enbins = 100.*pow(10., (npoints_array-2)*0.2)  # [GeV]
    
    if isLSTmode:
        if (zenith == Zeniths.hizd):  # MC values values computed at 59 zd angle, might by ~20% too optimistic >~1 TeV
            crabrate[0] = 0  # 50GeV, below threshold
            crabrate[1] = 0. # 79 GeV, below threshol
            crabrate[2] = 0. # 126 GeV below threshold
            crabrate[3] = 0. # 200 GeV, below threshold
            crabrate[4] = 0.6929 # 316 GeV
            crabrate[5] = 2.0531 
            crabrate[6] = 1.5028 
            crabrate[7] = 0.8963
            crabrate[8] = 0.5915
            crabrate[9] = 0.3513
            crabrate[10]= 0.2314
            crabrate[11]= 0.111  # 7 TeV
            crabrate[12]= 0.0488  # 12.6 TeV
            
            bgdrate[0] = 0 #  below threshold
            bgdrate[1] = 0.
            bgdrate[2] = 0.
            bgdrate[3] = 0. # 200 GeV below threshold
            bgdrate[4] = 0.349
            bgdrate[5] = 0.2163
            bgdrate[6] = 0.0393
            bgdrate[7] = 0.0126
            bgdrate[8] = 0.0037
            bgdrate[9] = 0.0013
            bgdrate[10]= 0.0018
            bgdrate[11]= 0.0005
            bgdrate[12]= 0.0006 # 12.6 TeV            
            ismc = True 
            # the MC sensitivites are often too optimistic. While at low zenith MAGIC+LST1 sensitivity more or less matches
            # between the data and MCs, at mid zenith there is some ~20% difference >~1 TeV. Very likely the same happens
            # at high zenith, so this is why we artificially increase background by 40% to emulate this effect
            # and make more fair MAGIC-only vs MAGIC+LST1 comparisons at high zenith
            bgdrate[ np.sqrt(enbins[:-1]*enbins[1:])> 1000] *= 1.4
            
        elif (zenith == Zeniths.midzd): # values from https://doi.org/10.1051/0004-6361/202346927 
        	crabrate[0] = 0 # 50GeV, below threshold
        	crabrate[1] = 0.59 
        	crabrate[2] = 2.43 
        	crabrate[3] = 2.65 
        	crabrate[4] = 2.03 
        	crabrate[5] = 1.171
        	crabrate[6] = 0.899
        	crabrate[7] = 0.806
        	crabrate[8] = 0.319
        	crabrate[9] = 0.185
        	crabrate[10]= 0.113
        	crabrate[11]= 0.084
        	crabrate[12]=0  # missing data
        
        	bgdrate[0] = 0 # below threshold
        	bgdrate[1] = 0.715
        	bgdrate[2] = 1.42
        	bgdrate[3] = 0.307
        	bgdrate[4] = 0.093
        	bgdrate[5] = 0.0229
        	bgdrate[6] = 0.0093
        	bgdrate[7] = 0.0096
        	bgdrate[8] = 0.00264
        	bgdrate[9] = 0.0007
        	bgdrate[10]= 0.00138
        	bgdrate[11]= 0.00148
        	bgdrate[12]=0  # missing data
    
        elif (zenith == Zeniths.lowzd):
            # those numbers are taken from MC, but they agree
            # with a small bunch of data available 
            crabrate[0] = 0.3795 # 50GeV
            crabrate[1] = 2.5808 # 80GeV
            crabrate[2] = 2.8257
            crabrate[3] = 1.6654
            crabrate[4] = 1.3289
            crabrate[5] = 1.2871
            crabrate[6] = 0.8834
            crabrate[7] = 0.7045
            crabrate[8] = 0.3908
            crabrate[9] = 0.1963
            crabrate[10]= 0.089
            crabrate[11]= 0.0368
            crabrate[12]=0  # missing data
            
            bgdrate[0] = 8.8870e-01
            bgdrate[1] = 1.6505e+00
            bgdrate[2] = 5.7688e-01
            bgdrate[3] = 5.0804e-02
            bgdrate[4] = 2.0521e-02
            bgdrate[5] = 2.1718e-02
            bgdrate[6] = 5.6106e-03
            bgdrate[7] = 3.9491e-03
            bgdrate[8] = 3.2053e-03
            bgdrate[9] = 1.8074e-03
            bgdrate[10]= 7.3253e-04
            bgdrate[11]= 1.1352e-05
            bgdrate[12]=0  # missing data
            ismc = True

    else:
        if isSUMT:
            if (zenith == Zeniths.lowzd):
                # 2018/19 Crab data analyzed with generic ST0307 SUMT MCs
                # 1st point  ~50GeV, last point ~12TeV
                crabrate[0] =1.39684 # ~50GeV 
                crabrate[1] =3.12657
                crabrate[2] =3.09145 
                crabrate[3] =2.40268 
                crabrate[4] =1.32915 
                crabrate[5] =0.86180
                crabrate[6] =0.51666 
                crabrate[7] =0.31533 
                crabrate[8] =0.16207 
                crabrate[9] =0.09279 
                crabrate[10]=0.04624 
                crabrate[11]=0.02345 
                crabrate[12]=0.00874 # ~12TeV
                
                bgdrate[0] =3.33321
                bgdrate[1] =3.24046 
                bgdrate[2] =1.32361 
                bgdrate[3] =0.406588 
                bgdrate[4] =0.091944 
                bgdrate[5] =0.032226 
                bgdrate[6] =0.007277 
                bgdrate[7] =0.003123 
                bgdrate[8] =0.001487 
                bgdrate[9] =0.001464 
                bgdrate[10]=0.001231 
                bgdrate[11]=0.001152 
                bgdrate[12]=0.000957
            else:
                print('Only low zenith sensitivity is available for SUMT :-(')
                # to be implemented
        else:
            if (zenith == Zeniths.hizd): # Crab results from 2.5 hrs of 55-62 zd data from 2016-2018. Dataset from Juliane van Scherpenberg.
                crabrate[0] = 0 # 50GeV, below threshold
                crabrate[1] = 0. # 79 GeV, below threshol
                crabrate[2] = 0. # 126 GeV below threshold
                crabrate[3] = 0. # 200 GeV, below threshold
                crabrate[4] = 0.503462 # 316 GeV
                crabrate[5] = 1.60232 
                crabrate[6] = 2.26558 
                crabrate[7] = 0.928094
                crabrate[8] = 0.698335
                crabrate[9] = 0.305662
                crabrate[10]= 0.173859
                crabrate[11]= 0.083892 # 7 TeV
                crabrate[12]= 0.069938 # 12.6 TeV
            
                bgdrate[0] = 0 # below threshold
                bgdrate[1] = 0.
                bgdrate[2] = 0.
                bgdrate[3] = 0. # 200 GeV below threshold
                bgdrate[4] = 0.750815
                bgdrate[5] = 0.597588
                bgdrate[6] = 0.564753
                bgdrate[7] = 0.089775
                bgdrate[8] = 0.0568584
                bgdrate[9] = 0.00954936
                bgdrate[10]= 0.00344762
                bgdrate[11]= 0.00147755
                bgdrate[12]= 0.00229841 # 12.6 TeV
                ismc = False
            elif (zenith == Zeniths.midzd): # Aleksic et al 2016 
                crabrate[0]=0 # below threshold !!!
                crabrate[1]=0.404836
                crabrate[2]=3.17608 
                crabrate[3]=2.67108 
                crabrate[4]=2.86307 
                crabrate[5]=1.76124
                crabrate[6]=1.43988 
                crabrate[7]=0.944385 
                crabrate[8]=0.673335 
                crabrate[9]=0.316263 
                crabrate[10]=0.200331 
                crabrate[11]=0.0991222 
                crabrate[12]=0.0289831
                
                bgdrate[0]=1.67777 # below threashold
                bgdrate[1]=2.91732 
                bgdrate[2]=2.89228 
                bgdrate[3]=0.542563 
                bgdrate[4]=0.30467 
                bgdrate[5]=0.0876449
                bgdrate[6]=0.0375621 
                bgdrate[7]=0.0197085 
                bgdrate[8]=0.0111295 
                bgdrate[9]=0.00927459 
                bgdrate[10]=0.00417356 
                bgdrate[11]=0.00521696 
                bgdrate[12]=0.000231865      
            elif (zenith == Zeniths.lowzd): # Aleksic et al 2016 
                # 0-30 deg values, 
                crabrate[0]=0.818446
                crabrate[1]=3.01248
                crabrate[2]=4.29046
                crabrate[3]=3.3699
                crabrate[4]=1.36207
                crabrate[5]=1.21791
                crabrate[6]=0.880268
                crabrate[7]=0.579754
                crabrate[8]=0.299179
                crabrate[9]=0.166192
                crabrate[10]=0.0931911
                crabrate[11]=0.059986
                crabrate[12]=0.017854
                
                bgdrate[0]=3.66424
                bgdrate[1]=4.05919
                bgdrate[2]=2.41479
                bgdrate[3]=0.543629
                bgdrate[4]=0.0660764
                bgdrate[5]=0.0270313
                bgdrate[6]=0.0132653
                bgdrate[7]=0.00592351
                bgdrate[8]=0.00266975
                bgdrate[9]=0.00200231
                bgdrate[10]=0.00141831
                bgdrate[11]=0.00458864
                bgdrate[12]=0.0016686
        # offset degradation is only applied if NOT in LST mode 
        crabrate *=offsetdegrad; 
        bgdrate *=offsetdegrad;
    
    return npoints_array, enbins, crabrate, bgdrate



#// Li & Ma eq 17. significance, function abridged from MARS 
def SignificanceLiMa(s, b, alpha):
    if b==0 and s==0:
        return 0

    b = FLT_MIN if b==0 else b #// Guarantee that a value is output even in the case b or s == 0
    s = FLT_MIN if s==0 else s #//   (which is not less allowed, possible, or meaningful than
                               #//    doing it with b,s = 0.001)
    sumsb = s+b

    if sumsb<0 or alpha<=0:
        return -1;
    l = s*np.log(s/sumsb*(alpha+1)/alpha)
    m = b*np.log(b/sumsb*(alpha+1)      )

    return -1 if l+m<0 else  np.sqrt((l+m)*2)


def Checks():
    if extension>1:
        print("Extension comparable to the size of the MAGIC camera cannot be simulated")
        return False
    if (numoff<=0) or (numoff>7):
        print("Number of OFF estimation regions must be in the range 1-7")
        return False
    if (extension >0.5) and (numoff>1):
        print("For large source extensions 1 OFF estimation region (numoff) should be used")
        return False
    if (isSUMT and  (zenith != Zeniths.lowzd)):
        print("Only low zenith sensitivity is available for SUMT :-(")
        return False
    if (isLSTmode and isSUMT):
        print("LST mode is not compatible with SUMT")        
        return False
    if (offsetdegrad>1.00001):
        print("No cheating! the performance degradation ({0}) should not be larger then 1".format(offsetdegrad))
        return False
    if pulsarmode:
        if (pulsarOnRange<=0 or pulsarOnRange>=1):
            print("Pulsar mode ON phase range is {0}, and it should be in range (0,1)".format(pulsarOnRange))
            return False
        if (pulsarOffRange<=0 or pulsarOffRange>=1):
            print("Pulsar mode OFF phase range is {0}, and it should be in range (0,1)".format(pulsarOffRange))
            return False
        if (redshift > 0):
            print("Do you really want to observe a pulsar at redshift of {0} ??".format(redshift))
    return True


## MAIN code mss.cpp
if not Checks():
     print("exiting")
    #return 0
else:
    npoints_array, enbins, crabrate, bgdrate = Prepare()

    
    en = []; sed = []; dsed = []; sigmas = []; detected = []; 
    nexc_all = 0; noff_all = 0;
    best_int_e = -1; best_int_sigma = -1;
    pulsarOnOffRatio = pulsarOnRange/pulsarOffRange

    for i, e1, e2 in reversed(list(zip(range(0,len(enbins)),enbins, enbins[1:]))):
        intcrab, error =  quad(Crab,    e1, e2)
        if redshift>0:
            intass,  error =  quad(lambda x: FluxObs(redshift,x,Assumed(x)),  e1, e2)
        else:
            intass,  error =  quad(Assumed, e1, e2)
        noff = bgdrate[i]*timeh*60 # number of off events
        noff *= (PSF*PSF+extension*extension)/(PSF*PSF) #larger integration cut due to extension
        dnoff = np.sqrt(noff/numoff) # error on the number of off events (computed from numoff regions)
        nexc = crabrate[i]*timeh*60*intass/intcrab
        dexc = np.sqrt(nexc + noff + dnoff*dnoff)  # error on the number of excess events 

        noffon=0
        if pulsarmode:
            noffon = noff*pulsarOnRange # number of bgd events in ON phase
            noff *= pulsarOffRange      # number of bgd events in OFF phase
            dnoff = np.sqrt(noff)*pulsarOnOffRatio #ignoring numoff for pulsars and scaling for the phase difference
            dexc = np.sqrt(nexc + noffon + dnoff*dnoff) # error on the number of excess events 
        nexc_all+=nexc
        noff_all+=noff

        # for tiny excesses (1.e-7 events) the function below can have numerical problems, and either way sigma should be 0 then
        sigma=0 
        if nexc>0.01:
            if (pulsarmode):
                sigma = SignificanceLiMa(nexc+noffon, noff, pulsarOnOffRatio)
                noff = noffon   # needed later for SBR 
            else:
                sigma = SignificanceLiMa(nexc+noff, noff*numoff, 1./numoff)

        detect = False
        if pulsarmode:
            if ((sigma>=5.0) and (nexc>minev)):
                detect = True
        else:
            if sigma>=5.0 and nexc/noff>minSBR and nexc>minev:
                detect = True
        print('{0:.1f}-{1:.1f} GeV: exc. = {2:.1f}+-{3:.1f} ev., SBR={4:.2f}%, sigma = {5:.1f}'.format(enbins[i], enbins[i+1], nexc, dexc, 100.*nexc/noff if noff>0 else np.nan, sigma)," DETECTION" if detect else " ")

        sigma_all = -1
        if (nexc_all>0.01):
            if (pulsarmode):
                sigma_all = SignificanceLiMa(nexc_all+noff_all*pulsarOnOffRatio, noff_all, pulsarOnOffRatio);
            else:
                sigma_all = SignificanceLiMa(nexc_all+noff_all, noff_all*numoff, 1./numoff);
                print('Integral significance > {0:.1f} GeV = {1:.2f} with {2:.1f} excess'.format(enbins[i],sigma_all,nexc_all),(' ' if pulsarmode else ', SBR={0:.2f}'.format(nexc_all/noff_all)))
            if ((nexc_all>minev) and (pulsarmode or (nexc_all>minSBR*noff_all) )):
                if (sigma_all>best_int_sigma):
                    best_int_sigma=sigma_all
                    best_int_e = enbins[i]

        #print(nexc, minerror, dexc)
        if nexc>minerror*dexc:
            tmpen = np.sqrt(enbins[i]*enbins[i+1])
            en.append(tmpen)
            if redshift>0:
                tmpsed = CalcSED(tmpen, FluxObs(redshift,tmpen,Assumed(tmpen))[0])
            else:
                tmpsed = CalcSED(tmpen, Assumed(tmpen))
            sed.append(float(tmpsed))
            dsed.append(float(tmpsed*dexc/nexc))
            sigmas.append(sigma)
            detected.append(detect)

    inst = "MAGIC+LST1" if isLSTmode else "MAGIC"

    # preparing reference SED graph
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set(xlim=(eplotmin, eplotmax), ylim=(yplotmin, yplotmax),
        xlabel='E [GeV]',
        ylabel='E$^2$ dN/dE [TeV cm$^{-2}$ s$^{-1}$]')
    x = np.logspace(np.log10(eplotmin), np.log10(eplotmax), 50)
    labeltext = "expected SED ($T_{obs}$ = "+str(timeh)+" h)"
    plt.plot(x, 1.e-6*x*x*Crab(x), '0.3', label="Crab (Aleksic et al 2016)" )
    if redshift>0:
        plt.plot(x,1.e-6*x*x*FluxObs(redshift,x,Assumed(x)), 'limegreen', label = srcname+" (Assumed, z={0:.2f})".format(redshift) )
    else:
        plt.plot(x, 1.e-6*x*x*Assumed(x), 'limegreen', label = srcname+" (Assumed)" )
    if len(en)>0: 
        ax.errorbar(en, sed, yerr = dsed, label = labeltext, color='0', fmt='o')
        handles, labels = ax.get_legend_handles_labels()
    else:
        handles, labels = ax.get_legend_handles_labels()
        scatter_proxy = mlines.Line2D([], [] , color='0', marker='o', linestyle='None', markersize=3)
        handles.append(scatter_proxy)
        labels.append(labeltext)

    #zdtag = "lowZd";
    #if (zenith == Zeniths.midzd): zdtag="midZd"
    #if (zenith == Zeniths.hizd): zdtag="highZd"
    zdtag =  zenith.value
    
    header=""
    header = inst + " " + zdtag + (", SUMT " if isSUMT else "")


    ax.legend(handles=handles, labels=labels, loc='upper right', title = header)
    ax.grid(True,which="both",axis = 'x', ls="--",color='0.95')
    ax.grid(True,which="major",axis = 'y', ls="--",color='0.95')
    if (drawsigma):
        for i in range(0,len(sigmas)):
            col = '0' if detected[i] else '0.75'
            ax.text(en[i], 2*sed[i], "{0:.1f}$\sigma$".format(sigmas[i]),     rotation=90, size=10, ha='left', va='bottom', color=col )
    ax.annotate('mss v'+version + 
                (', (MC based)' if ismc else '')+
                (f', offset degr.={offsetdegrad}' if offsetdegrad<0.99 else ''),
            xy=(0.02, 0.02), xycoords='axes fraction')
    if (drawsigma and (best_int_sigma>0)):
        ax.set(title = "$\sigma$ (>{0:.1f} GeV) = {1:.2f}".format( best_int_e, best_int_sigma))
    fig.savefig('SED'+srcname+'.pdf');
    plt.show()
    plt.close(fig)

