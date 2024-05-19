# Plot absorbance spectra from USB2000 spectrometer
# 19-May-2024 J.Beale

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt 
import bottleneck as bn

aMin = 0.001  # minimum valid absorbance
aMax = 3      # max valid absorbance
wSize= 11  # rolling average window size
pkStart = 370  # minimum wavelength peak of interest
pkStop = 800   # maximum wavelength peak of interest

# ------------------------------------------------------------

# replace all NaN values with 0
def rmN(A):
    A[np.isnan(A)] = 0
    return(A)

# boxcar-filter array A with window size w
# w should be an odd number, otherwise result has offset
def filt(A, w):
    b = bn.move_mean(A, window=w) # should left-shift result by (n-1)/2
    b = rmN(b)
    res = np.roll(b, -int((w-1)/2) )
    return(res)

# return array of repeat-counts of same values in input array A
def runlength(A):
    ret = np.diff(np.where( np.concatenate(([A[0]], A[:-1] != A[1:],[True])) ) [0])[::2]
    return (ret)

# find location of local maxima in A, using filter window size win
def findPeaks(A, win):
    Af =np.clip(np.concatenate([np.array([0]), np.diff(A)]),-0.02,0.02)
    Af = filt(Af, win)
    Af = filt(Af, win)

    s1 = np.sign(Af)
    sign1 = np.concatenate([(np.diff( s1 ) != 0)*1,np.array([0])])
    ppeak1 = (sign1==s1)
    return (ppeak1)

# Find and display peak location in units of 'pos' within range (start,stop)
def getPeaks(data, pos, wSize, start, stop):
    peaks = findPeaks(data, wSize)                # find peaks after highpass filtering
    idxStart = pos.searchsorted(start, 'right') - 1  # last index befor nm range of interest
    idxStop = pos.searchsorted(stop, 'right')  # first index after nm range of interest
    peaks[0:idxStart] = 0
    peaks[idxStop:] = 0
    pkIdx = np.argwhere(peaks)[:,0]
    #print("Peak count: %d" % (pkIdx.size))
    dCount = 0
    for i in pkIdx:  # display position of peaks in nm
        pk = pos[i]
        amp = data[i]
        if (pk >= start) and (pk <= stop):
            print("pk = %5.3f  size: %5.3f" % (pk, amp))
            dCount += 1
    print("Peaks: %d" % dCount)
    return(peaks)

# ------------------------------------------------------------
refFile1 = r"C:\Users\beale\oceanview\Exp-15-May-2024\Subt4__7__21-53-27-718.txt"
refFile2 = r"C:\Users\beale\oceanview\Exp-15-May-2024\Subt4__12__21-54-51-920.txt"
sFile1 = r"C:\Users\beale\oceanview\Exp-15-May-2024\Subt4__8__21-53-42-501.txt"

sFile2 = r"C:\Users\beale\oceanview\Exp-15-May-2024\Subt4__9__21-53-58-472.txt"
sFile3 = r"C:\Users\beale\oceanview\Exp-15-May-2024\Subt4__10__21-54-12-184.txt"
sFile4 = r"C:\Users\beale\oceanview\Exp-15-May-2024\Subt4__11__21-54-30-064.txt"
sFile4 = r"C:\Users\beale\oceanview\Exp-15-May-2024\Subt4__11__21-54-30-064.txt"
sFile5 = r"C:\Users\beale\oceanview\Absorbance_TreeMoss__1__22-31-11-978.txt"


# refFile1 = r"C:\Users\beale\oceanview\Moss\USB2E15011__0__22-34-17-306.txt"
# refFile2 = r"C:\Users\beale\oceanview\Moss\USB2E15011__2__22-34-46-049.txt"
# sFile1 = r"C:\Users\beale\oceanview\Moss\USB2E15011__1__22-34-32-111.txt"


df1 = pd.read_csv(refFile1,sep='\t',skiprows=(14),header=None)
df2 = pd.read_csv(refFile2,sep='\t',skiprows=(14),header=None)
dfs1 = pd.read_csv(sFile1,sep='\t',skiprows=(14),header=None)
dfs2 = pd.read_csv(sFile2,sep='\t',skiprows=(14),header=None)
dfs3 = pd.read_csv(sFile3,sep='\t',skiprows=(14),header=None)
dfs4 = pd.read_csv(sFile4,sep='\t',skiprows=(14),header=None)
dfs5 = pd.read_csv(sFile5,sep='\t',skiprows=(14),header=None)


nm = df1[0].to_numpy() # wavelengths in nm
r1 = df1[1].to_numpy() # reference 1
r2 = df2[1].to_numpy()

s1 = dfs1[1].to_numpy()  # sample #1
s2 = dfs2[1].to_numpy()  # sample #2
s3 = dfs3[1].to_numpy()  # sample #3
s4 = dfs4[1].to_numpy()  # sample #4
s5 = dfs5[1].to_numpy()  # sample #4

#print("Samples: %d" % np.size(s4))

rAvg = (r1 + r2)/2
abs1 = -np.log10(np.clip(s1/rAvg,aMin,aMax))  # absorbance units
abs2 = -np.log10(np.clip(s2/rAvg,aMin,aMax))
abs3 = -np.log10(np.clip(s3/rAvg,aMin,aMax))
abs4 = -np.log10(np.clip(s4/rAvg,aMin,aMax))
abs5 = -np.log10(np.clip(s5/rAvg,aMin,aMax))

abs1f = filt(abs1, wSize)
abs2f = filt(abs2, wSize)
abs3f = filt(abs3, wSize)
abs4f = filt(abs4, wSize)
abs5f = filt(abs5, wSize)

nmf = bn.move_mean(nm, window=wSize)  # nanometers, filtered
# ======================================================================

specF = abs3f  # find peaks and display this absorption spectrum

# ======================================================================
pkSearchWin = 41 # array element range, must be odd
specFF = filt(specF, pkSearchWin)
for i in range(2):
    specFF = filt(specFF, pkSearchWin)

specHP = specF - specFF
specHPF = filt(specHP, pkSearchWin*2)
specHPF = filt(specHPF, pkSearchWin*2)
specHPP = specHP - specHPF  # 2nd stage highpass filter
specHPP = np.clip(specHPP,0.002,10)  # clip all negative values to zero

specHP = specHPP
specHP = specHP * 5

ppeak2 = getPeaks(specHP, nm, wSize, pkStart, pkStop)   # get peaks after highpass filter

plt.plot(nm, specF)
# plt.plot(nm,ppeak2/2)  # show position of peaks 

pkI = ppeak2.nonzero()[0]
pk = specF[pkI] + 0.05    # offset in graph Y units, to place above spectrum plot line, not on it
nmPk = nm[pkI]
plt.plot(nmPk, pk, "kv")
for i in pkI:
    x = nm[i]
    y = specF[i] + 0.1
    s = ("%5.1f" % x)
    if (x > 420):
        plt.text(x,y,s,horizontalalignment='center')

plt.title("Sample 3: Copper Beech")
plt.xlabel('wavelength, nm')
plt.ylabel('absorption, OD')
plt.grid(axis='both')
plt.xlim(350, 800)
plt.show()
