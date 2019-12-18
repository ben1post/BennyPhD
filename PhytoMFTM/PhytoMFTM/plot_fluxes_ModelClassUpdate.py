from PhytoMFTM.ModelClasses_update import outarray, timedays
import numpy as np
import matplotlib.pyplot as plt
from pylab import cm
import matplotlib.colors as pltcol

#define global colors
#N
#ax2[i_plot].stackplot(timedays_ly, -PGains, labels=['Phyto Gains'], baseline='zero')
#ax2[i_plot].stackplot(timedays_ly, DRemin, NMixing, ZUnassimFeedNitrate,
#                      labels=['Remineralisation', 'Mixing', 'ZUnassimFeed'], baseline='zero')


cmap = cm.get_cmap('seismic', 100)    # PiYG
def HEXfromVal(val):
    rgb = cmap(val)[:3]
    return pltcol.rgb2hex(rgb)

cNMixing = HEXfromVal(0)

cPLinMort = HEXfromVal(5)
cPQuadMort = HEXfromVal(8)
cPMortality = HEXfromVal(6)
cPZooGrazed = HEXfromVal(15)
cPMixing = HEXfromVal(2)

cPGains = HEXfromVal(20)
cPNuts = HEXfromVal(18)
cPTemps = HEXfromVal(19)
cPLight = HEXfromVal(21)

cZGains = HEXfromVal(50)
cZLinMort = HEXfromVal(70)
cZQuadMort = HEXfromVal(75)
cZMixing = HEXfromVal(80)

cZUnassimFeedDetritus = HEXfromVal(52)

cDZooGrazed = HEXfromVal(55)
cDRemin = HEXfromVal(90)
cDMixing = HEXfromVal(100)


#P
#ax2[i_plot].stackplot(timedays_ly, -PLinMort, -PQuadMort, -PZooGrazed, -PMixing,
#                      labels=['Linear Mortality', 'Quad Mortality', 'Grazing', 'Mixing'], baseline='zero')
#ax2[i_plot].stackplot(timedays_ly, PNuts, PTemps, PLight,
#                      labels=['Nutrient Gains', 'Temp Dependency', 'Light Harvesting'], baseline='zero')
#Z
#ax2[i_plot].stackplot(timedays_ly, ZGains, labels=['Assimilated Grazing'], baseline='zero')
#ax2[i_plot].stackplot(timedays_ly, -ZLinMort, -ZQuadMort, -ZMixing, labels=['Mortality', 'Mixing'], baseline='zero')

#D
#ax2[i_plot].stackplot(timedays_ly, -DRemin, -DZooGrazed, -DMixing,
#                      labels=['D Remineralisation', 'D Zoograzing', 'Mixing'], baseline='zero')
#ax2[i_plot].stackplot(timedays_ly, ZUnassimFeedDetritus, ZLinMort, PMortality,
#                      labels=['Zoo UnassimFeeding', 'ZLinear Mort', 'Phyto Mortality'], baseline='zero')


def getOUT(out, index):
    out1 = [out[i + 1, index] - out[i, index] for i in range(len(out[:, 0]) - 1)]
    out2 = np.array(out1)
    out3 = np.concatenate([np.array(out[0, index]), out2], axis=None)
    return out3

def plotNfluxes(outarray, pfn, zn, i_plot, title):

    # PLOTTING
    timedays_ly = timedays[1:366]
    # truncate outarraySiNO to last year of 5:
    outarray_lyS = outarray[1460:1825]
    # color vectors
    colors = ['#808080','#d55e00', '#cc79a7', '#0072b2', '#009e73', '#009e73']
    alphas = [1., 0.8, 0.6, 0.4]
    lws = [1, 2.5, 4, 5.5]

    outindex = 4

    #Nitrate:
    N = outarray_lyS[:, 0]
    DRemin = getOUT(outarray_lyS, outindex+16)
    NMixing = getOUT(outarray_lyS, outindex+20)

    #Zoo:
    ZUnassimFeedNitrate = getOUT(outarray_lyS, outindex+14)

    #Phyto:
    PGains = getOUT(outarray_lyS, outindex+3)

    ax1[i_plot].set_title(title)

    ax1[i_plot].plot(timedays_ly, N, label='N', color='grey')
    ax1[i_plot].set_ylim(bottom=0)

    ax2[i_plot].stackplot(timedays_ly, -PGains, labels=['Phyto Gains'], baseline='zero', colors=[cPGains])
    ax2[i_plot].stackplot(timedays_ly, DRemin, NMixing,ZUnassimFeedNitrate, labels=['Remineralisation','Mixing','ZUnassimFeed'], baseline='zero', colors=[cDRemin,cNMixing,cZUnassimFeedDetritus])
    ax2[i_plot].plot(timedays_ly, DRemin+NMixing+ZUnassimFeedNitrate-PGains, label = 'Total Flux', color='black')
    ax2[i_plot].set_ylim(-0.3,0.3)


    ax1[i_plot].set_ylabel('Nitrate [µM N]')
    ax2[i_plot].set_ylabel('Nitrate Fluxes [µM N / day]')

    ax2[i_plot].legend(loc='lower right')
    ax1[i_plot].legend()

    ax2[i_plot].set_xlabel('Time [days]')
    ax2[i_plot].set_xlim(1,365)
    ax2[i_plot].set_xticks([1, 90, 180, 270, 365])

def plotSifluxes(outarray, pfn, zn, i_plot, title):

    # PLOTTING
    timedays_ly = timedays[1:366]
    # truncate outarraySiNO to last year of 5:
    outarray_lyS = outarray[1460:1825]
    # color vectors
    colors = ['#808080','#d55e00', '#cc79a7', '#0072b2', '#009e73', '#009e73']
    alphas = [1., 0.8, 0.6, 0.4]
    lws = [1, 2.5, 4, 5.5]

    # artist for legends
    FullArtist = plt.Line2D((0, 1), (0, 0), c=colors[4], alpha=alphas[1], lw=lws[0])

    outindex = 3+zn+pfn

    #Silicate:
    Si = outarray_lyS[:, 1]
    SiRemin = outarray_lyS[:, outindex + 2]
    SiMixing = outarray_lyS[:, outindex + 3]

    #Phyto:
    SiDrawdown = outarray_lyS[:, outindex + 6]

    ax1[i_plot].set_title(title)

    ax1[i_plot].plot(timedays_ly, Si, label='Si', color='grey')
    ax1[i_plot].set_ylim(bottom=0)

    ax2[i_plot].stackplot(timedays_ly, -SiDrawdown, labels=['Si Uptake'], baseline='zero')
    ax2[i_plot].stackplot(timedays_ly, SiRemin, SiMixing, labels=['Remineralisation','Mixing'], baseline='zero')
    ax2[i_plot].plot(timedays_ly, SiRemin+SiMixing-SiDrawdown, label = 'Total Flux', color='black')
    ax2[i_plot].set_ylim(-0.09,0.09)

    ax1[i_plot].set_ylabel('Silicate [µM Si]')
    ax2[i_plot].set_ylabel('Silicate Fluxes [µM Si / day]')

    ax2[i_plot].legend(loc='lower right')
    ax1[i_plot].legend()

    ax2[i_plot].set_xlabel('Time [days]')
    ax2[i_plot].set_xlim(1,365)
    ax2[i_plot].set_xticks([1, 90, 180, 270, 365])

def plotDetritusfluxes(outarray, pfn, zn, i_plot, title):

    # PLOTTING
    timedays_ly = timedays[1:366]
    # truncate outarraySiNO to last year of 5:
    outarray_lyS = outarray[1460:1825]
    # color vectors
    colors = ['#808080','#d55e00', '#cc79a7', '#0072b2', '#009e73', '#009e73']
    alphas = [1., 0.8, 0.6, 0.4]
    lws = [1, 2.5, 4, 5.5]

    outindex = 4

    #sum(ZooMortality) - NRemineralization - SiRemineralization - DetritusMixing + sum(UnassimilatedProduction) + sum([p[i].mortality(P[i]) for i in range(pfn)]) # Detritus

    #Detritus:
    Det = outarray_lyS[:, 3]
    ZUnassimFeedDetritus = getOUT(outarray_lyS, outindex+14) #+
    ZLinMort = getOUT(outarray_lyS, outindex+10) #+
    PMortality = getOUT(outarray_lyS, outindex+6) #+

    DRemin = getOUT(outarray_lyS, outindex+16) #-
    DZooGrazed = getOUT(outarray_lyS, outindex+17) #-
    DMixing = getOUT(outarray_lyS, outindex+18) #-


    ax1[i_plot].set_title(title)

    ax1[i_plot].plot(timedays_ly, Det, label='Det', color='grey')
    ax1[i_plot].set_ylim(bottom=0)

    ax2[i_plot].stackplot(timedays_ly, -DRemin, -DZooGrazed, -DMixing, labels=['D Remineralisation', 'D Zoograzing', 'Mixing'], baseline='zero')
    ax2[i_plot].stackplot(timedays_ly, ZUnassimFeedDetritus, ZLinMort, PMortality, labels=['Zoo UnassimFeeding', 'ZLinear Mort', 'Phyto Mortality'], baseline='zero')
    ax2[i_plot].plot(timedays_ly, -DRemin-DZooGrazed-DMixing+ZUnassimFeedDetritus+ZLinMort+PMortality, label = 'Total Flux', color='black')
    #ax2[i_plot].set_ylim(-0.09,0.09)


    ax1[i_plot].set_ylabel('Detritus [µM N]')
    ax2[i_plot].set_ylabel('Detritus Fluxes [µM N / day]')

    ax2[i_plot].legend(loc='lower right')
    ax1[i_plot].legend()

    ax2[i_plot].set_xlabel('Time [days]')
    ax2[i_plot].set_xlim(1,365)
    ax2[i_plot].set_xticks([1, 90, 180, 270, 365])

def plotPhyfluxes(outarray, pfn, zn, i_plot, title):

    # PLOTTING
    timedays_ly = timedays[1:366]
    # truncate outarraySiNO to last year of 5:
    outarray_lyS = outarray[1460:1825]
    # color vectors
    colors = ['#808080','#d55e00', '#cc79a7', '#0072b2', '#009e73', '#009e73']
    alphas = [1., 0.8, 0.6, 0.4]
    lws = [1, 2.5, 4, 5.5]


    outindex = 4

    # phy = [Gains[i] - PhytoGrazing[i] - PhytoMortality[i] - PhytoMixing[i] - PhytoSinking[i] for i in
    #       range(pfn)]  # Phytoplankton growth

    #Phytoplankton:
    P = outarray_lyS[:, 1]
    PGains = getOUT(outarray_lyS, outindex+3)  # +
    PLosses = getOUT(outarray_lyS, outindex+8)  # -

    PTempDepGrow = getOUT(outarray_lyS, outindex+0)  # +
    PNutUptake = getOUT(outarray_lyS, outindex+1)  # +
    PLightHarv = getOUT(outarray_lyS, outindex+2)  # +

    PLinMort = getOUT(outarray_lyS, outindex+4)  # -
    PQuadMort = getOUT(outarray_lyS, outindex+5)  # -
    PMortality = getOUT(outarray_lyS, outindex+6)  # -
    PZooGrazed = getOUT(outarray_lyS, outindex+7)  # -
    PMixing = getOUT(outarray_lyS, outindex+22)  # -

    PLight = PGains * PLightHarv / (PNutUptake + PLightHarv + PTempDepGrow)
    PNuts = PGains * PNutUptake / (PNutUptake + PLightHarv + PTempDepGrow)
    PTemps = PGains * PTempDepGrow / (PNutUptake + PLightHarv + PTempDepGrow)
    ax1[i_plot].set_title(title)

    ax1[i_plot].plot(timedays_ly, P, label='P', color='grey')
    ax1[i_plot].set_ylim(bottom=0)

    ax2[i_plot].stackplot(timedays_ly, -PLinMort, -PQuadMort, -PZooGrazed, -PMixing, labels=['Linear Mortality', 'Quad Mortality', 'Grazing', 'Mixing'], baseline='zero')
    ax2[i_plot].stackplot(timedays_ly, PNuts, PTemps, PLight, labels=['Nutrient Gains', 'Temp Dependency', 'Light Harvesting'], baseline='zero')
    ax2[i_plot].plot(timedays_ly, PGains-PLosses, label = 'Total Flux', color='black')
    #ax2[i_plot].set_ylim(-0.09,0.14)


    ax1[i_plot].set_ylabel('Phytoplankton [µM N]')
    ax2[i_plot].set_ylabel('Phytoplankton Fluxes [µM N / day]')

    ax2[i_plot].legend(loc='lower right')
    ax1[i_plot].legend()

    ax2[i_plot].set_xlabel('Time [days]')
    ax2[i_plot].set_xlim(1,365)
    ax2[i_plot].set_xticks([1, 90, 180, 270, 365])

def plotZoofluxes(outarray, pfn, zn, i_plot, title):

    # PLOTTING
    timedays_ly = timedays[1:366]
    # truncate outarraySiNO to last year of 5:
    outarray_lyS = outarray[1460:1825]
    # color vectors
    colors = ['#808080','#d55e00', '#cc79a7', '#0072b2', '#009e73', '#009e73']
    alphas = [1., 0.8, 0.6, 0.4]
    lws = [1, 2.5, 4, 5.5]

    outindex = 4

    # zoo = [ZooGrowth[i] - ZooMortality[i] - ZooMixing[i] for i in range(zn)]
    # Zooplankton losses due to mortality and mixing

    #Zooplankton:
    Z = outarray_lyS[:, 2]

    ZGains = getOUT(outarray_lyS, outindex+9)
    ZLinMort = getOUT(outarray_lyS, outindex+10)
    ZQuadMort = getOUT(outarray_lyS, outindex+11)
    ZMixing = getOUT(outarray_lyS, outindex+12)
    ZLosses = getOUT(outarray_lyS, outindex+13)

    ax1[i_plot].set_title(title)

    ax1[i_plot].plot(timedays_ly, Z, label='Z', color='grey')
    #ax1[i_plot].set_ylim(0, 0.7)

    ax2[i_plot].stackplot(timedays_ly, ZGains, labels=['Assimilated Grazing'], baseline='zero')
    ax2[i_plot].stackplot(timedays_ly, -ZLinMort, -ZQuadMort, -ZMixing, labels=['Linear Mortality','Quad Mortality', 'Mixing'], baseline='zero')
    ax2[i_plot].plot(timedays_ly, ZGains-ZLosses, label='Total Flux', color='black')
    #ax2[i_plot].set_ylim(-0.002,0.002)

    ax1[i_plot].set_ylabel('Zooplankton [µM N]')
    ax2[i_plot].set_ylabel('Zooplankton Fluxes [µM N / day]')

    ax2[i_plot].legend(loc='lower right')
    ax1[i_plot].legend()

    ax2[i_plot].set_xlabel('Time [days]')
    ax2[i_plot].set_xlim(1,365)
    ax2[i_plot].set_xticks([1, 90, 180, 270, 365])


# NITRATE     /// or Nitrogen????
f1, (ax1, ax2) = plt.subplots(2, 4, gridspec_kw = {'height_ratios':[1, 3]}, sharex='col', sharey='row')

plotNfluxes(outarray, 1, 1, 0, '1P1Z')
#plotSifluxes(out1P1Z, 1, 1, 1, '1P1Z')
plotPhyfluxes(outarray, 1, 1, 1, '1P1Z')
plotZoofluxes(outarray, 1, 1, 2, '1P1Z')
plotDetritusfluxes(outarray, 1, 1, 3, '1P1Z')

f1.align_ylabels()

plt.subplots_adjust(hspace=0)
#plt.tight_layout()

#plt.show()

#f1.savefig("fluxes01.png", dpi=1500)
