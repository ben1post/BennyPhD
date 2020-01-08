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


cmap = cm.get_cmap('cubehelix', 100)    # PiYG
def HEXfromVal(val):
    rgb = cmap(val)[:3]
    return pltcol.rgb2hex(rgb)

cNMixing = HEXfromVal(20)

cPLinMort = HEXfromVal(70)
cPQuadMort = HEXfromVal(75)
cPMortality = HEXfromVal(65)
cPZooGrazed = HEXfromVal(51)
cPMixing = HEXfromVal(20)

cPGains = HEXfromVal(25)
cPNuts = HEXfromVal(23)
cPTemps = HEXfromVal(22)
cPLight = HEXfromVal(34)

cZGains = HEXfromVal(50)
cZLinMort = HEXfromVal(70)
cZQuadMort = HEXfromVal(75)
cZMixing = HEXfromVal(20)

cZUnassimFeedDetritus = HEXfromVal(45)

cDZooGrazed = HEXfromVal(55)
cDRemin = HEXfromVal(80)
cDMixing = HEXfromVal(20)

cTOTALFLUX = HEXfromVal(3)


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
    out1x = out1[1:len(out1)]
    out2 = np.array(out1)
    #out3 = np.concatenate([np.array(out[0, index]), out2], axis=None)
    return out2

def plotNfluxes(outarray, pfn, zn, i_plot, title):

    # PLOTTING
    timedays_ly = timedays[1:366]
    # truncate outarraySiNO to last year of 5:
    outarray_ly = outarray[1460:1825]
    outarray_lyS = outarray[1459:1825]
    # color vectors
    colors = ['#808080','#d55e00', '#cc79a7', '#0072b2', '#009e73', '#009e73']
    alphas = [1., 0.8, 0.6, 0.4]
    lws = [1, 2.5, 4, 5.5]

    outindex = 4

    #Nitrate:
    N = outarray_ly[:, 0]
    DRemin = getOUT(outarray_lyS, outindex+16)
    NMixing = getOUT(outarray_lyS, outindex+20)

    #Zoo:
    ZUnassimFeedNitrate = getOUT(outarray_lyS, outindex+14)

    #Phyto:
    PGains = getOUT(outarray_lyS, outindex+3)

    ax1[i_plot].set_title(title)

    ax1[i_plot].plot(timedays_ly, N, label='N', color='grey')
    ax1[i_plot].set_ylim(bottom=0)

    ax2[i_plot].stackplot(timedays_ly, -PGains, labels=['P Growth'], baseline='zero', colors=[cPGains])
    ax2[i_plot].stackplot(timedays_ly, DRemin, NMixing,ZUnassimFeedNitrate, labels=['D Remineralisation','N Mixing','Z Unassimilated Feeding'], baseline='zero', colors=[cDRemin,cNMixing,cZUnassimFeedDetritus])
    ax2[i_plot].plot(timedays_ly, DRemin+NMixing+ZUnassimFeedNitrate-PGains, label = 'total growth/loss-rate', color=cTOTALFLUX)
    ax2[i_plot].set_ylim(-0.3,0.3)


    ax1[i_plot].set_ylabel('Nitrate [µM N]')
    ax2[i_plot].set_ylabel('Nitrate growth/loss rates  [µM N / day]')

    ax2[i_plot].legend(loc='lower right', fontsize='x-small')
    ax1[i_plot].legend()

    ax2[i_plot].set_xlabel('Time [days]')
    ax2[i_plot].set_xlim(1,365)
    ax2[i_plot].set_xticks([1, 90, 180, 270, 365])

def plotDetritusfluxes(outarray, pfn, zn, i_plot, title):

    # PLOTTING
    timedays_ly = timedays[1:366]
    # truncate outarraySiNO to last year of 5:
    outarray_ly = outarray[1460:1825]
    outarray_lyS = outarray[1459:1825]
    # color vectors
    colors = ['#808080','#d55e00', '#cc79a7', '#0072b2', '#009e73', '#009e73']
    alphas = [1., 0.8, 0.6, 0.4]
    lws = [1, 2.5, 4, 5.5]

    outindex = 4

    #sum(ZooMortality) - NRemineralization - SiRemineralization - DetritusMixing + sum(UnassimilatedProduction) + sum([p[i].mortality(P[i]) for i in range(pfn)]) # Detritus

    #Detritus:
    Det = outarray_ly[:, 3]
    ZUnassimFeedDetritus = getOUT(outarray_lyS, outindex+14) #+
    ZLinMort = getOUT(outarray_lyS, outindex+10) #+
    PMortality = getOUT(outarray_lyS, outindex+6) #+

    DRemin = getOUT(outarray_lyS, outindex+16) #-
    DZooGrazed = getOUT(outarray_lyS, outindex+17) #-
    DMixing = getOUT(outarray_lyS, outindex+18) #-


    ax1[i_plot].set_title(title)

    ax1[i_plot].plot(timedays_ly, Det, label='Det', color='grey')
    ax1[i_plot].set_ylim(bottom=0)

    ax2[i_plot].stackplot(timedays_ly, -DRemin, -DZooGrazed, -DMixing, labels=['D Remineralisation', 'Z Grazing', 'D Mixing/Sinking'], baseline='zero', colors=[cDRemin,cDZooGrazed,cDMixing])
    ax2[i_plot].stackplot(timedays_ly, ZUnassimFeedDetritus, ZLinMort, PMortality, labels=['Z Unassimilated Feeding', 'Z Linear Mortality', 'P Mortality'], baseline='zero', colors=[cZUnassimFeedDetritus,cZLinMort,cPMortality])
    ax2[i_plot].plot(timedays_ly, -DRemin-DZooGrazed-DMixing+ZUnassimFeedDetritus+ZLinMort+PMortality, label = 'total growth/loss-rate', color=cTOTALFLUX)
    #ax2[i_plot].set_ylim(-0.09,0.09)


    ax1[i_plot].set_ylabel('Detritus [µM N]')
    ax2[i_plot].set_ylabel('Detritus growth/loss rates  [µM N / day]')

    ax2[i_plot].legend(loc='lower right', fontsize='x-small')
    ax1[i_plot].legend()

    ax2[i_plot].set_xlabel('Time [days]')
    ax2[i_plot].set_xlim(1,365)
    ax2[i_plot].set_xticks([1, 90, 180, 270, 365])
    ax2[i_plot].set_ylim(-0.3,0.3)

def plotPhyfluxes(outarray, pfn, zn, i_plot, title):

    # PLOTTING
    timedays_ly = timedays[1:366]
    # truncate outarraySiNO to last year of 5:
    outarray_ly = outarray[1460:1825]
    outarray_lyS = outarray[1459:1825]
    # color vectors
    colors = ['#808080','#d55e00', '#cc79a7', '#0072b2', '#009e73', '#009e73']
    alphas = [1., 0.8, 0.6, 0.4]
    lws = [1, 2.5, 4, 5.5]


    outindex = 4

    # phy = [Gains[i] - PhytoGrazing[i] - PhytoMortality[i] - PhytoMixing[i] - PhytoSinking[i] for i in
    #       range(pfn)]  # Phytoplankton growth

    #Phytoplankton:
    P = outarray_ly[:, 1]
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

    ax2[i_plot].stackplot(timedays_ly, -PLinMort, -PQuadMort, -PZooGrazed, -PMixing, labels=['P Linear Mortality', 'P Quadratic Mortality', 'Z Grazing', 'P Mixing'], baseline='zero', colors=[cPLinMort,cPQuadMort,cPZooGrazed,cPMixing])
    ax2[i_plot].stackplot(timedays_ly, PGains, labels=['P Growth'], baseline='zero', colors=[cPGains])
    ax2[i_plot].plot(timedays_ly, PGains-PLosses, label='total growth/loss-rate', color=cTOTALFLUX)
    #ax2[i_plot].set_ylim(-0.09,0.14)


    ax1[i_plot].set_ylabel('Phytoplankton biomass [µM N]')
    ax2[i_plot].set_ylabel('Phytoplankton growth/loss rates [µM N / day]')

    ax2[i_plot].legend(loc='lower right', fontsize='x-small')
    ax1[i_plot].legend()

    ax2[i_plot].set_xlabel('Time [days]')
    ax2[i_plot].set_xlim(1,365)
    ax2[i_plot].set_xticks([1, 90, 180, 270, 365])
    ax2[i_plot].set_ylim(-0.3,0.3)

def plotPhyLim(outarray, pfn, zn, i_plot, title):
    # PLOTTING
    timedays_ly = timedays[1:366]
    # truncate outarraySiNO to last year of 5:
    outarray_lyS = outarray[1459:1825]
    # color vectors

    outindex = 4

    PTempDepGrow = getOUT(outarray_lyS, outindex+0)  # +
    PNutUptake = getOUT(outarray_lyS, outindex+1)  # +
    PLightHarv = getOUT(outarray_lyS, outindex+2)  # +

    print(PTempDepGrow)

    ax1[i_plot].set_title(title)

    ax1[i_plot].plot(timedays_ly, PTempDepGrow, label='VpT', color='grey')
    ax1[i_plot].set_ylim(bottom=0)

    ax2[i_plot].plot(timedays_ly, PLightHarv, label='Light Lim', color=cPLight)

    ax2[i_plot].plot(timedays_ly, PNutUptake, label='Nut Lim', color=cPNuts)
    ax2[i_plot].set_ylim(bottom=0)

    ax1[i_plot].set_ylabel('Maximum photosynth. rate [$d^{-1}$]')
    ax2[i_plot].set_ylabel('Phytoplankton Limiting terms')

    ax2[i_plot].legend(loc='lower right', fontsize='x-small')
    ax1[i_plot].legend()

    ax2[i_plot].set_xlabel('Time [days]')
    ax2[i_plot].set_xlim(1,365)
    ax2[i_plot].set_xticks([1, 90, 180, 270, 365])

def plotZoofluxes(outarray, pfn, zn, i_plot, title):

    # PLOTTING
    timedays_ly = timedays[1:366]
    # truncate outarraySiNO to last year of 5:
    outarray_ly = outarray[1460:1825]
    outarray_lyS = outarray[1459:1825]
    # color vectors
    colors = ['#808080','#d55e00', '#cc79a7', '#0072b2', '#009e73', '#009e73']
    alphas = [1., 0.8, 0.6, 0.4]
    lws = [1, 2.5, 4, 5.5]

    outindex = 4

    # zoo = [ZooGrowth[i] - ZooMortality[i] - ZooMixing[i] for i in range(zn)]
    # Zooplankton losses due to mortality and mixing

    #Zooplankton:
    Z = outarray_ly[:, 2]

    ZGains = getOUT(outarray_lyS, outindex+9)
    ZLinMort = getOUT(outarray_lyS, outindex+10)
    ZQuadMort = getOUT(outarray_lyS, outindex+11)
    ZMixing = getOUT(outarray_lyS, outindex+12)
    ZLosses = getOUT(outarray_lyS, outindex+13)

    ax1[i_plot].set_title(title)

    ax1[i_plot].plot(timedays_ly, Z, label='Z', color='grey')
    #ax1[i_plot].set_ylim(0, 0.7)

    ax2[i_plot].stackplot(timedays_ly, ZGains, labels=['Z Assimilated Grazing'], baseline='zero', colors=[cZGains])
    ax2[i_plot].stackplot(timedays_ly, -ZLinMort, -ZQuadMort, -ZMixing, labels=['Z Linear Mortality','Z Quadratic Mortality', 'Z Mixing'], baseline='zero', colors=[cZLinMort,cZQuadMort,cZMixing])
    ax2[i_plot].plot(timedays_ly, ZGains-ZLosses, label='total growth/loss-rate', color=cTOTALFLUX)
    #ax2[i_plot].set_ylim(-0.002,0.002)

    ax1[i_plot].set_ylabel('Zooplankton biomass [µM N]')
    ax2[i_plot].set_ylabel('Zooplankton growth/loss rates [µM N / day]')

    ax2[i_plot].legend(loc='lower right', fontsize='x-small')
    ax1[i_plot].legend()

    ax2[i_plot].set_xlabel('Time [days]')
    ax2[i_plot].set_xlim(1,365)
    ax2[i_plot].set_xticks([1, 90, 180, 270, 365])
    ax2[i_plot].set_ylim(-0.3,0.3)


# NITRATE     /// or Nitrogen????
f1, (ax1, ax2) = plt.subplots(2, 5, gridspec_kw = {'height_ratios':[1, 3]}, sharex='col')#, sharey='row')

plotNfluxes(outarray, 1, 1, 0, 'N')
#plotSifluxes(out1P1Z, 1, 1, 1, '1P1Z')
plotPhyfluxes(outarray, 1, 1, 1, 'P')
plotPhyLim(outarray, 1, 1, 2, 'P - growth limiting terms')
plotZoofluxes(outarray, 1, 1, 3, 'Z')
plotDetritusfluxes(outarray, 1, 1, 4, 'D')

f1.align_ylabels()

plt.subplots_adjust(hspace=0)
#plt.tight_layout()

plt.show()

#f1.savefig("fluxes01.png", dpi=1500)
