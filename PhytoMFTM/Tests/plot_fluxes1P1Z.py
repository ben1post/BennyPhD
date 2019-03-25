from Tests.runModel_2P2Zcheck import out1P1Z, timedays_model

import matplotlib.pyplot as plt

def plotNfluxes(outarray, pfn, zn, i_plot, title):

    # PLOTTING
    timedays_ly = timedays_model[1:366]
    # truncate outarraySiNO to last year of 5:
    outarray_lyS = outarray[1460:1825]
    # color vectors
    colors = ['#808080','#d55e00', '#cc79a7', '#0072b2', '#009e73', '#009e73']
    alphas = [1., 0.8, 0.6, 0.4]
    lws = [1, 2.5, 4, 5.5]

    # artist for legends
    FullArtist = plt.Line2D((0, 1), (0, 0), c=colors[4], alpha=alphas[1], lw=lws[0])

    outindex = 3+zn+pfn

    #Nitrate:
    N = outarray_lyS[:, 0]
    NRemin = outarray_lyS[:, outindex]
    NMixing = outarray_lyS[:, outindex + 1]

    DiffusiveMix = outarray_lyS[:, outindex]
    #Phyto:
    Gains = outarray_lyS[:, outindex + 5]



    ax1[i_plot].set_title(title)

    ax1[i_plot].plot(timedays_ly, N, label='N', color='grey')
    ax1[i_plot].set_ylim(bottom=0)

    ax2[i_plot].stackplot(timedays_ly, -Gains, labels=['Phyto Gains'], baseline='zero')
    ax2[i_plot].stackplot(timedays_ly, NRemin, NMixing, labels=['Remineralisation','Mixing'], baseline='zero')
    ax2[i_plot].plot(timedays_ly, NRemin+NMixing-Gains, label = 'Total Flux', color='black')
    ax2[i_plot].set_ylim(-0.09,0.09)


    ax1[i_plot].set_ylabel('Nitrate [µM N]')
    ax2[i_plot].set_ylabel('Nitrate Fluxes [µM N / day]')

    ax2[i_plot].legend(loc='lower right')
    ax1[i_plot].legend()

    ax2[i_plot].set_xlabel('Time [days]')
    ax2[i_plot].set_xlim(1,365)
    ax2[i_plot].set_xticks([1, 90, 180, 270, 365])

def plotSifluxes(outarray, pfn, zn, i_plot, title):

    # PLOTTING
    timedays_ly = timedays_model[1:366]
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
    timedays_ly = timedays_model[1:366]
    # truncate outarraySiNO to last year of 5:
    outarray_lyS = outarray[1460:1825]
    # color vectors
    colors = ['#808080','#d55e00', '#cc79a7', '#0072b2', '#009e73', '#009e73']
    alphas = [1., 0.8, 0.6, 0.4]
    lws = [1, 2.5, 4, 5.5]

    # artist for legends
    FullArtist = plt.Line2D((0, 1), (0, 0), c=colors[4], alpha=alphas[1], lw=lws[0])

    outindex = 3+zn+pfn

    #sum(ZooMortality) - NRemineralization - SiRemineralization - DetritusMixing + sum(UnassimilatedProduction) + sum([p[i].mortality(P[i]) for i in range(pfn)]) # Detritus

    #Detritus:
    Det = outarray_lyS[:, 2]
    ZooMortality = outarray_lyS[:, outindex + 14] #+
    NRemin = outarray_lyS[:, outindex] #-
    SiRemin = outarray_lyS[:, outindex + 2] #-
    DetMixing = outarray_lyS[:, outindex + 4] #-
    UnassimProd = outarray_lyS[:, outindex + 16] #+
    PMortality = outarray_lyS[:, outindex + 7] #+


    ax1[i_plot].set_title(title)

    ax1[i_plot].plot(timedays_ly, Det, label='Det', color='grey')
    ax1[i_plot].set_ylim(bottom=0)

    ax2[i_plot].stackplot(timedays_ly, -NRemin, -SiRemin, -DetMixing, labels=['N Remineralisation', 'Si Remineralisation', 'Mixing'], baseline='zero')
    ax2[i_plot].stackplot(timedays_ly, ZooMortality, UnassimProd, PMortality, labels=['Zoo Mortality', 'Unassimilated Feeding', 'Phyto Mortality'], baseline='zero')
    ax2[i_plot].plot(timedays_ly, -NRemin-SiRemin-DetMixing+ZooMortality+UnassimProd+PMortality, label = 'Total Flux', color='black')
    ax2[i_plot].set_ylim(-0.09,0.09)


    ax1[i_plot].set_ylabel('Detritus [µM N]')
    ax2[i_plot].set_ylabel('Detritus Fluxes [µM N / day]')

    ax2[i_plot].legend(loc='lower right')
    ax1[i_plot].legend()

    ax2[i_plot].set_xlabel('Time [days]')
    ax2[i_plot].set_xlim(1,365)
    ax2[i_plot].set_xticks([1, 90, 180, 270, 365])

def plotPhyfluxes(outarray, pfn, zn, i_plot, title):

    # PLOTTING
    timedays_ly = timedays_model[1:366]
    # truncate outarraySiNO to last year of 5:
    outarray_lyS = outarray[1460:1825]
    # color vectors
    colors = ['#808080','#d55e00', '#cc79a7', '#0072b2', '#009e73', '#009e73']
    alphas = [1., 0.8, 0.6, 0.4]
    lws = [1, 2.5, 4, 5.5]

    # artist for legends
    FullArtist = plt.Line2D((0, 1), (0, 0), c=colors[4], alpha=alphas[1], lw=lws[0])

    outindex = 3+zn+pfn

    # phy = [Gains[i] - PhytoGrazing[i] - PhytoMortality[i] - PhytoMixing[i] - PhytoSinking[i] for i in
    #       range(pfn)]  # Phytoplankton growth

    #Phytoplankton:
    P = sum([outarray_lyS[:, 3 + zn + i] for i in range(pfn)])
    Gains = outarray_lyS[:, outindex + 5]  # +

    PhytoMortality = outarray_lyS[:, outindex + 7]  # -
    PhytoGrazing = outarray_lyS[:, outindex + 8]  # -
    PhytoMixing = outarray_lyS[:, outindex + 9]  # -
    PhytoSinking = outarray_lyS[:, outindex + 10]  # -

    LightHrvstng = outarray_lyS[:, outindex + 11]  # +
    TempDepGro = outarray_lyS[:, outindex + 12]  # +

    TDGGains = Gains - P * Gains/TempDepGro
    LHGains = Gains - P * Gains/LightHrvstng
    NoLTGains = Gains - TDGGains - LHGains

    print(title,sum([sum(Gains),sum(-PhytoGrazing),sum(-PhytoMortality),sum(-PhytoMixing),sum(-PhytoSinking)]))

    ax1[i_plot].set_title(title)

    ax1[i_plot].plot(timedays_ly, P, label='P', color='grey')
    ax1[i_plot].set_ylim(bottom=0)

    ax2[i_plot].stackplot(timedays_ly, -PhytoMortality, -PhytoGrazing, -PhytoMixing, -PhytoSinking, labels=['Mortality', 'Grazing', 'Mixing', 'Sinking'], baseline='zero')
    ax2[i_plot].stackplot(timedays_ly, NoLTGains, TDGGains, LHGains, labels=['Nutrient Gains', 'Temp Dependency', 'Light Harvesting'], baseline='zero')
    ax2[i_plot].plot(timedays_ly, Gains-PhytoMortality-PhytoGrazing-PhytoMixing-PhytoSinking, label = 'Total Flux', color='black')
    ax2[i_plot].set_ylim(-0.09,0.14)


    ax1[i_plot].set_ylabel('Phytoplankton [µM N]')
    ax2[i_plot].set_ylabel('Phytoplankton Fluxes [µM N / day]')

    ax2[i_plot].legend(loc='lower right')
    ax1[i_plot].legend()

    ax2[i_plot].set_xlabel('Time [days]')
    ax2[i_plot].set_xlim(1,365)
    ax2[i_plot].set_xticks([1, 90, 180, 270, 365])

def plotZoofluxes(outarray, pfn, zn, i_plot, title):

    # PLOTTING
    timedays_ly = timedays_model[1:366]
    # truncate outarraySiNO to last year of 5:
    outarray_lyS = outarray[1460:1825]
    # color vectors
    colors = ['#808080','#d55e00', '#cc79a7', '#0072b2', '#009e73', '#009e73']
    alphas = [1., 0.8, 0.6, 0.4]
    lws = [1, 2.5, 4, 5.5]

    # artist for legends
    FullArtist = plt.Line2D((0, 1), (0, 0), c=colors[4], alpha=alphas[1], lw=lws[0])

    outindex = 3+zn+pfn

    # zoo = [ZooGrowth[i] - ZooMortality[i] - ZooMixing[i] for i in range(zn)]
    # Zooplankton losses due to mortality and mixing

    #Zooplankton:
    Z = sum([outarray_lyS[:, 3 + i] for i in range(zn)])

    ZooGrowth = outarray_lyS[:, outindex + 13]  # +
    ZooMortality = outarray_lyS[:, outindex + 14]  # -
    ZooMixing = outarray_lyS[:, outindex + 15]  # -

    print(title,sum([sum(ZooGrowth),sum(-ZooMixing),sum(-ZooMortality)]))

    ax1[i_plot].set_title(title)

    ax1[i_plot].plot(timedays_ly, Z, label='Z', color='grey')
    ax1[i_plot].set_ylim(0, 0.7)

    ax2[i_plot].stackplot(timedays_ly, ZooGrowth, labels=['Assimilated Grazing'], baseline='zero')
    ax2[i_plot].stackplot(timedays_ly, -ZooMortality, -ZooMixing, labels=['Mortality', 'Mixing'], baseline='zero')
    ax2[i_plot].plot(timedays_ly, ZooGrowth-ZooMortality-ZooMixing, label='Total Flux', color='black')
    ax2[i_plot].set_ylim(-0.002,0.002)

    ax1[i_plot].set_ylabel('Zooplankton [µM N]')
    ax2[i_plot].set_ylabel('Zooplankton Fluxes [µM N / day]')

    ax2[i_plot].legend(loc='lower right')
    ax1[i_plot].legend()

    ax2[i_plot].set_xlabel('Time [days]')
    ax2[i_plot].set_xlim(1,365)
    ax2[i_plot].set_xticks([1, 90, 180, 270, 365])


# NITRATE     /// or Nitrogen????
f1, (ax1, ax2) = plt.subplots(2, 5, gridspec_kw = {'height_ratios':[1, 3]}, sharex='col', sharey='row')

plotNfluxes(out1P1Z, 1, 1, 0, '1P1Z')
plotSifluxes(out1P1Z, 1, 1, 1, '1P1Z')
plotPhyfluxes(out1P1Z, 1, 1, 2, '1P1Z')
plotZoofluxes(out1P1Z, 1, 1, 3, '1P1Z')
plotDetritusfluxes(out1P1Z, 1, 1, 4, '1P1Z')

f1.align_ylabels()
plt.tight_layout()
plt.subplots_adjust(hspace=0)
#plt.show()

f1.savefig("fluxes01.png", dpi=1500)
