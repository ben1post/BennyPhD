from Tests.RunModel_2P2Zcheck import out1P1Z,out2P2Z

import matplotlib.pyplot as plt
import numpy as np

import matplotlib.gridspec as gridspec


def plotoutput(outarray, pfn, zn, i_plot):

    # PLOTTING
    timedays = timedays_model#[1:366]
    # truncate outarraySiNO to last year of 5:
    outarray_ly = outarray#[1460:1825]

    # color vectors
    colors = ['#808080','#d55e00', '#cc79a7', '#0072b2', '#009e73', '#009e73']
    alphas = [1., 0.8, 0.6, 0.4]
    lws = [1, 2.5, 4, 5.5]

    # artist for legends
    FullArtist = plt.Line2D((0, 1), (0, 0), c=colors[4], alpha=alphas[1], lw=lws[0])

    # Figure 1
    # N
    ax1[i_plot].plot(timedays, outarray_ly[:, 0], c=colors[1], lw=lws[0], alpha=alphas[0], label='Model')
    if i_plot == 0:
        ax1[i_plot].set_ylabel('Nitrate \n' '[µM]', multialignment='center', fontsize=10)
    ax1[i_plot].set_ylim(-0.1, 5)

    # Si
    ax2[i_plot].plot(timedays, outarray_ly[:, 1], c=colors[1], lw=lws[0], alpha=alphas[0])
    if i_plot == 0:
        ax2[i_plot].set_ylabel('Silicate \n' '[µM]', multialignment='center', fontsize=10)
    ax2[i_plot].set_ylim(-0.1, 12)


    # Phyto
    ax3[i_plot].plot(timedays, sum([outarray_ly[:, 3 + zn + i] for i in range(pfn)]), c=colors[4], lw=lws[1])
    [ax3[i_plot].plot(timedays, outarray_ly[:, 3 + zn + i], c=colors[i + 1]) for i in range(pfn)]
    if i_plot == 0:
        ax3[i_plot].set_ylabel('Phyto \n' '[µM N]', multialignment='center', fontsize=10)
    ax3[i_plot].set_ylim(-0.1, 0.8)

    ax3[i_plot].set_title('Phy Biomass & ChlA Data')

    # Z
    ax4[i_plot].plot(timedays, sum([outarray_ly[:, 3 + i] for i in range(zn)]), c=colors[4], lw=lws[1])
    [ax4[i_plot].plot(timedays, outarray_ly[:, 3 + i], c=colors[i + 1], lw=lws[0], alpha=alphas[0]) for i in range(zn)]
    if i_plot == 0:
        ax4[i_plot].set_ylabel('Zooplankton \n' '[µM N]', multialignment='center', fontsize=9)
    ax4[i_plot].tick_params('y', labelsize=10)

    ax4[i_plot].set_title('Zooplankton')
    ax4[i_plot].set_ylim(0, 0.62)

    # D
    ax5[i_plot].plot(timedays, outarray_ly[:, 3], c=colors[1], lw=lws[0], alpha=alphas[0])
    if i_plot == 0:
        ax5[i_plot].set_ylabel('Detritus \n' '[µM N]', multialignment='center', fontsize=9)

    ax5[i_plot].set_title('Detritus')
    ax5[i_plot].set_ylim(0,0.15)

    ax5[i_plot].set_xlabel('Day in year', fontsize=14)
    # Legend


def plotfluxes(outarray, pfn, zn, i_plot, title):

    # PLOTTING
    timedays_ly = timedays_model[1:366]
    # truncate outarraySiNO to last year of 5:
    outarray_lyS = outarray[1460:1825]
    #print(outarray)
    # color vectors
    colors = ['#808080','#d55e00', '#cc79a7', '#0072b2', '#009e73', '#009e73']
    alphas = [1., 0.8, 0.6, 0.4]
    lws = [1, 2.5, 4, 5.5]

    # artist for legends
    FullArtist = plt.Line2D((0, 1), (0, 0), c=colors[4], alpha=alphas[1], lw=lws[0])
    """
    y0 = NRemineralization + NMixing   # Nitrate upwelling and remineralisation
    y1 = SiRemineralization + SiMixing # Silicate upwelling and remineralisation
    y2 = sum(ZooMortality) - NRemineralization - SiRemineralization - DetritusMixing


    y0 = y0 - sum([P[i ] *Gains[i] for i in range(pfn)]) # Nitrate drawdown
    y1 = y1 - sum([p[i].silicatedrawdown(P[i] ,Gains[i]) for i in range(pfn)]) # Silicate drawdown
    y2 = y2 + sum(UnassimilatedProduction) + sum([p[i].mortality(P[i]) for i in range(pfn)]) # Detritus

    Gains = [p[i].gains(N_Uptake[i], Si_Uptake[i], LightHarvesting[i], TemperatureDepGrowth)
             for i in range(pfn)]
    Losses = [p[i].losses(int_MLD, PhytoGrazed[i], DiffusiveMixing) for i in range(pfn)]

    zoo = [ZooGrowth[i ] -ZooMortality[i ] -ZooMixing[i] for i in range(zn)]    # Zooplankton losses due to mortality and mixing
    phy = [P[i] * (Gains[i] - Losses[i]) for i in range(pfn)]  # Phytoplankton growth

    outputlist = [0 DiffusiveMixing, 1 ActiveMixing, 2 NRemineralization, 3 SiRemineralization,
                  4 DetritusMixing, 5 NMixing, 6 SiMixing, 7 sum(ZooMortality), 8 sum(N_Uptake), 9 sum(Si_Uptake),
                  10 sum(LightHarvesting), 11 TemperatureDepGrowth, 12 sum(Gains), 13 sum(PhytoGrazed), 14 sum(ZooGrazeGrowth),
                  15 sum(Losses), 16 sum(ZooMixing), 17 sum(ZooGrowth), 18 sum(UnassimilatedProduction)]
    outputlist[0] = N         - outputlist[0]
    outputlist[1] = ActiveMixing            - outputlist[1]
    outputlist[2] = NRemineralization     - outputlist[2]
    outputlist[3] = Si*SiRemineralization   - outputlist[3]
    outputlist[4] = NMixing                 - outputlist[4]
    outputlist[5] = SiMixing                - outputlist[5]
    outputlist[6] = sum(ZooMortality)       - outputlist[6]
    outputlist[7] = sum(N_Uptake)           - outputlist[7]
    outputlist[8] = sum(Si_Uptake)          - outputlist[8]
    outputlist[9] = sum(LightHarvesting)    - outputlist[9]
    outputlist[10] = TemperatureDepGrowth   - outputlist[10]
    outputlist[11] = sum([P[i ] *Gains[i] for i in range(pfn)])             - outputlist[11]
    outputlist[12] = sum(PhytoGrazed)       - outputlist[12]
    outputlist[13] = sum(ZooItots)          - outputlist[13]
    outputlist[14] = sum(Losses)            - outputlist[14]
    outputlist[15] = sum(ZooMixing)         - outputlist[15]
    outputlist[16] = sum(ZooGrowth)         - outputlist[16]
    outputlist[17] = sum(UnassimilatedProduction) - outputlist[17]
    outputlist[18] = sum(ZooMortality)      - outputlist[18]
    out = [y0, y1, y2] + zoo + phy + outputlist
    """
    outindex = 3+zn+pfn

    #Nitrate:
    N = outarray_lyS[:, outindex]
    NRemin = outarray_lyS[:, outindex + 2]
    NMixing = outarray_lyS[:, outindex + 4]

    DiffusiveMix = outarray_lyS[:, outindex]
    #Phyto:
    Psum = sum([outarray_lyS[:, 3 + zn + i] for i in range(pfn)])
    Gains = outarray_lyS[:, outindex + 11]


    Labels1 = ['PGains','NRemin','NMixing']

    #ax1[i_plot].stackplot(timedays_ly, -PGains, NRemin, NMixing, labels=Labels1)

    ax1[i_plot].set_title(title)

    ax1[i_plot].plot(timedays_ly, N, label='N', color='grey')
    ax1[i_plot].set_ylim(bottom=0)

    ax2[i_plot].stackplot(timedays_ly, -Gains, labels=['Phyto Gains'], baseline='zero')
    ax2[i_plot].stackplot(timedays_ly, NRemin, NMixing, labels=['Remineralisation','Mixing'], baseline='zero')
    ax2[i_plot].plot(timedays_ly, NRemin+NMixing-Gains, label = 'Total Flux', color='black')
    ax2[i_plot].set_ylim(-0.09,0.09)

    if i_plot == 0:
        ax1[i_plot].set_ylabel('Nitrate [µM N]')
        ax2[i_plot].set_ylabel('Nitrate Fluxes [µM N / day]')
    if i_plot == 1:
        ax2[i_plot].legend(loc='lower right')

    ax2[i_plot].set_xlabel('Time [days]')
    ax2[i_plot].set_xlim(1,365)
    ax2[i_plot].set_xticks([1, 90, 180, 270, 365])
    Labels2 = ['PGains']

    #ax2[i_plot].plot(Psum)#, labels=Labels2)
    #ax2[i_plot].legend(loc='lower right')


# plt.subplots_adjust(hspace=0.01)
#f1, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 2, sharex='col', sharey='row')


f1, (ax1, ax2) = plt.subplots(2, 2, gridspec_kw = {'height_ratios':[1, 3]}, sharex='col', sharey='row')

timedays_model = np.arange(0., 5 * 365., 1.0)


plotfluxes(out1P1Z, 1, 1, 0, '1P1Z')
#plotoutput(out2P2Z, 2, 2, 1)
#plotoutput(out3P3Z, 3, 3, 2)
#plotoutput(out4P4Z, 4, 4, 3)
plotfluxes(out2P2Z, 2, 2, 1, '2P2Z')

# f1.set_figheight(15)
f1.align_ylabels()
plt.tight_layout()
plt.subplots_adjust(hspace=0)
plt.show()

#print(out1P1Z)
# f1.savefig("foo2.pdf", bbox_inches='tight')

