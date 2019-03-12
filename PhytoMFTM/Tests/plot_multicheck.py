from Tests.RunModel_2P2Zcheck import out1P1Z,out5P5Z #out2P2Z,out3P3Z,out4P4Z,

import matplotlib.pyplot as plt
import numpy as np


def plotoutput(outarray, pfn, zn, i_plot):

    # PLOTTING
    timedays = timedays_model#[1:366]
    # truncate outarraySiNO to last year of 5:
    outarray_ly = outarray#[1460:1825]

    # color vectors
    #colors = ['#edc951', '#dddddd', '#00a0b0', '#343436', '#cc2a36']
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




def plotfluxes(outarray, pfn, zn, i_plot):

    # PLOTTING
    timedays_ly = timedays_model[1:366]
    # truncate outarraySiNO to last year of 5:
    outarray_lyS = outarray[1460:1825]

    # color vectors
    #colors = ['#edc951', '#dddddd', '#00a0b0', '#343436', '#cc2a36']
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

    outputlist = [DiffusiveMixing, ActiveMixing, NRemineralization, SiRemineralization,
                  DetritusMixing, NMixing, SiMixing, sum(ZooMortality), sum(N_Uptake), sum(Si_Uptake),
                  sum(LightHarvesting), TemperatureDepGrowth, sum(Gains), sum(PhytoGrazed), sum(ZooGrazeGrowth),
                  sum(Losses), sum(ZooMixing), sum(ZooGrowth), sum(UnassimilatedProduction)]

    out = [y0, y1, y2] + zoo + phy + outputlist
    """
    outindex = 3+zn+pfn

    #Nitrate:
    N = outarray_lyS[:, 0]
    NRemin = outarray_lyS[:, outindex + 2]
    NMixing = outarray_lyS[:, outindex + 5]

    #Phyto:
    Psum = sum([outarray_lyS[:, 3 + zn + i] for i in range(pfn)])
    Gains = outarray_lyS[:, outindex + 12]
    PGains = Psum * Gains

    print(Gains)

    Labels = ['NRemin','NMixing','PGains']

    ax1[i_plot].stackplot(timedays_ly,  -PGains,NRemin*N,NMixing*N, labels=Labels)
    ax1[i_plot].legend(loc='lower right')

    """
    # Figure 1
    # N
    ax1[i_plot].plot(timedays_ly, outarray_lyS[:, 0], c=colors[1], lw=lws[0], alpha=alphas[0], label='Model')
    if i_plot == 0:
        ax1[i_plot].set_ylabel('Nitrate \n' '[µM]', multialignment='center', fontsize=10)
    ax1[i_plot].set_ylim(-0.1, 5)

    # Si
    ax2[i_plot].plot(timedays_ly, outarray_lyS[:, 1], c=colors[1], lw=lws[0], alpha=alphas[0])
    if i_plot == 0:
        ax2[i_plot].set_ylabel('Silicate \n' '[µM]', multialignment='center', fontsize=10)
    ax2[i_plot].set_ylim(-0.1, 12)


    # Phyto
    ax3[i_plot].plot(timedays_ly, sum([outarray_lyS[:, 3 + zn + i] for i in range(pfn)]), c=colors[4], lw=lws[1])
    [ax3[i_plot].plot(timedays_ly, outarray_lyS[:, 3 + zn + i], c=colors[i + 1]) for i in range(pfn)]
    if i_plot == 0:
        ax3[i_plot].set_ylabel('Phyto \n' '[µM N]', multialignment='center', fontsize=10)
    ax3[i_plot].set_ylim(-0.1, 0.8)

    ax3[i_plot].set_title('Phy Biomass & ChlA Data')

    # Z
    ax4[i_plot].plot(timedays_ly, sum([outarray_lyS[:, 3 + i] for i in range(zn)]), c=colors[4], lw=lws[1])
    [ax4[i_plot].plot(timedays_ly, outarray_lyS[:, 3 + i], c=colors[i + 1], lw=lws[0], alpha=alphas[0]) for i in range(zn)]
    if i_plot == 0:
        ax4[i_plot].set_ylabel('Zooplankton \n' '[µM N]', multialignment='center', fontsize=9)
    ax4[i_plot].tick_params('y', labelsize=10)

    ax4[i_plot].set_title('Zooplankton')
    ax4[i_plot].set_ylim(0, 0.62)

    # D
    ax5[i_plot].plot(timedays_ly, outarray_lyS[:, 3], c=colors[1], lw=lws[0], alpha=alphas[0])
    if i_plot == 0:
        ax5[i_plot].set_ylabel('Detritus \n' '[µM N]', multialignment='center', fontsize=9)

    ax5[i_plot].set_title('Detritus')
    ax5[i_plot].set_ylim(0,0.15)

    ax5[i_plot].set_xlabel('Day in year', fontsize=14)
    # Legend
    """


# plt.subplots_adjust(hspace=0.01)
#f1, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 2, sharex='col', sharey='row')
f1, (ax1, ax2) = plt.subplots(2, 2, sharex='col', sharey='row')


timedays_model = np.arange(0., 5 * 365., 1.0)


plotfluxes(out1P1Z, 1, 1, 0)
#plotoutput(out2P2Z, 2, 2, 1)
#plotoutput(out3P3Z, 3, 3, 2)
#plotoutput(out4P4Z, 4, 4, 3)
plotfluxes(out5P5Z, 5, 5, 1)

# f1.set_figheight(15)
# plt.tight_layout()
plt.show()

# f1.savefig("foo2.pdf", bbox_inches='tight')

