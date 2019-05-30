#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import pandas

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.mplot3d import Axes3D

from Tests.MODELfunctions_plotParams import all_params, all_params21, all_params22
from PhytoMFTM.ModelClasses import Plankton
from lmfit import Parameters, Parameter

colormap = pandas.DataFrame({"name" : ["MLD", "Si", "N", "Phyto", "MikroZ", "MesoZ", "D", "SST", "PAR"],
                             "color" : ["#1D71B8", "#FFDE00", "#575756", "#009640", "#EA5B0C", "#BE1622", "#B17F4A", "#E94E1B", "#F9B233"]})

colormap["name"] = colormap["name"].apply(lambda x: x.lower())
c = dict(zip(*colormap.values.T))
mcolors.get_named_colors_mapping().update(c)

# WHICH FUNCTIONS TO PLOT?

# 1. Light Harvesting via Steele
#   - PI curve
# 2. Light Harvesting via Smith
#   - PI curve
#   - x_H at different depth and different phytoplankton biomass
#   - integrated over various depth
# 3. Feeding Functions
#   - grazing at constant Z, varying P (2 types, changing relative concentrations)
#   - plot inter ZOO grazing as func of Z1 & Z2
# 4. Nutrient Uptake for each functional type
#   - simple curves N conc vs Uptake --> all plotted together
# 5. Nutrient mixing
#   - with variable N_0
#   - with constant N_0
#   - with current MLD forcing
#   - with smoothed out MLD forcing

# LET'S GET STARTED

print(list(all_params)[:])

z = Plankton(all_params, 'Zooplankton').init()
p = Plankton(all_params, 'Phytoplankton').init()


z2 = Plankton(all_params21, 'Zooplankton').init()
p2 = Plankton(all_params21, 'Phytoplankton').init()


z3 = Plankton(all_params22, 'Zooplankton').init()
p3 = Plankton(all_params22, 'Phytoplankton').init()

def steelefunc(I):
    # simple definition fo steeles PI curve function
    PM = 1
    Iopt = 30
    pi = PM * I / Iopt * np.exp(1 - I / Iopt)
    return pi

def smithfunc(I):
    # simple definition fo steeles PI curve function
    alpha = 0.15
    VpMax = 1
    pi = alpha * I * VpMax / np.sqrt(VpMax**2 + alpha**2 * I**2)
    return pi

def PIcurves():
    #f1, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex='col', sharey='row')
    # N / Si / Pdt / Pc / Pdn / Pn / Zmu / Zlambda / D
    #ax1.set_title('Steele PI curve')
    PAR1 = np.linspace(0.1, 60, 50)
    smith1 = smithfunc(PAR1)
    steele1 = steelefunc(PAR1)
    # Figure 1
    # N
    #ax1.plot(intMLD, intPAR, c="N", lw=lws[0], alpha=alphas[0], label='Model')
    # Legend
    #f1.align_ylabels()
    fig1, ax1 = plt.subplots()
    ax1.plot(PAR1,smith1, label= 'smith')
    ax1.plot(PAR1,steele1, label='steele')
    #plt.margins(x=0)
    ax1.set_xlabel('$PAR_0$ [E $m^{−2}$ $s^{−1}$]')
    ax1.set_ylabel('rate of photosynthetis')
    ax1.set_title('PI curves')
    fig1.tight_layout()
    fig1.legend()
    plt.savefig('PIcurves.png')
    fig1.show()

def steelePI():
    #f1, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex='col', sharey='row')
    # N / Si / Pdt / Pc / Pdn / Pn / Zmu / Zlambda / D
    #ax1.set_title('Steele PI curve')
    intMLD = np.linspace(0.1, 60, 50)
    intPAR = np.linspace(0.1, 60, 50)
    Steele = np.array([p[0].lightharvesting(i, j) for i in intMLD for j in intPAR])
    print(Steele)
    # Figure 1
    # N
    #ax1.plot(intMLD, intPAR, c="N", lw=lws[0], alpha=alphas[0], label='Model')
    X, Y = np.meshgrid(intMLD, intPAR)
    Z = Steele.reshape(50, 50)
    # Legend
    #f1.align_ylabels()
    plt.margins(x=0)
    plt.tight_layout()

    fig, ax = plt.subplots()
    cs = ax.contourf(X, Y, Z, 50)
    fig.colorbar(cs)
    ax.set_xlabel('$PAR_0$ [E $m^{−2}$ $s^{−1}$]')
    ax.set_ylabel('MLD depth [m]')
    ax.set_title('STEELE at kw = 0.1')
    plt.savefig('SteelePI.png')
    #plt.savefig('FirstNaiveOutputCARIACO.png')
    fig.show()

def smithPI():
    #f1, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex='col', sharey='row')
    # N / Si / Pdt / Pc / Pdn / Pn / Zmu / Zlambda / D
    #ax1.set_title('Steele PI curve')
    intMLD = np.linspace(0.1, 60, 50)
    intPAR = np.linspace(0.1, 60, 50)
    SMITH = np.array([p[0].smithpi(i, j, [1]) for i in intMLD for j in intPAR])
    print(len(SMITH))
    print(SMITH)
    # Figure 1
    # N
    #ax1.plot(intMLD, intPAR, c="N", lw=lws[0], alpha=alphas[0], label='Model')
    X, Y = np.meshgrid(intMLD, intPAR)
    Z = SMITH.reshape(50, 50)
    # Legend
    #f1.align_ylabels()
    plt.margins(x=0)
    plt.tight_layout()

    fig3, ax3 = plt.subplots()
    cs = ax3.contourf(X, Y, Z, 50)
    fig3.colorbar(cs)
    ax3.set_xlabel('$PAR_0$ [E $m^{−2}$ $s^{−1}$]')
    ax3.set_ylabel('MLD depth [m]')
    ax3.set_title('SMITH at kw = 0.04 and P = 1 µM N')
    #plt.savefig('FirstNaiveOutputCARIACO.png')
    plt.savefig('SmithPI.png')
    fig3.show()


def PhytoGrazed1():
    # f1, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex='col', sharey='row')
    # N / Si / Pdt / Pc / Pdn / Pn / Zmu / Zlambda / D
    P1 = np.linspace(0.01, 1., 50)
    P2 = np.linspace(0.01, 1., 50)

    def f(P1, P2):
        Frho = P1 ** 2 + P2 ** 2  # active switching coefficient
        Fp = P1 + P2  # total food available
        # now calculate grazing probability of either Zoo-Type
        Fq = P1 ** 2

        # print('p1', P[0], 'p2', P[1], 'Frho', Frho, 'Fp', Fp, 'Fq/Frho', Fq/Frho)

        Gj = 1 * (1 / Frho) * ((Fp ** 2) / (0.5 ** 2 + Fp ** 2))
        GjJ = Gj * Fq
        return GjJ

    # z.grazinprobability, 2 Z -> loops across all 4 phy
    # p.zoograzing2, 4 P -> loops across all 2 zoo
    # or
    # z.intake 2 Z -> loops across all 4 phy

    # Gj1 = np.array([z[0].grazingprobability(i) for i in Z1])

    # Gj2 = [z[1].grazingprobability([0, P2[i], P3[j], 0]) for i in range(len(P2)) for j in range(len(P3))]

    #Gj = [[z2[0].grazingprobability([i, j]), i, j] for i in P1 for j in P2]
    #gj1 = [p2[0].zoograzing2([gj], i) for gj,i,j in Gj]

    X, Y = np.meshgrid(P1, P2)
    Z = f(X, Y)
    #P1Gr = np.array(gj1)
    #Z = P1Gr.reshape(50, 50)

    fig3, ax3 = plt.subplots()
    cs = ax3.contourf(X, Y, Z, 50)
    cbar = fig3.colorbar(cs)
    cbar.set_ticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
    ax3.set_xlabel('P1 [µM N]')
    ax3.set_ylabel('P2 [µM N]')
    ax3.set_title('Grazing on P1 [µM $m^{-3} d^{-1}$]')
    plt.savefig('PhytoGrazed.png')
    fig3.show()

def ZooIntake1():
    #f1, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex='col', sharey='row')
    # N / Si / Pdt / Pc / Pdn / Pn / Zmu / Zlambda / D
    Z1 = 1
    P1 = np.linspace(0.01, 1., 50)
    P2 = np.linspace(0.01, 1., 50)

    # z.grazinprobability, 2 Z -> loops across all 4 phy
    # p.zoograzing2, 4 P -> loops across all 2 zoo
    # or
    # z.intake 2 Z -> loops across all 4 phy


    Z1Gr = np.array([z2[0].zoointake(z2[0].grazingprobability([i, j]), [i,j]) for i in P1 for j in P2])
    #P1Gr = np.array([p2[0].zoograzing2([z2[0].grazingprobability([i,j])], i, [1,1]) for i in P1 for j in P2])

    print(Z1Gr)
    # Figure 1
    # N
    #ax1.plot(intMLD, intPAR, c="N", lw=lws[0], alpha=alphas[0], label='Model')
    X, Y = np.meshgrid(P2, P1)
    Z = Z1Gr.reshape(50, 50)
    # Legend
    #f1.align_ylabels()
    #plt.margins(x=0)
    #plt.tight_layout()

    fig4, ax4 = plt.subplots()
    cs = ax4.contourf(X, Y, Z, 50)
    cbar = fig4.colorbar(cs)
    cbar.set_ticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
    ax4.set_xlabel('P1 [µM N]')
    ax4.set_ylabel('P2 [µM N]')
    ax4.set_title('Grazing of 1 Z [µM $m^{-3} d^{-1}$]')
    plt.savefig('ZooIntake1.png')
    fig4.show()

def grazingvaryingP1():
    #f1, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex='col', sharey='row')
    # N / Si / Pdt / Pc / Pdn / Pn / Zmu / Zlambda / D
    #ax1.set_title('Steele PI curve')
    P1 = np.linspace(0.0001, 2, 50)
    P2 = .5

    def f(P1, P2):
        Frho = P1 ** 2 + P2 ** 2  # active switching coefficient
        Fp = P1 + P2  # total food available
        # now calculate grazing probability of either Zoo-Type
        Fq = P1 ** 2

        # print('p1', P[0], 'p2', P[1], 'Frho', Frho, 'Fp', Fp, 'Fq/Frho', Fq/Frho)
        MM = ((Fp ** 2) / (0.5 ** 2 + Fp **2))
        Gj = 1 * (Fq / Frho) * ((Fp ** 2) / (1 ** 2  + Fp ** 2))
        #GjJ = Gj * Fq
        return [Gj, MM, Frho, Fp, Fq]

    GrazingProb = z2[0].grazingprobability([P1, P2])
    IntakeP1 = p2[0].zoograzing2([GrazingProb], P1)
    IntakeP2 = p2[1].zoograzing2([GrazingProb], P2)

    IntakeP1 = f(P1,P2)[0]
    IntakeP2 = f(P2,P1)[0]
    Total = IntakeP1 + IntakeP2

    OtherOut = f(P1,P2)

    # Figure 1
    # N
    #ax1.plot(intMLD, intPAR, c="N", lw=lws[0], alpha=alphas[0], label='Model')
    # Legend
    #f1.align_ylabels()
    fig5, (ax5, ax6) = plt.subplots(2,1)
    ax5.plot(P1, IntakeP1, label='P1 Intake')
    ax5.plot(P1, IntakeP2, label='P2 Intake')
    ax5.plot(P1, Total, label='Total')

    #plt.margins(x=0)
    ax5.set_xlabel('P1 [µM N]')
    ax5.set_ylabel('Intake [µM $m^{-3} d^{-1}$]')
    ax5.set_title('Grazing at P2 = 0.5 µM N')


    ax6.plot(P1, OtherOut[1], label='MM')
    ax6.plot(P1, OtherOut[2], label='Frho')
    ax6.plot(P1, OtherOut[3], label='Fp')
    ax6.plot(P1, OtherOut[4], label='Fq')

    #plt.savefig('FirstNaiveOutputCARIACO.png')
    fig5.tight_layout()
    ax5.legend()
    ax6.legend()
    plt.savefig('grazingvaryingP1.png')
    fig5.show()


def FASHAMgrazingvaryingP1():
    #f1, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex='col', sharey='row')
    # N / Si / Pdt / Pc / Pdn / Pn / Zmu / Zlambda / D
    #ax1.set_title('Steele PI curve')
    P1 = np.linspace(0.0001, 2, 50)
    P2 = 1.5

    def f(P1, P2):
        Frho = P1 ** 2 + P2 ** 2  # active switching coefficient
        Fp = P1 + P2  # total food available
        # now calculate grazing probability of either Zoo-Type
        Fq = P1 ** 2

        # print('p1', P[0], 'p2', P[1], 'Frho', Frho, 'Fp', Fp, 'Fq/Frho', Fq/Frho)
        MM = ((Fp ** 2) / (1 * Fp + Fp ** 2))
        Gj = 1 * ((Fq ** 2) / (1 * Fp + Fp ** 2))
        #GjJ = Gj * Fq
        return [Gj, MM, Frho, Fp, Fq]

    #GrazingProb = z2[0].grazingprobability([P1, P2])
    #IntakeP1 = p2[0].zoograzing2([GrazingProb], P1)
    #IntakeP2 = p2[1].zoograzing2([GrazingProb], P2)

    IntakeP1 = f(P1, P2)[0]
    IntakeP2 = f(P2, P1)[0]
    print(IntakeP1,IntakeP2)
    Total = IntakeP1 + IntakeP2

    OtherOut = f(P1,P2)

    # Figure 1
    # N
    #ax1.plot(intMLD, intPAR, c="N", lw=lws[0], alpha=alphas[0], label='Model')
    # Legend
    #f1.align_ylabels()
    fig5, (ax5, ax6) = plt.subplots(2,1)
    ax5.plot(P1, IntakeP1, label='P1 Intake')
    ax5.plot(P1, IntakeP2, label='P2 Intake')
    ax5.plot(P1, Total, label='Total')

    #plt.margins(x=0)
    ax5.set_xlabel('P1 [µM N]')
    ax5.set_ylabel('Intake [µM $m^{-3} d^{-1}$]')
    ax5.set_title('Grazing at P2 = 0.5 µM N')


    ax6.plot(P1, OtherOut[1], label='MM')
    ax6.plot(P1, OtherOut[2], label='Frho')
    ax6.plot(P1, OtherOut[3], label='Fp')
    ax6.plot(P1, OtherOut[4], label='Fq')

    ax5.set_title('FASHAM')
    #plt.savefig('FirstNaiveOutputCARIACO.png')
    fig5.tight_layout()
    ax5.legend()
    ax6.legend()
    #plt.savefig('grazingvaryingP1.png')
    fig5.show()


def ANDERSONgrazingvaryingP1():
    #f1, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex='col', sharey='row')
    # N / Si / Pdt / Pc / Pdn / Pn / Zmu / Zlambda / D
    #ax1.set_title('Steele PI curve')
    P1 = np.linspace(0.0001, 2, 50)
    P2 = 1.5

    def f(P1, P2):
        Imax = 1
        rhoP1 = 0.5
        rhoP2 = 0.5
        kZ = .5
        Z = 1
        Gj = (Imax * rhoP1 * P1**2)/(kZ**2 + rhoP1 * P1**2 + rhoP2 * P2**2) * Z


        Frho = P1 ** 2 + P2 ** 2  # active switching coefficient
        Fp = P1 + P2  # total food available
        # now calculate grazing probability of either Zoo-Type
        Fq = P1 ** 2

        # print('p1', P[0], 'p2', P[1], 'Frho', Frho, 'Fp', Fp, 'Fq/Frho', Fq/Frho)
        MM = ((Fp ** 2) / (1 * Fp + Fp ** 2))
        #Gj = 1 * ((Fq ** 2) / (1 * Fp + Fp ** 2))
        #GjJ = Gj * Fq
        return [Gj, MM, Frho, Fp, Fq]

    #GrazingProb = z2[0].grazingprobability([P1, P2])
    #IntakeP1 = p2[0].zoograzing2([GrazingProb], P1)
    #IntakeP2 = p2[1].zoograzing2([GrazingProb], P2)

    IntakeP1 = f(P1, P2)[0]
    IntakeP2 = f(P2, P1)[0]
    print(IntakeP1,IntakeP2)
    Total = IntakeP1 + IntakeP2

    OtherOut = f(P1,P2)

    # Figure 1
    # N
    #ax1.plot(intMLD, intPAR, c="N", lw=lws[0], alpha=alphas[0], label='Model')
    # Legend
    #f1.align_ylabels()
    fig5, (ax5, ax6) = plt.subplots(2,1)
    ax5.plot(P1, IntakeP1, label='P1 Intake')
    ax5.plot(P1, IntakeP2, label='P2 Intake')
    ax5.plot(P1, Total, label='Total')

    #plt.margins(x=0)
    ax5.set_xlabel('P1 [µM N]')
    ax5.set_ylabel('Intake [µM $m^{-3} d^{-1}$]')
    ax5.set_title('Grazing at P2 = 0.5 µM N')


    ax6.plot(P1, OtherOut[1], label='MM')
    ax6.plot(P1, OtherOut[2], label='Frho')
    ax6.plot(P1, OtherOut[3], label='Fp')
    ax6.plot(P1, OtherOut[4], label='Fq')

    ax5.set_title('ANDERSON')
    #plt.savefig('FirstNaiveOutputCARIACO.png')
    fig5.tight_layout()
    ax5.legend()
    ax6.legend()
    #plt.savefig('grazingvaryingP1.png')
    fig5.show()


def grazingcomparison():
    #f1, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex='col', sharey='row')
    # N / Si / Pdt / Pc / Pdn / Pn / Zmu / Zlambda / D
    #ax1.set_title('Steele PI curve')
    P1 = np.linspace(0.0001, 2, 50)
    P2 = 1.5

    def anderson(P1, P2):

        Imax = 1
        rhoP1 = 0.5
        rhoP2 = 0.5
        kZ = .5
        Z = 1

        Gj = (Imax * rhoP1 * P1**2)/(kZ**2 + rhoP1 * P1**2 + rhoP2 * P2**2) * Z
        return Gj

    def fasham(P1, P2):

        Vmax = 1
        ksat = .5
        rhoP1 = .5
        rhoP2 = .5

        Frho = rhoP1 * P1 ** 2 + rhoP2 * P2 ** 2  # active switching coefficient
        Fp = rhoP1 * P1 + rhoP2 * P2  # total food available
        # now calculate grazing probability of either Zoo-Type
        Fq = rhoP1 * P1 ** 2

        Gj = Vmax * (Fq / (ksat * Fp + Frho))

        return Gj


    def vallina(P1, P2):

        Vmax = 1
        ksat = .5
        beta = 2
        rhoP1 = 0.5
        rhoP2 = 0.5

        Frho = rhoP1 * P1 ** 2 + rhoP2 * P2 ** 2  # active switching coefficient
        Fp = rhoP1 * P1 + rhoP2 * P2  # total food available
        # now calculate grazing probability of either Zoo-Type
        Fq = rhoP1 * P1 ** 2

        Gj = Vmax * (Fq / Frho) * ((Fp ** beta) / (ksat ** beta + Fp ** beta))
        return Gj

    #GrazingProb = z2[0].grazingprobability([P1, P2])
    #IntakeP1 = p2[0].zoograzing2([GrazingProb], P1)
    #IntakeP2 = p2[1].zoograzing2([GrazingProb], P2)

    IP1Vallina = vallina(P1, P2)
    IP2Vallina = vallina(P2, P1)
    TotIPVallina = IP1Vallina + IP2Vallina

    IP1Anderson = anderson(P1, P2)
    IP2Anderson = anderson(P2, P1)
    TotIPAnderson = IP1Anderson + IP2Anderson

    IP1Fasham = fasham(P1, P2)
    IP2Fasham = fasham(P2, P1)
    TotIPFasham = IP1Fasham + IP2Fasham


    # Figure 1
    # N
    #ax1.plot(intMLD, intPAR, c="N", lw=lws[0], alpha=alphas[0], label='Model')
    # Legend
    #f1.align_ylabels()
    fig5, (ax5, ax6, ax7) = plt.subplots(3,1, sharex='col', sharey='row')
    ax5.plot(P1, IP1Vallina, label='P1 Intake')
    ax5.plot(P1, IP2Vallina, label='P2 Intake')
    ax5.plot(P1, TotIPVallina, label='Total')

    #plt.margins(x=0)
    #ax5.set_xlabel('P1 [µM N]')
    ax5.set_ylabel('Intake [µM $m^{-3} d^{-1}$]')
    ax5.set_title('VALLINA - Grazing at P2 = 0.5 µM N')

    ax6.plot(P1, IP1Anderson, label='P1 Intake')
    ax6.plot(P1, IP2Anderson, label='P2 Intake')
    ax6.plot(P1, TotIPAnderson, label='Total')

    #ax6.set_xlabel('P1 [µM N]')
    ax6.set_ylabel('Intake [µM $m^{-3} d^{-1}$]')
    ax6.set_title('ANDERSON')

    ax7.plot(P1, IP1Fasham, label='P1 Intake')
    ax7.plot(P1, IP2Fasham, label='P2 Intake')
    ax7.plot(P1, TotIPFasham, label='Total')

    ax7.set_xlabel('P1 [µM N]')
    ax7.set_ylabel('Intake [µM $m^{-3} d^{-1}$]')
    ax7.set_title('FASHAM')

    #plt.savefig('FirstNaiveOutputCARIACO.png')
    ax5.legend()
    ax6.legend()
    ax7.legend()
    #plt.savefig('grazingvaryingP1.png')
    fig5.tight_layout()
    fig5.show()


def GrazingTry2():
    # Axes3D import has side effects, it enables using projection='3d' in add_subplot

    fig6 = plt.figure()
    ax6 = fig6.add_subplot(111, projection='3d')
    P1x = P2y = np.arange(.01, 3.0, 0.01)
    P1, P2 = np.meshgrid(P1x, P2y)

    zs = np.array(p[0].zoograzing2(z[0].grazingprobability([np.ravel(P1), np.ravel(P2)]), np.ravel(P1)))
    Z = zs.reshape(P1.shape)

    ax6.plot_surface(P1, P2, Z)

    ax6.set_xlabel('P1 [µM N]')
    ax6.set_ylabel('P2 [µM N]')
    ax6.set_zlabel('Ingestion P1')
    #plt.savefig('GrazingTry2.png')
    fig6.show()



def PhytoGrazed2():
    # f1, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex='col', sharey='row')
    # N / Si / Pdt / Pc / Pdn / Pn / Zmu / Zlambda / D
    Z1 = 1
    P1 = np.linspace(0.01, 2., 50)
    P2 = np.linspace(0.01, 3., 50)

    Gj = [[z2[0].zoofeeding([i, j], func='vallina'), i, j] for i in P1 for j in P2]
    gj1 = [p2[0].zoograzing([gj], i, [Z1], j) for gj,i,j in Gj]
    Pone = [Gj[x][1] for x in range(50*50)]
    Pone2 = np.array(Pone)
    Pone3 = Pone2.reshape(50,50)
    Ptwo = [Gj[x][2] for x in range(50 * 50)]
    Ptwo2 = np.array(Ptwo)
    Ptwo3 = Ptwo2.reshape(50, 50)
    print(Pone3)
    X, Y = np.meshgrid(P1, P2)
    print(X)
    #Z = f(X, Y)
    P1Gr = np.array(gj1)
    Z = P1Gr.reshape(50, 50)

    fig3, ax3 = plt.subplots()
    cs = ax3.contourf(Pone3, Ptwo3, Z, 50)
    cbar = fig3.colorbar(cs)
    cbar.set_ticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
    ax3.set_xlabel('P1 [µM N]')
    ax3.set_ylabel('P2 [µM N]')
    ax3.set_title('Grazing on P1 [µM $d^{-1}$]')
    plt.savefig('PhytoGrazed.png')
    fig3.show()


def ZooIntake2(func):
    #f1, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex='col', sharey='row')
    # N / Si / Pdt / Pc / Pdn / Pn / Zmu / Zlambda / D
    Z1 = 1
    P1 = np.linspace(0.01, 1., 50)
    P2 = np.linspace(0.01, 1., 50)

    # z.grazinprobability, 2 Z -> loops across all 4 phy
    # p.zoograzing2, 4 P -> loops across all 2 zoo
    # or
    # z.intake 2 Z -> loops across all 4 phy

    Gj = [[z2[0].fullgrazing(z2[0].zoofeeding([i, j], [Z1], func = 'vallina'), [i,j], [Z1], Z1), i, j, Z1] for i in P1 for j in P2]
    Z1Gr = np.array([z2[0].fullgrazing(z2[0].zoofeeding([i, j], [Z1], func = func), [i,j], [Z1], Z1) for i in P1 for j in P2])
    #P1Gr = np.array([p2[0].zoograzing2([z2[0].grazingprobability([i,j])], i, [1,1]) for i in P1 for j in P2])

    print(Z1Gr)

    Pone = [Gj[x][1] for x in range(50 * 50)]
    Pone2 = np.array(Pone)
    Pone3 = Pone2.reshape(50, 50)
    Ptwo = [Gj[x][2] for x in range(50 * 50)]
    Ptwo2 = np.array(Ptwo)
    Ptwo3 = Ptwo2.reshape(50, 50)

    # Figure 1
    # N
    #ax1.plot(intMLD, intPAR, c="N", lw=lws[0], alpha=alphas[0], label='Model')
    Z = Z1Gr.reshape(50, 50)
    # Legend
    #f1.align_ylabels()
    #plt.margins(x=0)
    #plt.tight_layout()

    fig4, ax4 = plt.subplots()
    cs = ax4.contourf(Pone3, Ptwo3, Z, 50)
    cbar = fig4.colorbar(cs)
    cbar.set_ticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
    ax4.set_xlabel('P1 [µM N]')
    ax4.set_ylabel('P2 [µM N]')
    ax4.set_title('Grazing of 1 Z [µM $d^{-1}$]' + func)
    plt.savefig('ZooIntake1.png')
    fig4.show()

def ZooIntake3():
    #f1, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex='col', sharey='row')
    # N / Si / Pdt / Pc / Pdn / Pn / Zmu / Zlambda / D
    Z1 = np.linspace(0.01, 1., 50)
    Z2 = np.linspace(0.01, 1., 50) # 1
    P1 = 0
    P2 = 0

    # z.grazinprobability, 2 Z -> loops across all 4 phy
    # p.zoograzing2, 4 P -> loops across all 2 zoo
    # or
    # z.intake 2 Z -> loops across all 4 phy

    #Gj = [[z3[1].fullgrazing(z3[1].zoofeeding([i, P2], [j,Z2], func = 'vallina'), [i,P2], [j,Z2], Z2), i, j, Z2] for i in P1 for j in Z1]
    #Z1Gr = np.array([z3[1].fullgrazing(z3[1].zoofeeding([i, P2], [j,Z2], func = 'anderson'), [i,P2], [j,Z2], Z2) for i in P1 for j in Z1])

    #INTERZOOGRAZE
    Gj = [[1, j, i, 1] for i in Z2 for j in Z1]

    # THIS GIVES THE GRAZING (positive) on specific zooplankton type
    Z1Gr = np.array([z3[0].interzoograze([z3[x].zoofeeding([P1, P2], [j, i], func = 'anderson') for x in range(2)], [j, i], j) for i in Z2 for j in Z1])
    #P1Gr = np.array([p2[0].zoograzing2([z2[0].grazingprobability([i,j])], i, [1,1]) for i in P1 for j in P2])


    print(Z1Gr)

    Pone = [Gj[x][1] for x in range(50 * 50)]
    Pone2 = np.array(Pone)
    Pone3 = Pone2.reshape(50, 50)
    Ptwo = [Gj[x][2] for x in range(50 * 50)]
    Ptwo2 = np.array(Ptwo)
    Ptwo3 = Ptwo2.reshape(50, 50)

    # Figure 1
    # N
    #ax1.plot(intMLD, intPAR, c="N", lw=lws[0], alpha=alphas[0], label='Model')
    Z = Z1Gr.reshape(50, 50)
    # Legend
    #f1.align_ylabels()
    #plt.margins(x=0)
    #plt.tight_layout()

    fig4, ax4 = plt.subplots()
    cs = ax4.contourf(Pone3, Ptwo3, Z, 50)
    cbar = fig4.colorbar(cs)
    #cbar.set_ticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
    ax4.set_xlabel('Z2 [µM N]')
    ax4.set_ylabel('Z1 [µM N]')
    ax4.set_title('Grazing of 1 Z [µM $d^{-1}$]')
    #plt.savefig('ZooIntake1.png')
    fig4.show()

#PIcurves()

#steelePI()

#smithPI()

#PhytoGrazed1()

#ZooIntake1()

#grazingvaryingP1()

#FASHAMgrazingvaryingP1()

#ANDERSONgrazingvaryingP1()

#grazingcomparison()
#GrazingTry2()

ZooIntake2('fasham')
ZooIntake2('anderson')
ZooIntake2('vallina')

#ZooIntake3()

#PhytoGrazed2()