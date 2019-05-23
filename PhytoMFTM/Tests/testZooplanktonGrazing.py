from PhytoMFTM.ModelClasses import Plankton
from Tests.runModel_CARIACO_Prelim import all_params

def rungrazecalc(pfn,zn):
    # number of phytoplankton func types
    all_params.add('pfun_num', value=pfn, vary=False)
    # number of zooplankton groups
    all_params.add('zoo_num', value=zn, vary=False)

    pfn = all_params['pfun_num'].value
    zn = all_params['zoo_num'].value

    P = [2 for i in range(pfn)]
    Z = [2/zn for i in range(zn)]
    N = 1
    Si = 1

    z = Plankton(all_params, 'Zooplankton').init()
    p = Plankton(all_params, 'Phytoplankton').init()

    # Phytoplankton related processes
    # Nutrient uptake
    N_Uptake = [p[i].n_uptake(N) for i in range(pfn)]
    Si_Uptake = [p[i].si_uptake(Si) for i in range(pfn)]

    # Light and Temperature
    LightHarvesting = [1 for i in range(pfn)]
    TemperatureDepGrowth = [1 for i in range(pfn)]
    # Phytoplankton Growth
    Gains = [p[i].gains(N_Uptake[i], Si_Uptake[i], LightHarvesting[i], TemperatureDepGrowth[i], P[i]) for i in range(pfn)]
    SilicateDrawdown = [p[i].silicatedrawdown(Gains[i]) for i in range(pfn)]

    # Phytoplankton related processes
    # Nutrient uptake
    N_Uptake = [p[i].n_uptake(N) for i in range(pfn)]
    Si_Uptake = [p[i].si_uptake(Si) for i in range(pfn)]

    # Zooplankton Grazing:
    Rj = [z[i].ressourcedensity(P) for i in range(zn)] # total available ressource density for Zooplankton 1
    #Rj = [sum(Rji) for i in range(zn)]
    Itot = [z[i].itot(Rj[i]) for i in range(zn)]

    PhytoGrazed = [p[i].zoograzing(Itot, Rj, P[i], Z) for i in range(pfn)] #returns phyto grazed per type
    # p1 [grazed by z1, z2], p2 [grazed by z1,z2]

    AssimilatedGrazing = [z[i].assimgrazing(Itot[i], Z[i]) for i in range(zn)]
    UnassimilatedGrazing = [z[i].unassimilatedgrazing(Itot[i], Z[i]) for i in range(zn)]


    ZooMixing = [Z[i] * 0 for i in range(zn)]
    ZooMortality = [z[i].zoomortality(Z[i]) for i in range(zn)]


    # Phytoplankton losses
    PhytoMortality = [p[i].mortality(P[i]) for i in range(pfn)]
    PhytoSinking = [p[i].sinking(20, P[i]) for i in range(pfn)]
    PhytoMixing = [P[i] * 0 for i in range(pfn)]

    zoo = [AssimilatedGrazing[i] - ZooMortality[i] - ZooMixing[i] for i in range(zn)]    # Zooplankton losses due to mortality and mixing
    phy = [Gains[i] - PhytoGrazed[i] - PhytoMortality[i] - PhytoMixing[i] - PhytoSinking[i] for i in range(pfn)]  # Phytoplankton growth
    print('Zoo',sum(zoo),'Phy',sum(phy),'AssimGraz',sum(AssimilatedGrazing), 'ZooMortality',sum(ZooMortality),'ZooMixing',sum(ZooMixing))



    P = [P[i]+phy[i] for i in range(pfn)]
    Z = [Z[i]+zoo[i] for i in range(zn)]
    print(pfn,zn)
    print('Rj',Rj,'Itot',Itot)
    print("AssimilatedGrazing", AssimilatedGrazing, sum(AssimilatedGrazing))
    print("PhytoGraze", PhytoGrazed, sum(PhytoGrazed))
    print('P',sum(P),'grz',sum(PhytoGrazed),'Z',sum(Z),'ZG,UP',sum(AssimilatedGrazing+UnassimilatedGrazing))

    #print("PhytoGrazing", PhytoGrazing, sum(PhytoGrazing))

def rungrazecalcNEW(pfn,zn):
    # number of phytoplankton func types
    all_params.add('pfun_num', value=pfn, vary=False)
    # number of zooplankton groups
    all_params.add('zoo_num', value=zn, vary=False)

    pfn = all_params['pfun_num'].value
    zn = all_params['zoo_num'].value

    P = [1,2,3,4] #[1/4 for i in range(pfn)]
    Z = [1,1] #[1/2 for i in range(zn)]
    N = 1
    Si = 1

    z = Plankton(all_params, 'Zooplankton').init()
    p = Plankton(all_params, 'Phytoplankton').init()

    # Phytoplankton related processes
    # Nutrient uptake
    N_Uptake = [p[i].n_uptake(N) for i in range(pfn)]
    Si_Uptake = [p[i].si_uptake(Si) for i in range(pfn)]

    # Light and Temperature
    LightHarvesting = [1 for i in range(pfn)]
    TemperatureDepGrowth = [1 for i in range(pfn)]
    # Phytoplankton Growth
    Gains = [p[i].gains(N_Uptake[i], Si_Uptake[i], LightHarvesting[i], TemperatureDepGrowth[i], P[i]) for i in range(pfn)]
    SilicateDrawdown = [p[i].silicatedrawdown(Gains[i]) for i in range(pfn)]

    # Phytoplankton related processes
    # Nutrient uptake
    N_Uptake = [p[i].n_uptake(N) for i in range(pfn)]
    Si_Uptake = [p[i].si_uptake(Si) for i in range(pfn)]

    # Zooplankton Grazing:
    Gj = [z[i].grazingprobability(P) for i in range(zn)] # total available ressource density for Zooplankton 1
    Itot = [z[i].zoointake(Gj[i], P) for i in range(zn)]

    PhytoGrazed = [p[i].zoograzing2(Gj, P[i], Z) for i in range(pfn)] #returns phyto grazed per type
    # p1 [grazed by z1, z2], p2 [grazed by z1,z2]

    AssimilatedGrazing = [z[i].assimgrazing(Itot[i], Z[i]) for i in range(zn)]
    UnassimilatedGrazing = [z[i].unassimilatedgrazing(Itot[i], Z[i]) for i in range(zn)]


    ZooMixing = [Z[i] * 0 for i in range(zn)]
    ZooMortality = [z[i].zoomortality(Z[i]) for i in range(zn)]


    # Phytoplankton losses
    PhytoMortality = [p[i].mortality(P[i]) for i in range(pfn)]
    PhytoSinking = [p[i].sinking(20, P[i]) for i in range(pfn)]
    PhytoMixing = [P[i] * 0 for i in range(pfn)]

    zoo = [AssimilatedGrazing[i] - ZooMortality[i] - ZooMixing[i] for i in range(zn)]    # Zooplankton losses due to mortality and mixing
    phy = [Gains[i] - PhytoGrazed[i] - PhytoMortality[i] - PhytoMixing[i] - PhytoSinking[i] for i in range(pfn)]  # Phytoplankton growth
    print('Zoo',sum(zoo),'Phy',sum(phy),'AssimGraz',sum(AssimilatedGrazing), 'ZooMortality',sum(ZooMortality),'ZooMixing',sum(ZooMixing))



    P = [P[i]+phy[i] for i in range(pfn)]
    Z = [Z[i]+zoo[i] for i in range(zn)]
    print(pfn,zn)
    print('Gj',Gj,'Itot',Itot)
    print("AssimilatedGrazing", AssimilatedGrazing, sum(AssimilatedGrazing))
    print("PhytoGraze", PhytoGrazed, sum(PhytoGrazed))
    print('P',sum(P),'grz',sum(PhytoGrazed),'Z',sum(Z),'ZG,UP',sum(AssimilatedGrazing+UnassimilatedGrazing))

    #print("PhytoGrazing", PhytoGrazing, sum(PhytoGrazing))


rungrazecalcNEW(4,2)
#rungrazecalcNEW(1,2)
#rungrazecalcNEW(2,1)
#rungrazecalcNEW(2,2)
#rungrazecalcNEW(4,1)
#rungrazecalcNEW(4,2)