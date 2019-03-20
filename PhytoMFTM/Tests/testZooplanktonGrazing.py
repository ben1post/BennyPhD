
from PhytoMFTM.ModelClasses import Plankton
from Tests.RunModel_2P2Zcheck import all_params

def rungrazecalc(pfn,zn):
    # number of phytoplankton func types
    all_params.add('pfun_num', value=pfn, vary=False)
    # number of zooplankton groups
    all_params.add('zoo_num', value=zn, vary=False)

    pfn = all_params['pfun_num'].value
    zn = all_params['zoo_num'].value

    P = [1+i for i in range(pfn)]
    Z = [2/zn for i in range(zn)]

    z = Plankton(all_params, 'Zooplankton').init()
    p = Plankton(all_params, 'Phytoplankton').init()


    GrazeZoo = [z[i].zoograzing(P=P, Z=Z[i]) for i in range(zn)]  # pass all phytoplankton

    ZooItots = [x[0] for x in GrazeZoo]
    RperZ = [x[1] for x in GrazeZoo]

    ZooItots1 = [sum(ZooItots) for x in GrazeZoo]
    RperZ1 = [sum(RperZ) for x in GrazeZoo]  ## <<- this could work, also check edit made to ZOOGRAZE formulation, can multiply by num of func types



    PhytoGrazed = [p[i].grazedphyto(Itot=ZooItots, P=P[i], R=RperZ) for i in range(pfn)]  # pass all Itots + bm of the specific P


    print(pfn,zn)
    print("ZooItots",ZooItots, sum(ZooItots))
    print("PhytoGraze",PhytoGrazed, sum(PhytoGrazed))


rungrazecalc(1,1)
rungrazecalc(2,1)
rungrazecalc(2,2)
rungrazecalc(2,3)
rungrazecalc(3,3)
rungrazecalc(30,30)