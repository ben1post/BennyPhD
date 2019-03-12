
from PhytoMFTM.ModelClasses import Plankton
from Tests.RunModel_2P2Zcheck import all_params

# number of phytoplankton func types
all_params.add('pfun_num', value=6, vary=False)
# number of zooplankton groups
all_params.add('zoo_num', value=3, vary=False)

pfn = all_params['pfun_num'].value
zn = all_params['zoo_num'].value

P = [2/pfn for i in range(pfn)]
Z = [1/zn for i in range(zn)]

z = Plankton(all_params, 'Zooplankton').init()
p = Plankton(all_params, 'Phytoplankton').init()


GrazeZoo = [z[i].zoograzing(P=P) for i in range(zn)]  # pass all phytoplankton
ZooItots = [x[0] for x in GrazeZoo]
RperZ = [x[1] for x in GrazeZoo]
ZooItots1 = [sum(ZooItots) for x in GrazeZoo]
RperZ1 = [sum(RperZ) for x in GrazeZoo]  ## <<- this could work, also check edit made to ZOOGRAZE formulation, can multiply by num of func types

print(sum(RperZ))
print(GrazeZoo)
print(ZooItots)
print(RperZ)

PhytoGrazed = [p[i].grazedphyto(Itot=ZooItots, P=P[i], R=RperZ1) for i in range(pfn)]  # pass all Itots + bm of the specific P


print(pfn,zn)
print("ZooItots",ZooItots, sum(ZooItots))
print("PhytoGraze",PhytoGrazed, sum(PhytoGrazed))
