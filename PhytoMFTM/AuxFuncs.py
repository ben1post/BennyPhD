#!/usr/bin/env python
# -*- coding: utf-8 -*-


# For ParDict:

def sliceparams(pardict, parprefix):
    """Function to extract functioanl type parameters from Parameters object,
    using prefix in key"""
    return {k: v.value for k, v in pardict.items() if k.startswith(parprefix)}


def sliceoffparams(pardict, parprefix):
    """Function to remove e.g. functioanl type parameters from Parameters object,
    using prefix in key"""
    return {k: v.value for k, v in pardict.items() if not k.startswith(parprefix)}


def extractparam(pardict, parsuffix):
    """Function to extract certain single parameter from sliced (!) Parameters object,
    using final characters of key"""
    return next(v for k, v in pardict.items() if k.endswith(parsuffix))


def checkreplaceparam(stdpardict, functypepardict, parsuffix):
    """Function to check that certain single parameter from sliced (!) Parameters object,
    using final characters of key"""
    try:
        ftpara = next(v for k, v in functypepardict.items() if k.endswith(parsuffix))
        return ftpara
    except StopIteration:
        return next(v for k, v in stdpardict.items() if k.endswith(parsuffix))



# For Grazing:

# FUNCTIONS to handle multiple type grazing
# (perhaps implement this within classes, for now outside)

def feedingmatrix(P, Z, pn, zn):
    l = [[phy, zoo] for phy in range(pn) for zoo in range(zn)]
    return l


def grazing(P, Z, pn, zn, zclass):
    l = [[phy, zoo] for phy in range(pn) for zoo in range(zn)]
    all_graze = [[phy, zoo, zclass[zoo].grazing(P[phy], Z[zoo])] for phy, zoo in l]
    return all_graze


def zoogrowth_zoo(GrzingMat, P, pn, zn):
    pairwise = [graz * P[phy] for phy, zoo, graz in GrzingMat]
    print(pairwise)
    sumperzootype = [sum([pairwise[zoo + (phy * zn)] for phy in range(pn)]) for zoo in range(zn)]
    return sumperzootype  # values in list


def zoograzeloss_phy(GrzingMat, pn, zn):
    pairwise = [graz for phy, zoo, graz in GrzingMat]
    sumperphytype = [sum([pairwise[zoo + (phy * zn)] for zoo in range(zn)]) for phy in range(pn)]
    return sumperphytype  # values in list

