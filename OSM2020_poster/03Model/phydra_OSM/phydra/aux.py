#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""These are helper functions to handle the model parameters"""

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
        try:
            return next(v for k, v in stdpardict.items() if k.endswith(parsuffix))
        except StopIteration:
            return 'Parameter not defined or not passed to class'
            # raise Exception('Parameter {} is not found in Parameters passed to SV class'.format(parsuffix))
