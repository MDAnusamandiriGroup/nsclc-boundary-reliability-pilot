#!/usr/bin/env python3
"""Construct two observationally equivalent Poisson count models.
A global biological amplification and a technical capture multiplier have
the same expected counts. No count-only analysis can tell which generated
those observations without more assumptions or independent measurements.
"""
import numpy as np
rng=np.random.default_rng(20261023)
n,g=500,80
region=np.repeat([0,1],n//2)
mu=rng.lognormal(1.5,0.65,g)
factor=np.where(region==1,2.0,1.0)
lam_tech=factor[:,None]*mu[None,:]            # capture efficiency doubles
lam_biol=factor[:,None]*mu[None,:]            # biological global RNA doubles
assert np.array_equal(lam_tech,lam_biol)
Y=rng.poisson(lam_tech)
print("Equal conditional Poisson intensities:",np.array_equal(lam_tech,lam_biol))
print("Median total counts by spatial region:",
      np.median(Y.sum(axis=1)[region==0]),
      np.median(Y.sum(axis=1)[region==1]))
print("Interpretation: biological versus technical global scaling is not identifiable from these observed counts alone.")
