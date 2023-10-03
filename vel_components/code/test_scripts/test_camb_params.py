import camb
import numpy as np

omch2_planck = 0.1199
ombh2_planck = 0.02222
omh2_planck = 0.14212
h_planck = 0.6726
ns_planck = 0.9652
sigma8_planck = 0.830
As_default = 2e-09
sigma8_default = 0.7917199
As_rescale = (sigma8_planck / sigma8_default) ** 2.

zs = [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1.5, 5., 20., 49.]
zs.reverse()

pars = camb.CAMBparams()
pars.set_cosmology(H0=h_planck*100., ombh2=ombh2_planck, omch2=omch2_planck)
pars.InitPower.set_params(ns=ns_planck, As=As_default*As_rescale)
pars.set_matter_power(redshifts=zs, kmax=2.0)

pars.NonLinear = camb.model.NonLinear_none
results = camb.get_results(pars)
sigma8s = np.array(results.get_sigma8())
fsigma8s = np.array(results.get_fsigma8())

for i, z in enumerate(zs):
	H = results.hubble_parameter(z)
	print("z={}".format(z), "sigma8={}".format(sigma8s[i]), "fsigma8={}".format(fsigma8s[i]), "H={}".format(H))