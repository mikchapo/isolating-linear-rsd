import numpy as np

fake_cat = np.empty((6*80 + 1, 15))

x0 = 50.
y0 = 50.
z0 = 50.
fake_cat[0, :] = (0, -1, 5.e13, 200., 0, 500., 0, 100, x0, y0, z0, (np.random.rand()-0.5)*600., (np.random.rand()-0.5)*600., (np.random.rand()-0.5)*600., -1)

logDsep = (np.log10(100) - np.log10(0.01)) / 80.

for i in range(80):
  fake_cat[6*i+2, :] = (0, -1, 5.e13, 200., 0, 500., 0, 100, x0 + 10.**((i+0.5)*logDsep-2.), y0, z0, (np.random.rand()-0.5)*600., (np.random.rand()-0.5)*600., (np.random.rand()-0.5)*600., -1)
  fake_cat[6*i+4, :] = (0, -1, 5.e13, 200., 0, 500., 0, 100, x0, y0 + 10.**((i+0.5)*logDsep-2.), z0, (np.random.rand()-0.5)*600., (np.random.rand()-0.5)*600., (np.random.rand()-0.5)*600., -1)
  fake_cat[6*i+6, :] = (0, -1, 5.e13, 200., 0, 500., 0, 100, x0, y0, z0 + 10.**((i+0.5)*logDsep-2.), (np.random.rand()-0.5)*600., (np.random.rand()-0.5)*600., (np.random.rand()-0.5)*600., -1)

  if 10.**((i+0.5)*logDsep-2.) < 50.:
    fake_cat[6*i+1, :] = (0, -1, 5.e13, 200., 0, 500., 0, 100, x0 - 10.**((i+0.5)*logDsep-2.), y0, z0, (np.random.rand()-0.5)*600., (np.random.rand()-0.5)*600., (np.random.rand()-0.5)*600., -1)
    fake_cat[6*i+3, :] = (0, -1, 5.e13, 200., 0, 500., 0, 100, x0, y0 - 10.**((i+0.5)*logDsep-2.), z0, (np.random.rand()-0.5)*600., (np.random.rand()-0.5)*600., (np.random.rand()-0.5)*600., -1)
    fake_cat[6*i+5, :] = (0, -1, 5.e13, 200., 0, 500., 0, 100, x0, y0, z0 - 10.**((i+0.5)*logDsep-2.), (np.random.rand()-0.5)*600., (np.random.rand()-0.5)*600., (np.random.rand()-0.5)*600., -1)

  else:
    fake_cat[6*i+1, :] = (0, -1, 5.e13, 200., 0, 500., 0, 100, 1100. - 10.**((i+0.5)*logDsep-2.), y0, z0, (np.random.rand()-0.5)*600., (np.random.rand()-0.5)*600., (np.random.rand()-0.5)*600., -1)
    fake_cat[6*i+3, :] = (0, -1, 5.e13, 200., 0, 500., 0, 100, x0, 1100. - 10.**((i+0.5)*logDsep-2.), z0, (np.random.rand()-0.5)*600., (np.random.rand()-0.5)*600., (np.random.rand()-0.5)*600., -1)
    fake_cat[6*i+5, :] = (0, -1, 5.e13, 200., 0, 500., 0, 100, x0, y0, 1100. - 10.**((i+0.5)*logDsep-2.), (np.random.rand()-0.5)*600., (np.random.rand()-0.5)*600., (np.random.rand()-0.5)*600., -1)

np.savetxt("../../../projects/rrg-wperciva/mj3chapm/P2/aemulus_mocks/halos/pbc_test_cat.dat", fake_cat, header="#ID DescID M200b Vmax Vrms R200b Rs Np X Y Z VX VY VZ Parent_ID")
