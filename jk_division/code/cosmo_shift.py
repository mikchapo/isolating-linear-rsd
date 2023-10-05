from astropy.io import ascii as ap_ascii
from astropy.table import Table
import numpy as np

caps = ["NGC", "SGC"]
cats = ["dat", "rand"]
cols = ["X", "Y", "Z"]

for cap in caps:
    print("Starting %s" % cap)
    for cat in cats:
        print("Starting %s" % cat)
        mc = Table.read("../archive/mc/output/eBOSS_LRG_pip_v7_2_%s_jk_200.%s" % (cap, cat), format="ascii")
        ec = Table.read("../data/eBOSS_LRG_%s_pip_v7_2.%s" % (cap, cat), format="ascii")
        for col in cols:
            print("Starting %s" % col)
            mc[col][:] = ec[col][:]
        ap_ascii.write(mc, output="../output/eBOSS_LRG_pip_v7_2_%s_jk_200.%s" % (cap, cat), delimiter="\t", format="commented_header")
