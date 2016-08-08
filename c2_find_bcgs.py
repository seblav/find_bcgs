import core
import numpy as np

cl_cat = core.catalog()
cl_cat.LoadCatalog('/Users/seb-air/Google Drive/Astro/PhD/BCG_find/catalog/c2_clusters.csv')


for cl, redshift in zip(cl_cat.data['tag'], cl_cat.data['redshift']):
    if 0. < float(redshift) < 1.2:
        cl_obj = core.Cluster(cl)
        print cl_obj.tag

        if cl_obj.dec < -4.:
            try:
                cl_obj.FindBCG(overwrite=True, max_rad=350.)
                cl_obj.FindBCG(overwrite=True, max_rad=350., xbrightest=3)

                out = '/Users/seb-air/Google Drive/Astro/PhD/BCG_find/image/c2/' + cl_obj.tag + '.pdf'

                if cl_obj.cutout in ['y'] and cl_obj.status in ['confirmed'] and float(cl_obj.BCGInfo().data['Z']) > 0.:
                    cl_obj.MakeImage(out_path=out, show_contours=True, show_sources=True, show_members=True, show_BCG=True,
                                     show_scalebar=True, show_xBCGs=True)
            except:
                print 'Error while building catalog for ' + cl_obj.tag
                pass