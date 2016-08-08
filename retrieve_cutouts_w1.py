import core
import numpy as np

c2_info = np.genfromtxt('/Users/seb-air/Desktop/c2_query_filter.csv', delimiter=',', comments='#', dtype=str)
c2_q_field = c2_info[:, 0]
c2_q_ra = c2_info[:, 1]
c2_q_dec = c2_info[:, 2]

print c2_info

c2_cat = core.catalog()
c2_cat.LoadCatalog('/Users/seb-air/Google Drive/Astro/PhD/BCG_find/catalog/c2_clusters.csv')

out_dir = '/Users/seb-air/Desktop/c2_cutouts/'

for j in range(len(c2_cat.data['tag'])):
    xray_ra = float(c2_cat.data['ra'][j])
    xray_dec = float(c2_cat.data['dec'][j])
    cluster = c2_cat.data['tag'][j]

    dist = []
    for i in range(len(c2_q_field)):
        field = c2_q_field[i]
        field_ra = float(c2_q_ra[i])
        field_dec = float(c2_q_dec[i])

        dist.append(np.sqrt((xray_ra - field_ra)**2. + (xray_dec - field_dec)**2.))

    arg_cluster = np.argmin(dist)

    cluster_field = c2_q_field[arg_cluster]

    print cluster_field

    core.RetrieveW1Cutouts(output=out_dir + cluster + '_cutout.fits', ra=xray_ra, dec=xray_dec, radius=(5./60.), field=cluster_field, filter='R')