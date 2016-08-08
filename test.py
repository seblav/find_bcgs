import core
import numpy as np

# print core.LoadClusterCatalog()
#
# print core.GetClusterInfo('n0000')

# cluster_test = core.Cluster('n0019')
#
# #cluster_test.BuildSourcesCatalog(rad=1000.)
#
# #src = cluster_test.SourcesCatalog()
#
# #print cluster_test.SourcesCatalog()['ra']
#
# #cluster_test.BuildMembersCatalog(rad=1000., reject_stars=True, class_max=0.5, photoz_sigma=1., photoz_bias_corr=True)
# #
# # print cluster_test.MembersCatalog()['ra']
#
# z_ind = np.argmin(cluster_test.MembersCatalog()['Z'])
#
# # print z_ind, cluster_test.MembersCatalog()['Z'][z_ind], cluster_test.MembersCatalog()['ra'][z_ind], cluster_test.MembersCatalog()['dec'][z_ind]
#
# # with open('/Users/seb-air/Desktop/regions.reg', 'wb') as out:
# #     for i in range(len(cluster_test.SourcesCatalog()['ra'])):
# #         out.write('fk5; circle %7.5f %7.5f 3" # color=red\n' % (float(cluster_test.SourcesCatalog()['ra'][i]), float(cluster_test.SourcesCatalog()['dec'][i])))
# #
# #     for i in range(len(cluster_test.MembersCatalog()['ra'])):
# #         out.write('fk5; circle %7.5f %7.5f 5" # color=green width=3\n' % (float(cluster_test.MembersCatalog()['ra'][i]), float(cluster_test.MembersCatalog()['dec'][i])))
# #
# #     out.write('fk5; circle %7.5f %7.5f 8" # color=yellow width=5\n' % (float(cluster_test.MembersCatalog()['ra'][z_ind]), float(cluster_test.MembersCatalog()['dec'][z_ind])))
# #
# # out.close()
#
# # print cluster_test.contours
# #
# # cluster_test.FindBCG()
#
#
# # cat = core.LoadClusterCatalog()
# #
# # with open('/Users/seb-air/Desktop/cluster_pos.txt', 'wb') as out:
# #     for i in range(len(cat['RA'])):
# #         out.write('%s %s 5arcmin\n' % (cat['RA'][i], cat['DEC'][i]))
# #
# # out.close()
#
# # print cluster_test.BCG()['ra']
#
# cluster_test.MakeImage(show_contours=True, show_sources=True, show_members=True, show_BCG=True, show_scalebar=True)

#cluster_test = core.Cluster('n0019')

# cluster_test.BuildSourcesCatalog(rad=1000.)
#
# cluster_test.BuildMembersCatalog(rad=1000., reject_stars=True, class_max=0.5, photoz_sigma=1., photoz_bias_corr=True)

#cluster_test.FindBCG(max_rad=1000.)

#cluster_test.MakeImage(show_contours=True, show_sources=True, show_members=True, show_BCG=True, show_scalebar=True)

#
# cl_cat = core.LoadClusterCatalog()
#
# for cl, redshift in zip(cl_cat['tag'], cl_cat['redshift']):
#     if redshift > 0.:
#         cl_obj = core.Cluster(cl)
#         print cl_obj.tag
#
#         if not cl_obj.sources:
#             cl_obj.BuildSourcesCatalog()
#         if not cl_obj.members:
#             cl_obj.BuildMembersCatalog(rad=1000., reject_stars=True, class_max=0.5, photoz_sigma=1., photoz_bias_corr=True)
#         if not cl_obj.bcg:
#             cl_obj.FindBCG()
#
#         if cl_obj.cutout in ['y']:
#             cl_obj.MakeImage(show_contours=True, show_sources=True, show_members=True, show_BCG=True, show_scalebar=True)

# cl_obj = core.Cluster('n0239')
# print cl_obj.tag
# print cl_obj.sources
# print cl_obj.members
# print cl_obj.redshift
#
# if not cl_obj.sources:
#     cl_obj.BuildSourcesCatalog()
# if not cl_obj.members:
#     cl_obj.BuildMembersCatalog(rad=1000., reject_stars=True, class_max=0.5, photoz_sigma=1., photoz_bias_corr=True)
# if not cl_obj.bcg:
#     cl_obj.FindBCG()
#
# if cl_obj.cutout in ['y']:
#     cl_obj.MakeImage(show_contours=True, show_sources=True, show_members=True, show_BCG=True, show_scalebar=True)

import numpy as np
# cat_test = core.catalog()
#
# cat_test.header = np.vstack((['test', 'hello'], [12, 'hi']))
#
# cat_test.columns = ['ra', 'dec', 'yooo']

# cat_test.data = np.vstack(([35, 36, 37], [-4, -7, -8]))

# print cat_test.data.shape
#
# print cat_test.data
# cat_test.AppendCatalog(append_line=[38, -10])
#
# print cat_test.data

# # cat_test.WriteCatalog(filename='/Users/seb-air/Desktop/test_2.txt', write_columns=True, write_header=True)
# cat_test.LoadCatalog(filename='/Users/seb-air/Desktop/test_2.txt')
# print cat_test.data
# cat_test.AppendCatalog(append_column=[-100, -200], column_name='more')
# print cat_test.data
# cat_test.AppendCatalog(append_line=['a', 'b'])
#
# cat_test.WriteCatalog(filename='/Users/seb-air/Desktop/test_3.txt', write_columns=True, write_header=True)


# cat = np.genfromtxt('/Users/seb-air/Desktop/test_3.txt', dtype=None, delimiter=',', names=True, comments='#', skip_header=2)
#
# print cat.

# print cat_test.data
#
# cat_test.AppendCatalog(append_line=[38, -10])


# print cat_test.header.shape
# print cat_test.data
#
# cat_test.WriteCatalog(filename='/Users/seb-air/Desktop/test.txt', write_columns=True, write_header=True)
#
# cat_test.LoadCatalog(filename='/Users/seb-air/Desktop/test.txt')
#
# #print cat_test.header.keys()
# #print cat_test.data[cat_test.data.keys()[0]]
# #print cat_test.columns
#
# # cat_test.WriteCatalog(filename='/Users/seb-air/Desktop/test_2.txt', write_columns=True, write_header=True)
# print len(cat_test.data.keys())
#
# cat_test.AppendCatalog(append_line=[38, -10])
# cat_test.AppendCatalog(append_column=[11,22,33,44], column_name='test')
#
# print cat_test.data
#
# print cat_test.header['test']
# print cat_test.data['ra']
#
# cat_test.WriteCatalog(filename='/Users/seb-air/Desktop/test_3.txt', write_columns=True, write_header=True)

#
# cat_test.WriteCatalog(filename='/Users/seb-air/Desktop/test_2.txt', write_columns=True, write_header=True)

# cl_obj = core.Cluster('n0019')

#cl_obj.BuildSourcesCatalog(rad=1000., overwrite=True)

# src = cl_obj.SourcesCatalog()
#
# print src.data['ra']
# print cl_obj.SourcesCatalog().data['ra']

#cl_obj.BuildMembersCatalog(rad=1000., reject_stars=True, class_max=0.5, photoz_sigma=1., photoz_bias_corr=True)

# core.RetrieveW1Cutouts(1,1,1,1,'W1+2+3','R')
#
# c2_info = np.genfromtxt('/Users/seb-air/Desktop/c2_query_filter.csv', delimiter=',', comments='#', dtype=str)
# c2_q_field = c2_info[:, 0]
# c2_q_ra = c2_info[:, 1]
# c2_q_dec = c2_info[:, 2]
#
# print c2_info
#
# c2_cat = core.catalog()
# c2_cat.LoadCatalog('/Users/seb-air/Google Drive/Astro/PhD/BCG_find/catalog/c2_clusters.csv')
#
# out_dir = '/Users/seb-air/Desktop/c2_cutouts/'
#
# for j in range(len(c2_cat.data['tag'])):
#     xray_ra = float(c2_cat.data['ra'][j])
#     xray_dec = float(c2_cat.data['dec'][j])
#     cluster = c2_cat.data['tag'][j]
#
#     dist = []
#     for i in range(len(c2_q_field)):
#         field = c2_q_field[i]
#         field_ra = float(c2_q_ra[i])
#         field_dec = float(c2_q_dec[i])
#
#         dist.append(np.sqrt((xray_ra - field_ra)**2. + (xray_dec - field_dec)**2.))
#
#     arg_cluster = np.argmin(dist)
#
#     cluster_field = c2_q_field[arg_cluster]
#
#     print cluster_field
#
#     core.RetrieveW1Cutouts(output=out_dir + cluster + '_cutout.fits', ra=xray_ra, dec=xray_dec, radius=(5./60.), field=cluster_field, filter='R')




# cl_cat = core.catalog()
# cl_cat.LoadCatalog('/Users/seb-air/Google Drive/Astro/PhD/BCG_find/catalog/c2_clusters.csv')
#
# for cl, redshift in zip(cl_cat.data['tag'], cl_cat.data['redshift']):
#     if redshift > 0.:
#         cl_obj = core.Cluster(cl)
#         print cl_obj.tag
#
#         out = '/Users/seb-air/Desktop/c2_images/' + cl_obj.tag + '.pdf'
#
#         if cl_obj.cutout in ['y'] and cl_obj.status in ['confirmed']:
#             cl_obj.MakeImage(out_path=out, show_contours=True, show_sources=False, show_members=False, show_BCG=False, show_scalebar=True)

# cl_obj = core.Cluster('n0019')
# # cl_obj.BuildMembersCatalog(reject_stars=True, photoz_bias_corr=True, overwrite=True)
#
# # cat = cl_obj.MembersCatalog()
# #print cat.data['ra']
#
# cl_obj.FindBCG(overwrite=True)
#
# print cl_obj.BCGInfo().data['Z']

# cl_cat = core.catalog()
# cl_cat.LoadCatalog('/Users/seb-air/Google Drive/Astro/PhD/BCG_find/catalog/c2_clusters.csv')
#
# for cl, redshift in zip(cl_cat.data['tag'][50:], cl_cat.data['redshift'][50:]):
#     if 0. < float(redshift) < 1.2:
#         cl_obj = core.Cluster(cl)
#         print cl_obj.tag
#
#         if cl_obj.dec < -4.:
#             try:
#                 cl_obj.BuildSourcesCatalog(overwrite=True)
#                 cl_obj.BuildMembersCatalog(overwrite=True, reject_stars=True)
#                 cl_obj.FindBCG(overwrite=True)
#
#                 out = '/Users/seb-air/Google Drive/Astro/PhD/BCG_find/image/' + cl_obj.tag + '.pdf'
#
#                 if cl_obj.cutout in ['y'] and cl_obj.status in ['confirmed'] and float(cl_obj.BCGInfo().data['Z']) > 0.:
#                     cl_obj.MakeImage(out_path=out, show_contours=True, show_sources=True, show_members=True, show_BCG=True, show_scalebar=True)
#             except:
#                 print 'Error while building catalog for ' + cl_obj.tag
#                 pass

#
# c1_info = np.genfromtxt('/Users/seb-air/Google Drive/Astro/PhD/BCG_find/c1_query_filter.csv', delimiter=',', comments='#', dtype=str)
# c1_q_field = c1_info[:, 0]
# c1_q_ra = c1_info[:, 1]
# c1_q_dec = c1_info[:, 2]
#
# print c1_info
#
# c1_cat = core.catalog()
# c1_cat.LoadCatalog('/Users/seb-air/Google Drive/Astro/PhD/BCG_find/catalog/c1_clusters.csv')
#
# out_dir = '/Users/seb-air/Desktop/c1_cutouts/'
#
# for j in range(len(c1_cat.data['tag'])):
#     xray_ra = float(c1_cat.data['ra'][j])
#     xray_dec = float(c1_cat.data['dec'][j])
#     cluster = c1_cat.data['tag'][j]
#
#     dist = []
#     for i in range(len(c1_q_field)):
#         field = c1_q_field[i]
#         field_ra = float(c1_q_ra[i])
#         field_dec = float(c1_q_dec[i])
#
#         dist.append(np.sqrt((xray_ra - field_ra)**2. + (xray_dec - field_dec)**2.))
#
#     arg_cluster = np.argmin(dist)
#
#     cluster_field = c1_q_field[arg_cluster]
#
#     print cluster_field
#
#     core.RetrieveW1Cutouts(output=out_dir + cluster + '_cutout.fits', ra=xray_ra, dec=xray_dec, radius=(5./60.), field=cluster_field, filter='R')


import core
import numpy as np

cl_cat = core.catalog()
cl_cat.LoadCatalog('/Users/seb-air/Google Drive/Astro/PhD/BCG_find/catalog/c1_clusters.csv')

for cl, redshift in zip(cl_cat.data['tag'], cl_cat.data['redshift']):
    if 0. < float(redshift) < 1.2:
        cl_obj = core.Cluster(cl)
        print cl_obj.tag

        if cl_obj.dec < -4.:
            cl_obj.FindBCG(overwrite=True, max_rad=350.)
            cl_obj.FindBCG(overwrite=True, max_rad=350., xbrightest=3)

            out = '/Users/seb-air/Google Drive/Astro/PhD/BCG_find/image/c1/' + cl_obj.tag + '.pdf'

            if cl_obj.cutout in ['y'] and cl_obj.status in ['confirmed'] and float(cl_obj.BCGInfo().data['Z']) > 0.:
                cl_obj.MakeImage(out_path=out, show_contours=True, show_sources=True, show_members=True, show_BCG=True, show_scalebar=True, show_xBCGs=True)
