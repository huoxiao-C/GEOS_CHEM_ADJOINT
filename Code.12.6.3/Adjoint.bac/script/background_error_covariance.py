

import xarray as xr

import numpy as np





def Barometri_formula(pressure):

 '''



 :param pressure:

 :return:

 '''



 p0 = 101325

 cp = 1004.68506

 g = 9.80665

 M = 0.02896968

 R0 = 8.314462618

 L = 0.00976

 T0 = 288.16



 h = (1 - np.exp(np.log(pressure / p0) * R0 / (cp * M))) * (cp * T0) / g

 return h



def exp_dis(dis, cons, index):

    temp = np.exp(-np.power(dis, index)/np.power(cons, index))

    return temp


root_dir = '/share/nas1_share1/huoxiao_share/GEOSChem_Simulation_test/'
gc_nc_data_24 = xr.open_dataset(root_dir+'GEOSChem.SpeciesConc.20171002_0000z.nc4')

met_nc_data_24 = xr.open_dataset(root_dir+'GEOSChem.StateMet.20171002_0000z.nc4')

gc_nc_data_48 = xr.open_dataset(root_dir+'GEOSChem.SpeciesConc.20171003_0000z.nc4')

met_nc_data_48 = xr.open_dataset(root_dir+'GEOSChem.StateMet.20171003_0000z.nc4')



# vert = []

# height = np.zeros(47 )

# pressure = []

# corre = np.zeros((47,47))

# diff_list = []

# diff_layer = []

# for j in range(47):

#     height[j] = \

#         Barometri_formula(np.mean(met_nc_data_48.variables['Met_PMID'][0, j, :, :] * 100))

#     for i in range(j, 47):

#

#         layer_ith = np.mean(gc_nc_data_48.variables['SpeciesConc_CO2'][0, j, :, :])-\

#                     np.mean(gc_nc_data_24.variables['SpeciesConc_CO2'][0, j, :, :])

#         layer_ithp = np.mean(gc_nc_data_48.variables['SpeciesConc_CO2'][0, i, :, :])-\

#                     np.mean(gc_nc_data_24.variables['SpeciesConc_CO2'][0, i, :, :])

#         diff = abs(layer_ith*layer_ithp/(layer_ith*layer_ith))

#         if i == j:

#             min = layer_ith*layer_ithp/(layer_ith*layer_ith)

#

#         else:

#             if diff<min:

#                 min = diff

#         if diff >= min:

#             diff = min

#         # if i != 46:

#         #     diff_layer.append(np.mean(gc_nc_data_48.variables['SpeciesConc_CO2'][0, i, ...])*10**6)

#         corre[j, i] = diff

#         #

#         # pressure.append(np.mean(met_nc_data_48.variables['Met_PMID'][0, i, :, :]*100))

#         # if i == 0:

#         #    corre.append(1)

#         # else:

#         #

#         #    corre.append(exp_dis((height[-1]-height[-2]), 150, 1.2))

#         #    if i == 1:

#         #        print('test:', exp_dis((height[-1]-height[-2]), 150, 1.2))

#            # corre.append(exp_dis((abs(i-50)*200), 1000, 1))

#            # print(exp_dis((abs(i-50)*200), 200, 1))

#

# gaus_scatter = np.zeros(47*24)

# nms_scatter = np.zeros(47*24)

# for i in range(47):

#     for j in range(0, i):

#         corre[i, j] = corre[j, i]

# eval, evect = np.linalg.eig(corre)

# diag = np.diag(abs(eval))

# corre = np.dot(np.dot(evect, diag),np.linalg.inv(evect))

# gaus_corre = np.zeros((47, 47))

# for i in range(47):

#     for j in range(47):

#         gaus_corre[i, j] = exp_dis(abs(height[j]-height[i]), 2000, 2.2)

# index = 0

# print(height)

# # for i in range(47):

# #     for j in range(0, i+1):

# #         gaus_scatter[index] =

#

# plt.pcolor(corre)

# plt.show()

# plt.pcolor(gaus_corre)

# plt.show()

#



lat = gc_nc_data_24.variables['lat'][:]

lon = gc_nc_data_24.variables['lon'][:]



L = 1

Dist = 1

lat_num = 91

lon_num = 144

lev = 47



Sx = np.zeros((lon_num, lon_num))

Sy = np.zeros((lat_num, lat_num))

Sz = np.zeros((lev, lev))

#factorization

#x dimension

for row in range(lon_num):

    for col in range(lon_num):

        Sx[row, col] = np.exp(-abs((row-col)*Dist)**2/(2*L**2))

Sx_h = np.linalg.cholesky(Sx)

#y dimension

for row in range(lat_num):

    for col in range(lat_num):

        Sy[row, col] = np.exp(-abs((row-col)*Dist)**2/(2*L**2))

Sy_h = np.linalg.cholesky(Sy)

#z dimension

height = np.zeros(47)

for j in range(lev):

    height[j] = \
         Barometri_formula(np.mean(met_nc_data_48.variables['Met_PMID'][0, j, :, :] * 100))

for row in range(lev):

    for col in range(lev):

        Sz[row, col] = exp_dis(abs(height[row]-height[col]), 2000, 1.4)

Sz_h = np.linalg.cholesky(Sz)







D = np.sqrt(np.abs(gc_nc_data_48.variables['SpeciesConc_CO2'][...]-gc_nc_data_24.variables['SpeciesConc_CO2'][...]))

xr_out = xr.Dataset(data_vars={'Sx_h': (('lon', 'lon'), Sx_h),

                               'Sy_h': (('lat', 'lat'), Sy_h),

                                'Sz_h': (('lev', 'lev'), Sz_h),

                               'D': (('time','lev', 'lat', 'lon'), D)},

                    coords={'time': (np.arange(1)), 'lev':(np.arange(47)),'lat':(lat), 'lon':(lon)})
xr_out.to_netcdf('./BEC.nc4')
