import xarray as xr
import numpy as np
import sys
import datetime
from scipy import linalg as lg

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
    temp = np.exp(-np.power(dis, index) / np.power(cons, index))

    return temp


def main(year, month, day, hour):
    time_str = year + month + day + '_' + hour + '00'
    print('Background error covariance calculation in ' + time_str)
    begin_time_date = datetime.datetime.strptime(time_str, '%Y%m%d_%H%M')
    end_time_date = begin_time_date + datetime.timedelta(days=1)

    root_dir = '/share/nas1_share1/huoxiao_share/GEOSChem_Simulation_test/'
    gc_nc_data_24 = xr.open_dataset(root_dir + 'GEOSChem.SpeciesConc.' +
                                    '{:0>4d}{:0>2d}{:0>2d}_{:0>2d}00z.nc4'.
                                    format(begin_time_date.year, begin_time_date.month,
                                           begin_time_date.day, begin_time_date.hour))

    met_nc_data_24 = xr.open_dataset(root_dir + 'GEOSChem.StateMet.' +
                                     '{:0>4d}{:0>2d}{:0>2d}_{:0>2d}00z.nc4'.
                                     format(begin_time_date.year, begin_time_date.month,
                                            begin_time_date.day, begin_time_date.hour))

    gc_nc_data_48 = xr.open_dataset(root_dir + 'GEOSChem.SpeciesConc.' +
                                    '{:0>4d}{:0>2d}{:0>2d}_{:0>2d}00z.nc4'.
                                    format(end_time_date.year, end_time_date.month,
                                           end_time_date.day, end_time_date.hour))

    met_nc_data_48 = xr.open_dataset(root_dir + 'GEOSChem.StateMet.' +
                                     '{:0>4d}{:0>2d}{:0>2d}_{:0>2d}00z.nc4'.
                                     format(end_time_date.year, end_time_date.month,
                                            end_time_date.day, end_time_date.hour))

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

    Lx = 5
    Ly = 2

    Dist = 1

    lat_num = 91

    lon_num = 144

    lev = 47

    Sx = np.zeros((lon_num, lon_num))

    Sy = np.zeros((lat_num, lat_num))

    Sz = np.zeros((lev, lev))

    # factorization

    # x dimension

    for row in range(lon_num):

        for col in range(lon_num):
            Sx[row, col] = np.exp(-abs((row - col) * Dist)  / (Lx ))

    Sx_h = np.linalg.cholesky(Sx)
    Sx_inv = np.linalg.inv(Sx)
    # y dimension

    for row in range(lat_num):

        for col in range(lat_num):
            Sy[row, col] = np.exp(-abs((row - col) * Dist) / (Ly ))

    Sy_h = np.linalg.cholesky(Sy)
    Sy_inv = np.linalg.inv(Sy)
    # z dimension

    height = np.array([72.180, 63.053, 54.834, 47.135, 40.166, 34.024, 28.654, 25.307,
                       23.020, 20.836, 18.727, 17.243, 16.222, 15.198, 14.170, 13.134,
                       12.086, 11.021, 9.936, 8.846, 7.943, 7.237, 6.585, 5.980, 5.413,
                       4.879, 4.375, 3.896, 3.439, 3.074, 2.792, 2.517, 2.249, 1.988,
                       1.759, 1.584, 1.436, 1.290, 1.146, 1.004, 0.864, 0.726, 0.589,
                       0.454, 0.320, 0.189, 0.058]) * 1000
    gc_pressure = np.array([1005.650, 990.408, 975.122, 959.837, 944.553, 929.268, 913.984, 898.701,
                            883.418, 868.135, 852.852, 837.570, 819.743, 796.822, 771.354, 745.890, 720.429,
                            694.969, 663.146, 624.967, 586.793, 548.628, 510.475, 472.335, 434.212, 396.112,
                            358.038, 313.966, 267.087, 226.745, 192.587, 163.661, 139.115, 118.250, 100.514,
                            85.439, 67.450, 48.282, 34.272, 24.080, 14.542, 6.685, 2.864, 1.134, 0.414, 0.139,
                            0.038
                            ]) * 100
    height.sort()
    for row in range(lev):
        for col in range(lev):
            Sz[row, col] = exp_dis(abs(gc_pressure[row] - gc_pressure[col]), 12000, 1)

    Sz_h = np.linalg.cholesky(Sz)

    Sz_inv = np.linalg.inv(Sz)
    print(np.max(Sx_inv), np.max(Sy_inv), np.max(Sz_inv), np.min(Sx_inv), np.min(Sy_inv), np.min(Sz_inv))
#    Sx_inv /= 1.0 * 10 ** 6
#    Sy_inv /= 100.0
#    Sz_inv /= 100.0
    print(np.max(Sx_inv), np.max(Sy_inv), np.max(Sz_inv), np.min(Sx_inv), np.min(Sy_inv), np.min(Sz_inv))
#    print(Sx_inv[58, 49])
#    Sx_inv = lg.inv(Sx)
    print('Sx', np.max(np.dot(Sx, Sx_inv)), np.max(np.dot(Sx_inv, Sx)))
    print('Sy', np.max(np.dot(Sy, Sy_inv)), np.max(np.dot(Sy_inv, Sy)))
    print('Sz', np.max(np.dot(Sz, Sz_inv)), np.max(np.dot(Sz_inv, Sz)))
    D = np.sqrt(
        np.abs(gc_nc_data_48.variables['SpeciesConc_CO2'][...] - gc_nc_data_24.variables['SpeciesConc_CO2'][...]))

    D_inv = np.square(gc_nc_data_48.variables['SpeciesConc_CO2'][...] - gc_nc_data_24.variables['SpeciesConc_CO2'][...])
    xr_out = xr.Dataset(data_vars={'Sx_h': (('lon', 'lon'), Sx),
                                   'Sy_h': (('lat', 'lat'), Sy),
                                   'Sz_h': (('lev', 'lev'), Sz),
                                   'Sx_inv': (('lon', 'lon'), Sx_inv),
                                   'Sy_inv': (('lat', 'lat'), Sy_inv),
                                   'Sz_inv': (('lev', 'lev'), Sz_inv),
                                   'D': (('time', 'lev', 'lat', 'lon'), D),
                                   'D_inv': (('time', 'lev', 'lat', 'lon'), D)},

                        coords={'time': (np.arange(1)), 'lev': (np.arange(47)), 'lat': (lat), 'lon': (lon)})
    xr_out.to_netcdf('~/GEOS_CHEM_ADJOINT_V1_INC/rundirs/merra2_2x25_CO2/BEC.nc4')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
