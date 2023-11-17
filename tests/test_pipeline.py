from jwstmiri.pipeline import pipelinestages
from jwstmiri.pipeline.pipelinestages import checkdata

import matplotlib.pyplot as plt
import glob
from astropy.io import fits
import numpy as np

pipelinestages.pipeline('/Users/helenakuehnle/Dateien/PhD/testpipeline/',[0], [5], ['MRS'], ['LONG','MEDIUM'], [2,3], '/Users/helenakuehnle/Dateien/PhD/testpipeline/output', True,0)

#checkdata('/Users/helenakuehnle/Dateien/PhD/testpipeline/output/pipelined/stage3/cubes')
output_cubes = '/Users/helenakuehnle/Dateien/PhD/testpipeline/output/pipelined/stage3/cubes/'
output_cubes2 = '/Users/helenakuehnle/Dateien/PhD/Analysis/WISE_J0855/pipelined/stage3/cubes_nodsub/'

data_str = [d for d in glob.glob(output_cubes + "*_x1d.fits")]
data_str = sorted(data_str)
data_all = []
data_str3d = [d for d in glob.glob(output_cubes + "*_s3d.fits")]
data_str3d = sorted(data_str3d)
data_all = []
data_str2 = [d for d in glob.glob(output_cubes2 + "*short_x1d.fits")]
data_str2 = sorted(data_str2)
data_all2 = []
data_str23d = [d for d in glob.glob(output_cubes2 + "*short_s3d.fits")]
data_str23d = sorted(data_str23d)
data_all23d = []


for i in range(len(data_str)):
    data = data_str[i]
    print(data)
    data = fits.getdata(data)
    dataz = np.array(list(zip(*data)))
    data_all.append(dataz)

    data2 = data_str2[i]
    print(data2)
    data2 = fits.getdata(data2)
    dataz2 = np.array(list(zip(*data2)))
    data_all2.append(dataz2)

plt.figure(figsize=(12, 4))
for i in range(len(data_str)):
    plt.plot(data_all[i][0], data_all[i][1],color='red')
    plt.plot(data_all2[i][0], data_all2[i][1],color='black')
plt.plot(data_all[0][0], data_all[0][1],color='red',label='jwst 1.12.5')
plt.plot(data_all2[0][0], data_all2[0][1],color='black',label='jwst 1.12.2')
plt.legend()
plt.ylim([-0.0001, 0.012])
plt.xlim([3, 24])
# plt.plot(x,y,linewidth=0.5,color='k')
plt.xlabel('Wavelength [µm]')
plt.ylabel('Flux [Jy]')


plt.figure()
for i in range(len(data_str)):
    plt.plot(data_all[i][0],data_all[i][1],color='red')
plt.ylim([-0.0001, 0.012])
plt.xlim([3, 24])
# plt.plot(x,y,linewidth=0.5,color='k')
plt.xlabel('Wavelength [µm]')
plt.ylabel('Flux [Jy]')
plt.title('jwst 1.12.5')
plt.figure()
plt.xlabel('Wavelength [µm]')
plt.ylabel('Flux [Jy]')
plt.title('jwst 1.12.2')
plt.ylim([-0.0001, 0.012])
plt.xlim([3, 24])
for i in range(len(data_str)):
    plt.plot(data_all2[i][0],data_all2[i][1],color='black')
plt.show()

with fits.open(data_str[1]) as hdu:
    print(hdu.info)
    x_coord = hdu[1].header['EXTR_X']
    y_coord = hdu[1].header['EXTR_Y']

with fits.open(data_str3d[1]) as hdu3d:
    cubesum = np.nansum(hdu3d['sci'].data, axis=0)
print(data_str3d[0])
plt.figure()
plt.imshow(cubesum,origin = 'lower')#,vmax=1000,vmin=0)
plt.colorbar()
plt.plot(x_coord,y_coord,'rx')
print(y_coord)
print(x_coord)

with fits.open(data_str2[1]) as hdu:
    print(hdu.info)
    x_coord = hdu[1].header['EXTR_X']
    y_coord = hdu[1].header['EXTR_Y']

with fits.open(data_str23d[1]) as hdu3d:
    cubesum = np.nansum(hdu3d['sci'].data, axis=0)
print(data_str23d[0])
plt.figure()
plt.imshow(cubesum,origin = 'lower')#,vmax=1000,vmin=0)
plt.colorbar()
plt.plot(x_coord,y_coord,'rx')
print(y_coord)
print(x_coord)
plt.show()