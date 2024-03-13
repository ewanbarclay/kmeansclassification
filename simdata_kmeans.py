#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from scipy.spatial import ConvexHull
import pandas as pd


import glob
import math
import mpmath
import numpy as np
from numpy.random import default_rng
import pdb
import matplotlib.pyplot as plt
import scipy
import scipy.special
import scipy.stats
import astropy
from astropy.io import fits
from astropy.nddata import StdDevUncertainty
from astropy.table import Table, join
import astropy.units as u
print(astropy.__version__)
# import spectrum
# print(spectrum.__version__)
# help(spectrum.Spectrum.__init__)
from specutils import Spectrum1D, SpectrumCollection
from specutils.manipulation import FluxConservingResampler
#dir(SpectrumCollection)


# # Data

# In[2]:


cat = Table.read('/its/research/astro/gama/loveday/Data/4MOST/OpR2exGalSample_170820_orig.csv')

for survey in ('s7', 's8'):
    sel = cat['SURVEY'] == survey
    print(survey, astropy.table.unique(cat[sel], keys='TEMPLATE')['TEMPLATE'])


# In[3]:


df = pd.read_csv('OpR2exGalSample_170820_orig.csv', usecols=['RULESET','SURVEY'])  

df


# In[5]:


df_filtered = df[df['SURVEY'] == 's7']

df_filtered

data_drop2 = df_filtered.drop_duplicates(subset=['RULESET'])
data_drop2.RULESET.values.tolist()


# In[4]:


data_drop = df.drop_duplicates(subset=['RULESET'])
data_drop.RULESET.values.tolist()


# In[6]:


dt = Table.read('/its/research/astro/gama/loveday/Data_local/4most/OpR2/spectab.fits')  

dt


# In[8]:


get_ipython().run_line_magic('cd', "'/its/research/astro/gama/loveday/Data_local/4most/OpR2'")
cat = Table.read('/research/astro/gama/loveday/Data/4MOST/OpR2exGalSample_170820_orig.csv')
# imin = np.argmin(cat['REDSHIFT'])
# print(cat[imin])
norm = True  # True to normalise spectra before combining
spectab = Table.read('spectab.fits')
fluxcon = FluxConservingResampler()
wmin, wmax = 3700, 9500  # Observed wavelength range
rng = np.random.default_rng()
nsum = 50  # Number of spectra to sum
for survey, tmplist in zip(('s7', 's8'), (('Pass', 'SF'), ('ELG', 'LRG', 'qso_BL'))):
    for tmp in tmplist:
        sel = (cat['SURVEY'] == survey) & (np.char.find(cat['TEMPLATE'], tmp) > -1) 
        nobj = len(cat[sel])
        zmin = np.min(cat[sel]['REDSHIFT'])
        zmax = np.max(cat[sel]['REDSHIFT'])
        lmin, lmax = wmin/(1 + zmax), wmax/(1 + zmin)  # Restframe wavelength range
        disp_grid = np.arange(lmin, lmax, 3) * u.AA
        print(len(disp_grid), nobj)
        flux = np.zeros((nobj, len(disp_grid)))
        ivar = np.zeros((nobj, len(disp_grid)))
        match = join(cat[sel], spectab, keys='CNAME', join_type='left')
        print(survey, tmp, len(cat[sel]), len(match), zmin, zmax, lmin, lmax)
        for i in range(len(match)):
            specname = match['SPECNAME'][i]
            z = match['REDSHIFT'][i]
            spect = Table.read(specname)
            flux_unit = spect['FLUX'].unit
            if norm:
                sumflux = np.sum(spect['FLUX'].data[0])
                spect['FLUX'].data[0] /= sumflux
                spect['ERR_FLUX'].data[0] /= sumflux
            spec = Spectrum1D(spectral_axis=spect['WAVE'].data[0] * spect['WAVE'].unit / (1 + z), 
                              flux=spect['FLUX'].data[0] * flux_unit,
                              uncertainty=StdDevUncertainty(spect['ERR_FLUX'].data[0]))
            specr = fluxcon(spec, disp_grid)
            flux[i] = specr.flux.value
            ivar[i] = specr.uncertainty.array**(-2)


# In[6]:


get_ipython().run_line_magic('cd', "'/its/research/astro/gama/loveday/Data_local/4most/OpR2'")
cat = Table.read('/research/astro/gama/loveday/Data/4MOST/OpR2exGalSample_170820_orig.csv')
# imin = np.argmin(cat['REDSHIFT'])
# print(cat[imin])
norm = True  # True to normalise spectra before combining
spectab = Table.read('spectab.fits')
fluxcon = FluxConservingResampler()
wmin, wmax = 3700, 9500  # Observed wavelength range
rng = np.random.default_rng()
nsum = 50  # Number of spectra to sum
for survey, tmplist in zip(('s7', 's8'), (('Pass', 'SF'), ('ELG', 'LRG', 'qso_BL'))):
    for tmp in tmplist:
     #   sel = (cat['SURVEY'] == survey) * (np.char.find(cat['TEMPLATE'], tmp) > -1)
        sel = (cat['SURVEY'] == survey) * (np.char.find(cat['TEMPLATE'], tmp) > -1) * cat['SPECNAME'].str.contains('LJ1')
        nobj = len(cat[sel])
        zmin = np.min(cat[sel]['REDSHIFT'])
        zmax = np.max(cat[sel]['REDSHIFT'])
        lmin, lmax = wmin/(1 + zmax), wmax/(1 + zmin)  # Restframe wavelength range
        disp_grid = np.arange(lmin, lmax, 3) * u.AA
        print(len(disp_grid), nsum)
        flux = np.zeros((nsum, len(disp_grid)))
        ivar = np.zeros((nsum, len(disp_grid)))
        match = join(cat[sel], spectab, keys='CNAME', join_type='left')
        print(survey, tmp, len(cat[sel]), len(match), zmin, zmax, lmin, lmax)
        speclist = rng.choice(len(match), nsum)
        for i in range(nsum):
            ispec = speclist[i]
            specname = match['SPECNAME'][ispec]
            z = match['REDSHIFT'][ispec]
            spect = Table.read(specname)
            flux_unit = spect['FLUX'].unit
            if norm:
                sumflux = np.sum(spect['FLUX'].data[0])
                spect['FLUX'].data[0] /= sumflux
                spect['ERR_FLUX'].data[0] /= sumflux
            spec = Spectrum1D(spectral_axis=spect['WAVE'].data[0] * spect['WAVE'].unit / (1 + z), 
                              flux=spect['FLUX'].data[0] * flux_unit,
                              uncertainty=StdDevUncertainty(spect['ERR_FLUX'].data[0]))
            specr = fluxcon(spec, disp_grid)

#             plt.clf()
#             plt.plot(spec.spectral_axis, spec.flux)
#             plt.plot(spec.spectral_axis, spec.uncertainty.array)
#             plt.plot(specr.spectral_axis, specr.flux)
#             plt.plot(specr.spectral_axis, specr.uncertainty.array**-0.5)
#             plt.xlabel('Wavelength [A]')
#             plt.ylabel('Flux [adu]')
#             plt.title(f'{specname} {z}')
#             plt.show()

            flux[i, :] = specr.flux
            ivar[i, :] = specr.uncertainty.array  # Note that rebinning changes uncvertainty to inverse variance
        flux = np.ma.masked_invalid(flux)
        ivar = np.ma.masked_invalid(ivar)
        flux_mean = np.ma.average(flux, weights=ivar, axis=0)

        plt.clf()
        plt.plot(disp_grid, flux_mean)
        plt.xlabel('Wavelength [A]')
        plt.ylabel(f'Flux [{flux_unit}]')
        plt.title(survey + tmp)
        plt.show()


# In[ ]:





# In[7]:


#Truncate wavelength range
def truncate_spectra_data(data, lower_bound, upper_bound):
    truncated_data = [(x, y) for x, y in data if lower_bound <= x <= upper_bound]
    return truncated_data

spectra_data = []

for survey, tmplist in zip(('s7', 's8'), (('Pass', 'SF'), ('ELG', 'LRG', 'qso_BL'))):
    for tmp in tmplist:
        sel = (cat['SURVEY'] == survey) * (np.char.find(cat['TEMPLATE'], tmp) > -1)
        nobj = len(cat[sel])
        zmin = np.min(cat[sel]['REDSHIFT'])
        zmax = np.max(cat[sel]['REDSHIFT'])
        lmin, lmax = wmin/(1 + zmax), wmax/(1 + zmin)  # Restframe wavelength range
        disp_grid = np.arange(lmin, lmax, 3) * u.AA
        flux = np.zeros((nsum, len(disp_grid)))
        ivar = np.zeros((nsum, len(disp_grid)))
        match = join(cat[sel], spectab, keys='CNAME', join_type='left')
        #print(survey, tmp, len(cat[sel]), len(match), zmin, zmax, lmin, lmax)
        speclist = rng.choice(len(match), nsum)
        for i in range(nsum):
            ispec = speclist[i]
            specname = match['SPECNAME'][ispec]
            z = match['REDSHIFT'][ispec]
            spect = Table.read(specname)
            flux_unit = spect['FLUX'].unit
            if norm:
                sumflux = np.sum(spect['FLUX'].data[0])
                spect['FLUX'].data[0] /= sumflux
                spect['ERR_FLUX'].data[0] /= sumflux
            spec = Spectrum1D(spectral_axis=spect['WAVE'].data[0] * spect['WAVE'].unit / (1 + z), 
                              flux=spect['FLUX'].data[0] * flux_unit,
                              uncertainty=StdDevUncertainty(spect['ERR_FLUX'].data[0]))
            specr = fluxcon(spec, disp_grid)

            flux[i, :] = specr.flux
            ivar[i, :] = specr.uncertainty.array  # Note that rebinning changes uncvertainty to inverse variance
        flux = np.ma.masked_invalid(flux)
        ivar = np.ma.masked_invalid(ivar)
        flux_mean = np.ma.average(flux, weights=ivar, axis=0)
        
        disp_grid = np.asarray(disp_grid)
        flux_mean = np.asarray(flux_mean)
        
        wavelength_flux = zip(disp_grid,flux_mean)
        lower_bound = 2500
        upper_bound = 4500
        truncated_data = truncate_spectra_data(wavelength_flux, lower_bound, upper_bound)
        #print(truncated_data)
        spectra_data.append(truncated_data)
        
        #wavelength = [x for x,y in truncated_data]
        #wavelength = np.asarray(wavelength)
        #print(wavelength)
        #flux = [y for x,y in truncated_data]
        #flux = np.asarray(flux)
        #print(flux)
        
####   BIN THE DATA
print(spectra_data)


# In[8]:


# Define the wavelength range
min_wavelength = 2500
max_wavelength = 4500

# Define the width of each bin
bin_width = 3

# Calculate the number of bins
num_bins = (max_wavelength - min_wavelength) // bin_width

print("Number of bins:", num_bins)


# In[9]:


# Determine the range of wavelengths across all spectra
min_wavelength = min(min(wavelength for wavelength, flux in spectrum) for spectrum in spectra_data)
max_wavelength = max(max(wavelength for wavelength, flux in spectrum) for spectrum in spectra_data)

print(min_wavelength)
# Define the width of each bin
bin_width = 3

# Calculate the number of bins
#num_bins = (max_wavelength - min_wavelength) // bin_width

#print("Number of bins:", num_bins)


# Define the width of each bin
bin_width = 3

# Calculate the number of bins
num_bins = int((max_wavelength - min_wavelength) / bin_width) + 1  # Adding 1 to include the upper bound

# Create the common wavelength grid
common_wavelengths = np.linspace(min_wavelength, max_wavelength, num_bins)

# Bin the fluxes onto the common wavelength grid for each spectrum
binned_fluxes = []
for spectrum in spectra_data:
    wavelengths, fluxes = zip(*spectrum)
    binned_flux, _ = np.histogram(wavelengths, bins=num_bins, range=(min_wavelength, max_wavelength), weights=fluxes)
    binned_fluxes.append(binned_flux)

print("Common Wavelengths:", common_wavelengths)
print("Binned Fluxes:", binned_fluxes)



#bin_edges = np.linspace(min_wavelength, max_wavelength, num_bins + 1)

# Bin the flux values onto the common wavelength grid
#binned_flux_values = []
#for wavelengths, flux_values in spectra_data:
 #   binned_flux, _ = np.histogram(wavelengths, bins=bin_edges, weights=flux_values)
  #  binned_flux_values.append(binned_flux)
X = np.array(binned_fluxes)


# In[ ]:


# Perform k-means clustering
kmeans = KMeans(n_clusters=2, random_state=42)  # You can change the number of clusters (n_clusters) as per your requirement
kmeans.fit(X)

# Get cluster centers
cluster_centers = kmeans.cluster_centers_

bin_edges = np.linspace(min_wavelength, max_wavelength, num_bins + 1)
# Plot the cluster centers
plt.figure(figsize=(10, 6))
for i, center in enumerate(cluster_centers):
    plt.plot(bin_edges[:-1], center, label=f'Cluster {i+1}')
plt.xlabel('Wavelength')
plt.ylabel('Flux')
plt.title('Cluster Centers of K-Means Clustering')
plt.legend()
plt.show()


# In[10]:


#X_transposed = X.T
# Perform PCA
pca = PCA()  # retains all principal components
X_pca = pca.fit_transform(X)

print(X.shape)
print(len(X_pca))
#print(X_transposed)

# Get variance ratio
variance_ratio = pca.explained_variance_ratio_
print(len(variance_ratio))

# Create a DataFrame to store the explained variance
variance_df = pd.DataFrame({'Principal Component': range(1, len(variance_ratio) + 1),
                            'Explained Variance Ratio': variance_ratio})

print(variance_df)


# In[ ]:


# Perform PCA
pca = PCA(n_components=2)  # Specify the number of components you want to keep
X_pca = pca.fit_transform(X)

# Plot PCA results
plt.figure(figsize=(8, 6))
plt.scatter(X_pca[:, 0], X_pca[:, 1], alpha=0.7)
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.title('PCA Plot of Spectra Data')
plt.grid(True)
plt.show()


# In[ ]:



# Perform PCA on the truncated spectra
pca = PCA() 
galaxy_pca = pca.fit_transform(X)

# Calculate the amount of variance explained by each principal component
variance_explained = pca.explained_variance_ratio_

# Create a table showing the amount of variance in the first 10 principal components
variance_table = pd.DataFrame({'Principal Component': range(1,667),
                               'Variance Explained': variance_explained})

print("Variance Explained by All Principal Components:")
print(variance_table)# PCA 


# In[ ]:


print(disp_grid)


# In[ ]:



n_samples = 1000
n_wavelengths = 200
galaxy_spectra = np.random.rand(n_samples, n_wavelengths)
print(galaxy_spectra)

wavelength_range = slice(50, 150)
galaxy_spectra_truncated = galaxy_spectra[:, wavelength_range]
print(len(galaxy_spectra_truncated))
        


# In[ ]:





# # PCA 

# In[ ]:


# Perform PCA on the truncated spectra
pca = PCA() 

galaxy_pca = pca.fit_transform(galaxy_spectra_truncated)
print(len(galaxy_pca))
print(galaxy_spectra_truncated.shape)
# Calculate the amount of variance explained by each principal component
variance_explained = pca.explained_variance_ratio_
print(len(variance_explained))
# Create a table showing the amount of variance in the first 10 principal components
variance_table = pd.DataFrame({'Principal Component': range(1,101),
                               'Variance Explained': variance_explained})

print("Variance Explained by All Principal Components:")
print(variance_table)# PCA 


# ### Plot of Variance in each Principal Component

# In[ ]:


# Plot the variance explained by each principal component
plt.figure(figsize=(8, 6))
plt.bar(range(1, 101), variance_explained, color='skyblue')
plt.title('Variance Explained by Each Principal Component')
plt.xlabel('Principal Component')
plt.ylabel('Variance Explained')
#plt.xticks(range(1, 101))
plt.grid(False)
plt.show()


# ### PCA for just first 2 components and Plot

# In[ ]:


# Perform PCA on the truncated spectra
pca = PCA(n_components=2)  # Consider only the first two principal components
galaxy_pca = pca.fit_transform(galaxy_spectra_truncated)
print(len(galaxy_pca))
print(len(galaxy_spectra_truncated))
variance_explained = pca.explained_variance_ratio_
print(len(variance_explained))
# Create a table showing the amount of variance in the first 10 principal components
variance_table = pd.DataFrame({'Principal Component': range(1, 3),
                               'Variance Explained': variance_explained})

print("Variance Explained by the First 2 Principal Components:")
print(variance_table)

# Plot the data as functions of the first two principal components
plt.figure(figsize=(8, 6))
plt.scatter(galaxy_pca[:, 0], galaxy_pca[:, 1], alpha=0.8)
plt.title('PCA')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.grid(True)
plt.show()


# In[ ]:


import numpy as np
from sklearn.decomposition import PCA
Y = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
pca = PCA(n_components=4)
pca.fit(Y)
PCA(n_components=2)
print(pca.explained_variance_ratio_)

print(pca.singular_values_)


# In[ ]:




