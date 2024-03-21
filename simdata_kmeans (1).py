#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from scipy.spatial import ConvexHull
import pandas as pd

import os
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


# In[4]:


df_filtered = df[df['SURVEY'] == 's7']

df_filtered

data_drop2 = df_filtered.drop_duplicates(subset=['RULESET'])
data_drop2.RULESET.values.tolist()


# In[5]:


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


# In[11]:


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


# # Start

# In[7]:


# Directory containing the spectra
directory = '/its/research/astro/gama/loveday/Data_local/4most/OpR2/20200907/singlespec'

files = os.listdir(directory)

# Print each file name
for file in files:
    print(file)


# In[8]:


# Directory containing the spectra
directory = ['20200907','20200908','20200909','20200910','20200912','20200913']
#directory = '/its/research/astro/gama/loveday/Data_local/4most/OpR2/20200907/singlespec'

files_with_LJ1 = []
fluxes = []
# Count how many spectra have 'LJ1' in their names
for i in directory:
    for file in os.listdir(f'/its/research/astro/gama/loveday/Data_local/4most/OpR2/{i}/singlespec'):
        if 'LJ1' in file and file.endswith('.fits'):
            files_with_LJ1.append(file)
     #       filepath = os.path.join(f'/its/research/astro/gama/loveday/Data_local/4most/OpR2/{i}/singlespec', file)
      #      with fits.open(filepath) as hdul:
       #         data = hdul[1].data  
        #        flux = data['FLUX']
         #       fluxes.append(flux)
                
                
num_files_with_LJ1 = len(files_with_LJ1)
print(f"Number of spectra with 'LJ1' in their names: {num_files_with_LJ1}")
#print(files_with_LJ1[20000:20050])

#flat_list = np.concatenate(fluxes).tolist()
#print(flat_list)


# In[26]:


all_fluxes = []
directories = ['/its/research/astro/gama/loveday/Data_local/4most/OpR2/20200907/singlespec',
               '/its/research/astro/gama/loveday/Data_local/4most/OpR2/20200908/singlespec',
               '/its/research/astro/gama/loveday/Data_local/4most/OpR2/20200909/singlespec',
               '/its/research/astro/gama/loveday/Data_local/4most/OpR2/20200910/singlespec',
               '/its/research/astro/gama/loveday/Data_local/4most/OpR2/20200912/singlespec',
               '/its/research/astro/gama/loveday/Data_local/4most/OpR2/20200913/singlespec']
# Load flux values for each spectrum
for directory in directories:
    files_with_LJ1 = [os.path.join(directory, file) for file in os.listdir(directory) if 'LJ1' in file and file.endswith('.fits')]
    for file_path in files_with_LJ1:
        with fits.open(file_path) as hdul:
            flux = hdul[1].data['FLUX'].flatten()
            all_fluxes.append(flux)
            wavelength = hdul[1].data['WAVE'].flatten()
wavelengths = []
wavelengths.append(wavelength)
data = np.array(all_fluxes)

# Perform PCA
pca = PCA(n_components=30)  # Specify the number of principal components
principal_components = pca.fit_transform(data)

# Accessing explained variance ratio
variance = pca.explained_variance_ratio_
#print("Explained Variance Ratio:", explained_variance_ratio)

# Plotting the principal components
plt.scatter(principal_components[:, 0], principal_components[:, 1])
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.title('Principal Component Analysis (PCA)')
plt.grid(True)
plt.show()


# In[27]:


# Create a DataFrame to store the explained variance
variance_df = pd.DataFrame({'Principal Component': range(1, len(variance) + 1),
                            'Explained Variance Ratio': variance})

variance_df


# In[29]:


# Plot the variance explained by each principal component
plt.figure(figsize=(8, 6))
plt.bar(range(1, 31), variance[:30], color='skyblue')
plt.title('Variance Explained by Each Principal Component')
plt.xlabel('Principal Component')
plt.ylabel('Variance Explained')
#plt.xticks(range(1, 101))
plt.grid(False)
plt.show()


# In[30]:


Cumalative_var = pca.explained_variance_ratio_.cumsum()
# Create a DataFrame to store the explained variance
C_variance_df = pd.DataFrame({'Principal Component': range(1, len(Cumalative_var) + 1),
                            'Explained Variance Ratio': Cumalative_var})

C_variance_df


# In[31]:


# Compute the cumulative explained variance ratio
Cumulative_variance = np.cumsum(variance)

# Plot the cumulative variance explained by each component
plt.plot(np.arange(1, len(Cumulative_variance) + 1), Cumulative_variance, marker='o')
plt.xlabel('Number of Components')
plt.ylabel('Cumulative Variance %')
plt.title('Cumulative Variance Explained by Principal Components')
plt.grid(True)
plt.axhline(y=0.8, color='grey', linestyle='--')
plt.text(140, 0.81, '80% cut-off threshold (PC29)', color = 'black', fontsize=12)
plt.show()


# In[18]:


# Keep only the first 16 principal components
selected_principal_components = principal_components[:, :16]

# Create a 4x4 grid of subplots for each principal component
fig, axes = plt.subplots(4, 4, figsize=(12, 12))

# Plot each principal component against the original data
for i, ax in enumerate(axes.flat):
    ax.scatter(data[:, 0], data[:, 1], c=selected_principal_components[:, i], cmap='viridis', marker='o')
    ax.set_title(f'PC {i+1}')
    ax.set_xlabel('Feature 1')
    ax.set_ylabel('Feature 2')
    ax.grid(True)

# Adjust layout to prevent overlap of labels
plt.tight_layout()
plt.show()


# In[32]:


#fig = plt.figure(dpi=150)
#grid = gridspec.GridSpec(3,1, height_ratios=(2,0.4, 0.4))
#fig.suptitle(r"$ ^{87}Rb$ $F_{g}=1$ Spectral Lines")
#ax1 = fig.add_subplot(grid[0])
#ax2 = fig.add_subplot(grid[1])
#ax3 = fig.add_subplot(grid[2])
#grid.update(hspace=0)


num_components = 5

# Scale the principal component data
scaler = StandardScaler()
scaled_principal_components = scaler.fit_transform(principal_components)

# Create subplots
fig, axes = plt.subplots(num_components, 1, figsize=(8, 12), dpi=150, sharex=True)

# Plot each scaled principal component
for i in range(num_components):
    axes[i].plot(wavelengths, scaled_principal_components[:, i], label=f'PC {i+1}')
    axes[i].set_ylabel('Scaled Flux')
    axes[i].legend()

# Set x-axis label for the last subplot
axes[num_components-1].set_xlabel('Wavelength')

# Adjust layout to prevent overlap of labels
plt.tight_layout()

plt.show()


# In[37]:


from matplotlib import gridspec

num_components = 5

# Scale the principal component data
scaler = StandardScaler()
scaled_principal_components = scaler.fit_transform(principal_components)
print(scaled_principal_components.shape)
# Create a grid of subplots
fig = plt.figure(figsize=(8, 12), dpi=150)
grid = gridspec.GridSpec(num_components, 1, height_ratios=[2]*num_components)

# Plot each scaled principal component in a subplot within the grid
for i in range(num_components):
    ax = fig.add_subplot(grid[i])
    ax.plot(scaled_principal_components[:, i], label=f'PC {i+1}')
    ax.set_ylabel('Scaled Flux')
    ax.legend()

# Set x-axis label for the last subplot
ax.set_xlabel('Wavelength')

# Adjust layout to prevent overlap of labels
plt.tight_layout()

plt.show()


# In[38]:


from sklearn.preprocessing import StandardScaler

def transform_pca(X, n):

    pca = PCA(n_components=n)
   # pca.fit(X)
    X_new = pca.inverse_transform(pca.fit_transform(X))

    return X_new

rows = 4
cols = 4
comps = 1
scaler = StandardScaler()
X_scaled = scaler.fit_transform(data)

fig, axes = plt.subplots(rows, 
                         cols, 
                         figsize=(12,12), 
                         sharex=True, 
                         sharey=True)


for row in range(rows):
    for col in range(cols):
        
        X_new = transform_pca(X_scaled, comps)
        axes[row, col].scatter(X_scaled[:, 0], X_scaled[:, 1], c='gray', alpha=0.3)
            
            # Plot transformed data (principal components) in black
        axes[row, col].scatter(X_new[:, 0], X_new[:, 1], c='black')
            
        axes[row, col].set_title(f'PCA Components: {comps}')
        
        comps += 1
        #except:
         #   pass
plt.xlim(0,10)
plt.ylim(0,20)
plt.tight_layout()
#plt.savefig('pcavisualize_2.png', dpi=300)


# # K-means

# In[49]:


# Perform k-means clustering on the PCA-transformed data
n_clusters = 5  # Number of clusters
kmeans = KMeans(n_clusters=n_clusters)
kmeans.fit(principal_components[:, :29])

# Count the number of data points in each cluster
cluster_counts = np.bincount(kmeans.labels_)

# Print the number of data points in each cluster
for i, count in enumerate(cluster_counts):
    print(f"Cluster {i}: {count} data points")


# In[50]:


# Plot the data with clusters
plt.figure(figsize=(8, 6))
plt.scatter(principal_components[:, 0], principal_components[:, 1], c=kmeans.labels_, cmap='Set2', alpha=0.8)
plt.scatter(kmeans.cluster_centers_[:, 0], kmeans.cluster_centers_[:, 1], marker='x', s=100, c='red', label='Centroids')


# In[11]:


wavelengths = []
fluxes = []
directories = ['/its/research/astro/gama/loveday/Data_local/4most/OpR2/20200907/singlespec',
               '/its/research/astro/gama/loveday/Data_local/4most/OpR2/20200908/singlespec',
               '/its/research/astro/gama/loveday/Data_local/4most/OpR2/20200909/singlespec',
               '/its/research/astro/gama/loveday/Data_local/4most/OpR2/20200910/singlespec',
               '/its/research/astro/gama/loveday/Data_local/4most/OpR2/20200912/singlespec',
               '/its/research/astro/gama/loveday/Data_local/4most/OpR2/20200913/singlespec']
# Loop through each directory
for directory in directories:
    # Search for FITS files containing 'LJ1' in their names
    files_with_LJ1 = [file for file in os.listdir(directory) if 'LJ1' in file and file.endswith('.fits')]
    
    # Print the number of spectra with 'LJ1' in their name for each directory
   # print(f"Number of spectra with 'LJ1' in their name in {directory}: {len(files_with_LJ1)}")
    
    # Loop through each spectrum in the current directory
    for file in files_with_LJ1[:2]:
        with fits.open(os.path.join(directory, file)) as hdul:
            data = hdul[1].data  # Assuming the data is in the first extension
            wavelength = data['WAVE'].flatten()  # Flatten the wavelength array
            flux = data['FLUX'].flatten()  # Flatten the flux array
            
            # Store flattened wavelength and flux values for each spectrum
            wavelengths.append(wavelength)
            fluxes.append(flux)
            
            
# Print the number of spectra with 'LJ1' in their name for each directory
print(f"Number of spectra with 'LJ1' in their name in {directory}: {len(files_with_LJ1)}")

flat_list = np.concatenate(fluxes).tolist()
print(flat_list)


# In[10]:


# Plot the first 10 spectra

fluxes = []

for i in directory:
    for file in files_with_LJ1[:50]:
        with fits.open(os.path.join(f'/its/research/astro/gama/loveday/Data_local/4most/OpR2/{i}/singlespec', file)) as hdul:
    #hdul = fits.open(file)
    #data = fits.getdata(f'{image_dir}/{f}_sci.fits') # read FITS file data into numpy array

            data = hdul[1].data  # Assuming the data is in the first extension
        #wavelength = data['WAVE'] # Assuming the wavelength is stored in the 'WAVE' column
        #print(wavelength)
            flux = data['FLUX'] # Assuming the flux is stored in the 'FLUX' column
        #print(len(flux))
        #plt.plot(wavelength, flux, label=file)
        #for i in wavelength:
       # wavelengths.append(wavelength)
        #for i in flux:
            fluxes.append(flux)

#print(wavelengths)
#print(fluxes)
#plt.plot(wavelengths,fluxes)  

#fluxes = np.asarray(fluxes)

flat_list = np.concatenate(fluxes).tolist()
print(flat_list)


# In[19]:


# Read data from 'spectab.fits'
#spectab = Table.read('/its/research/astro/gama/loveday/Data_local/4most/OpR2/20200907')

# Directory containing the spectra
directory = ['20200907','20200908','20200909','20200910','20200912','20200913']
#directory = '/its/research/astro/gama/loveday/Data_local/4most/OpR2/20200907/singlespec'
wavelengths = []
fluxes = []
# Count how many spectra have 'LJ1' in their names
for i in directory:
    files_with_LJ1 = [file for file in os.listdir(f'/its/research/astro/gama/loveday/Data_local/4most/OpR2/{i}/singlespec') if 'LJ1' in file and file.endswith('.fits')]

#num_files_with_LJ1 = len(files_with_LJ1)
#print(f"Number of spectra with 'LJ1' in their names: {num_files_with_LJ1}")

# Iterate over each FITS file
#for fits_file in files_with_LJ1:
    # Open the FITS file
   # with fits.open(os.path.join(directory, fits_file)) as hdul:
        # Print the column names (headings) of the first HDU (extension)
 #   print(f"\nColumns for {fits_file}:")
  #  print(hdul[1].columns.names)
    
    
# Plot the first 10 spectra
#for file in files_with_LJ1[:1000]:
    with fits.open(os.path.join(directory, file)) as hdul:
    #hdul = fits.open(file)
    #data = fits.getdata(f'{image_dir}/{f}_sci.fits') # read FITS file data into numpy array

        data = hdul[1].data  # Assuming the data is in the first extension
        wavelength = data['WAVE'] # Assuming the wavelength is stored in the 'WAVE' column
        #print(wavelength)
        flux = data['FLUX'] .flatten() # Assuming the flux is stored in the 'FLUX' column
        #print(len(flux))
        #plt.plot(wavelength, flux, label=file)
        #for i in wavelength:
       # wavelengths.append(wavelength)
        #for i in flux:
        fluxes.append(flux)

#print(wavelengths)
#print(fluxes)
#plt.plot(wavelengths,fluxes)  

fluxes = np.asarray(fluxes)
#print(flux)
# Customize the plot
#plt.xlabel('Wavelength')
#plt.ylabel('Flux')
#plt.title('Flux vs. Wavelength for the First 10 Spectra with LJ1')
#plt.legend()
#plt.grid(True)

# Show the plot
#plt.show()


# In[20]:


pca = PCA(n_components=50)  # retains all principal components
flux_pca = pca.fit_transform(fluxes)

#print(X.shape)
#print(len(flux_pca))
#print(X_transposed)
print(fluxes.shape)
# Get variance ratio
variance_ratio = pca.explained_variance_ratio_
print(len(variance_ratio))

# Create a DataFrame to store the explained variance
variance_df = pd.DataFrame({'Principal Component': range(1, len(variance_ratio) + 1),
                            'Explained Variance Ratio': variance_ratio})

print(variance_df)


# In[17]:


flist = ('20200907/singlespec/qmost_14031061-0228192_20200907_100150_LR1.fits',
        '20200907/singlespec/qmost_14031061-0228192_20200907_100150_LG1.fits',
        '20200907/singlespec/qmost_14031061-0228192_20200907_100150_LB1.fits',
        '20200907/singlespec/qmost_14031061-0228192_20200907_100150_LJ1.fits',
        '20200907/singlespec/qmost_14025043-0246181_20200907_100128_LR1.fits',
        '20200907/singlespec/qmost_14025043-0246181_20200907_100128_LG1.fits',
        '20200907/singlespec/qmost_14025043-0246181_20200907_100128_LB1.fits',
        '20200907/singlespec/qmost_14025043-0246181_20200907_100128_LJ1.fits',
        '20200907/singlespec/qmost_14025043-0246181_20200907_100150_LR1.fits',
        '20200907/singlespec/qmost_14025043-0246181_20200907_100150_LG1.fits',
        '20200907/singlespec/qmost_14025043-0246181_20200907_100150_LB1.fits',
        '20200907/singlespec/qmost_14025043-0246181_20200907_100150_LJ1.fits',
        '20200907/singlespec/qmost_14023696-0249329_20200907_100104_LR1.fits',
        '20200907/singlespec/qmost_14023696-0249329_20200907_100104_LG1.fits',
        '20200907/singlespec/qmost_14023696-0249329_20200907_100104_LB1.fits',
        '20200907/singlespec/qmost_14023696-0249329_20200907_100104_LJ1.fits',
        '20200907/singlespec/qmost_14023696-0249329_20200907_100128_LR1.fits',
        '20200907/singlespec/qmost_14023696-0249329_20200907_100128_LG1.fits',
        '20200907/singlespec/qmost_14023696-0249329_20200907_100128_LB1.fits',
        '20200907/singlespec/qmost_14023696-0249329_20200907_100128_LJ1.fits')
allfiles = glob.glob('20200907/singlespec/*')
get_ipython().run_line_magic('cd', "'/research/astro/gama/loveday/4most/OpR2'")
print(allfiles)
for f in allfiles:
   # with fits.open(f) as hdul:
#         hdul.verify()
#         print(hdul[0].header['OBJECT'])
    spec = Table.read(f)
    print(spec)
    plt.clf()
    plt.plot(spec['WAVE'].data[0], spec['FLUX'].data[0])
    if 'LJ1' in f:
        plt.plot(spec['WAVE'].data[0], spec['ERR_FLUX'].data[0])
    else:
        plt.plot(spec['WAVE'].data[0], spec['ERR_IVAR'].data[0]**-0.5)
# #         plt.plot(spec['WAVE'].data[0], spec['SENSFUNC'].data[0])


            
    plt.xlabel('Wavelength [A]')
    plt.ylabel('Flux [adu]')
    plt.title(f)

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




