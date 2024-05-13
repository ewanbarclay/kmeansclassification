#!/usr/bin/env python
# coding: utf-8

# In[1]:


import csv
import re
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
from specutils import Spectrum1D, SpectrumCollection
from specutils.manipulation import FluxConservingResampler
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score
from sklearn.metrics import silhouette_samples, silhouette_score
from matplotlib import cm
import tarfile


# In[2]:


# Function to extract flux data from FITS file
def extract_flux(fits_filename):
    with fits.open(fits_filename) as hdul:
        flux = hdul[1].data['FLUX'].flatten()
    return flux


# In[3]:


# Function to extract flux data from FITS file
def extract_flux_wave(fits_filename):
    with fits.open(fits_filename) as hdul:
        flux = hdul[1].data['FLUX'].flatten() 
        wavelengths = hdul[1].data['WAVE'].flatten()
    return flux, wavelengths


# In[4]:


all_fluxes = []
data_type = []
root_folder = '/its/home/ecb41/4MOST/files2/Galaxies'


# In[5]:


lj1_redshifts = {}

# Open the text file and extract LJ1 .fits filenames and their corresponding redshift values
with open('/its/mnt/lustre/scratch/astro/loveday/4most/OpR2.5/opr25_metadata_05042022.txt', 'r') as file:
    reader = csv.reader(file, delimiter=' ')
    for row in reader:
        if "LJ1" in row[0] and row[0].endswith(".fits"):
            other_fits = [f for f in row[1:] if f.endswith(".fits")][0]
          #  print("Other fits:", other_fits)  
            # Search for 'redshift' 
            redshift_match = re.search(r'redshift([\d.]+)\.fits', other_fits)
            if redshift_match:
                redshift_value = float(redshift_match.group(1))
                lj1_redshifts[row[0]] = redshift_value
                continue  # Move to the next row

           
            # Search for 'z'
            z_match = re.search(r'z([\d]{1,2}(?:pt\d{2})?)', other_fits)
            if z_match:
                redshift_str = z_match.group(1).replace('pt', '.')
                redshift_value = float(redshift_str)  
                lj1_redshifts[row[0]] = redshift_value
                continue  

            # Search for 'z-'
            z_minus_match = re.search(r'z-([\d]+\.\d{2})', other_fits)
            if z_minus_match:
                redshift_value = float(z_minus_match.group(1))
                lj1_redshifts[row[0]] = redshift_value
                continue 

print("lj1_redshifts:")
#for lj1_file, redshift_value in lj1_redshifts.items():
 #   print(f"{lj1_file}: {redshift_value}")
print(len(lj1_redshifts.items()))


# In[6]:


count = 0
all_spectra = []
redshift = []
# Iterate over each folder containing tar files
for folder_name in os.listdir(root_folder):
    folder_path = os.path.join(root_folder, folder_name)
    if not os.path.isdir(folder_path):
        continue
    
    # Find the tar files in the folder
    tar_files = [filename for filename in os.listdir(folder_path) if filename.endswith('.tar.gz')]
    if not tar_files:
        continue
    
    # Iterate over each tar file
    for tar_filename in tar_files:
        tar_filepath = os.path.join(folder_path, tar_filename)
        
        
        with tarfile.open(tar_filepath, 'r:gz') as tar:
            # Extract all FITS files from the tar file
            fits_filenames = [member for member in tar.getmembers() if member.name.endswith('.fits') and "LJ1" in member.name]
            total_fits_files = len(fits_filenames)
            processed_count = 0
            
            # Extract each FITS file
            for fits_member in fits_filenames:
                processed_count += 1
                # Check if the FITS file matches with LJ1 filenames in lj1_redshifts
                if fits_member.name in lj1_redshifts:
                    
                    # Match LJ1 filename with its corresponding redshift value
                    redshift_value = lj1_redshifts[fits_member.name]
                    #print(f"Redshift: {redshift_value}")
                    redshift.append(redshift_value)
                    if redshift_value <= 0.25 :
                    
                        tar.extract(fits_member)
                    
                        # Extract flux and wavelength data
                        flux, wavelengths = extract_flux_wave(fits_member.name)
                        
                        # Convert flux to float32
                        flux = flux.astype(np.float32)
                        
                        # Correct wavelength for redshift
                        rest_frame_wavelength = wavelengths / (1 + redshift_value)
                    
                        # Create Spectrum1D object
                        spectrum = Spectrum1D(spectral_axis=rest_frame_wavelength * u.Angstrom, flux= flux * u.Unit('adu'))
                    
                   
                        all_spectra.append(spectrum)
                        
                    
                        data_type.append(folder_name)
                    
                        count += 1
                        print(f"Extracted data from {fits_member.name}. Total extracted: {count}")
                    
                        # Delete the extracted FITS file -> save space
                        os.remove(fits_member.name)
                    
                # Check if all FITS files in the tar file have been processed
                if processed_count >= total_fits_files:
                    break


# In[7]:


root_folder2 = '/its/home/ecb41/4MOST/files2/Stars'


# In[8]:


count = 0

# Iterate over each folder containing tar files in root_folder2
for folder_name in os.listdir(root_folder2):
    folder_path = os.path.join(root_folder2, folder_name)
    if not os.path.isdir(folder_path):
        continue
    
    # Find the tar files in the folder
    tar_files = [filename for filename in os.listdir(folder_path) if filename.endswith('.tar.gz')]
    if not tar_files:
        continue
    
    # Iterate over each tar file
    for tar_filename in tar_files:
        tar_filepath = os.path.join(folder_path, tar_filename)
        
        with tarfile.open(tar_filepath, 'r:gz') as tar:
            # Extract all FITS files from the tar file
            fits_filenames = [member for member in tar.getmembers() if member.name.endswith('.fits')  and "LJ1" in member.name]
            total_fits_files = len(fits_filenames)
            processed_count = 0
            
            # Extract each FITS file
            for fits_member in fits_filenames:
                processed_count += 1
                tar.extract(fits_member)
                
                # Extract flux and wavelength data
                flux, wavelengths = extract_flux_wave(fits_member.name)
                
                # Convert flux to float32
                flux = flux.astype(np.float32)
                
                # Create Spectrum1D object
                spectrum = Spectrum1D(spectral_axis=wavelengths * u.Angstrom, flux=flux * u.Unit('adu'))
                
                # Append the spectrum to the list
                all_spectra.append(spectrum)
                
                data_type.append(folder_name)
                
                # Increment the count
                count += 1
                print(f"Extracted data from {fits_member.name}. Total extracted: {count}")
                
                # Delete the extracted FITS file -> save space
                os.remove(fits_member.name)
                    
                # Check if all FITS files in the tar file have been processed
                if processed_count >= total_fits_files:
                    break


# In[9]:


print(len(all_spectra))


# ### Redshift Distribution

# In[10]:


for _ in range(19789):
    redshift.append(0)


# In[11]:


print(len(redshift))
print(max(redshift))


# In[12]:


#plot distribution
plt.hist(redshift, bins=30, color='skyblue', edgecolor= 'none')
plt.xlabel('z')
plt.ylabel('#')
#plt.title('Distribution of Redshift Values')
plt.axvline(x=0.25, color='red', linestyle='--') #cut-off
plt.text(0.28,10100, 'Cut-off at z=0.25', color='red', rotation=90, verticalalignment='bottom', fontsize=14)
plt.grid(True)
plt.savefig('figures/redshift')
plt.show()


# ## Resample spectra to new wavelength range

# In[13]:


min_wavelength = 3800
max_wavelength = 7500
num_points = 3701

common_wavelength_grid = np.linspace(min_wavelength, max_wavelength, num=num_points) * u.Angstrom

resampled_spectra = []

fluxcon = FluxConservingResampler()

count = 0
# Resample each spectrum to common wavelength grid
for spectrum in all_spectra:
    resampled_spectrum = fluxcon(spectrum, common_wavelength_grid)
    resampled_spectra.append(resampled_spectrum)
    count+=1
    print(count)


# In[14]:


# Extract resampled fluxes from the resampled spectra
resampled_fluxes = [resampled_spectrum.flux.value for resampled_spectrum in resampled_spectra]
flux_data = np.array([arr.astype(np.float32) for arr in resampled_fluxes])
#flux_data = resampled_fluxes.astype(np.float32)
#print(resampled_fluxes.shape)
# Convert the list of flux arrays into a 2D numpy array
#flux_data = np.array(resampled_fluxes)
print(flux_data.shape)
#print(resampled_fluxes.shape)


# In[15]:


# Check for NaN values
nan_indices = np.isnan(flux_data)
num_nan_values = np.sum(nan_indices)

# Print the number of NaN values
print(f"Number of NaN values in flux data: {num_nan_values}")


# In[16]:


# Normalize flux data
scaler = StandardScaler()
normalized_flux = scaler.fit_transform(flux_data)


# ## PCA

# In[17]:


pca = PCA(n_components=30)  # Specify the number of principal components
principal_components = pca.fit_transform(normalized_flux)

# Accessing explained variance ratio
variance = pca.explained_variance_ratio_
#print("Explained Variance Ratio:", explained_variance_ratio)

# Plotting the principal components
plt.scatter(principal_components[:, 0], principal_components[:, 1],s=2)
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
#plt.xlim(-13e-17,5e-17)
#plt.ylim(-12e-17,5e-17)
#plt.title('Principal Component Analysis (PCA)')
plt.grid(True)
plt.savefig('figures/pca')
plt.show()


# In[19]:


# Define the range for PC1 and PC2
pc1_range = (300,310)  
pc2_range = (95, 105)  

# Find indices of data points within the specified range
indices_within_range = np.where(
    (principal_components[:, 0] >= pc1_range[0]) &
    (principal_components[:, 0] <= pc1_range[1]) &
    (principal_components[:, 1] >= pc2_range[0]) &
    (principal_components[:, 1] <= pc2_range[1])
)[0]

# Print positions in the list of spectra for data points within the range
for index in indices_within_range:
    print("Position in the list of spectra: {}".format(index))


# In[23]:


# Define the range for PC1 and PC2
pc1_range = (200,208) 
pc2_range = (-105, -100)  

# Find indices of data points within the specified range
indices_within_range = np.where(
    (principal_components[:, 0] >= pc1_range[0]) &
    (principal_components[:, 0] <= pc1_range[1]) &
    (principal_components[:, 1] >= pc2_range[0]) &
    (principal_components[:, 1] <= pc2_range[1])
)[0]

# Print positions in the list of spectra for data points within the range
for index in indices_within_range:
    print("Position in the list of spectra: {}".format(index))


# In[28]:


# Function to plot flux vs. wavelength
def plot_spectrum(wavelength, flux, filename, save_path):
    plt.plot(wavelength, flux)
    plt.xlabel('Wavelength / Å')
    plt.ylabel('Flux / $\mathrm{erg \, s^{-1} \, cm^{-2} \, \AA^{-1}}$')
    plt.savefig(save_path)
   # plt.title(filename)
    plt.show()

tar_file_path = '/its/home/ecb41/4MOST/files2/Stars/F/FGKTeff_6500-6700_singlespec.tar.gz'

# List of filenames to search for
target_filenames = ['20220905/singlespec/qmost_20552700-4134000_20220905_100104_LJ1.fits', '20220901/singlespec/qmost_20555350-4046000_20220901_100008_LJ1.fits']


spectra_data = []

# Open the tar file
with tarfile.open(tar_file_path, 'r:gz') as tar:
    # Iterate over members (files) in the tar file
    for member in tar.getmembers():
        # Check if the member is a FITS file and matches one of the target filenames
        if member.isfile() and member.name in target_filenames and member.name.endswith('.fits'):
            # Extract the FITS file
            fits_file = tar.extractfile(member)
            
            with fits.open(fits_file) as hdul:
                # Extract wavelength and flux data
                wavelength = hdul[1].data['WAVE'].flatten()
                flux = hdul[1].data['FLUX'].flatten()
                # Store the spectra data
                spectra_data.append((wavelength, flux, member.name))
            # Break if both target files are found
            if len(spectra_data) == len(target_filenames):
                break

for i, data in enumerate(spectra_data, start=1):
    save_path = f'figure_{i}.png'  # Define the save path for the figure
    plot_spectrum(data[0], data[1], data[2], save_path)


# In[131]:


# Variance in each PC
variance_df = pd.DataFrame({'Principal Component': range(1,11),
                            'Explained Variance Ratio': variance[:10]})

variance_df


# In[132]:


print(variance_df.to_latex())


# In[133]:


# Compute the cumulative explained variance ratio
Cumulative_variance = pca.explained_variance_ratio_.cumsum()

# Plot the cumulative variance explained by each component
plt.figure(figsize=(8, 6))
plt.plot(np.arange(1,11), Cumulative_variance[:10], marker='o')
plt.xlabel('Principal Component')
plt.ylabel('Cumulative Variance %')
#plt.title('Cumulative Variance Explained by Principal Components')
plt.grid(True)
#plt.axhline(y=0.9, color='grey', linestyle='--')
#plt.text(140, 0.81, '80% cut-off threshold', color = 'black', fontsize=12)
plt.savefig("figures/Cumulativevariance")
plt.show()


# ### Plot a eigenspectrum as function of PCs

# In[135]:


min_wavelength = 3800
max_wavelength = 7500
num_points = 3701

common_wavelength_grids = np.linspace(min_wavelength, max_wavelength, num=num_points)

# Retrieve the coefficients of PC1 to PC5 from the PCA results
pc_coefficients = pca.components_[:5]

# Create a figure with subplots for each principal component
subfig = 3
fig, axs = plt.subplots(5, 1, figsize=(10,9))


synthetic_spectrum_pc = np.zeros_like(common_wavelength_grids)
    
    # Multiply each coefficient by the corresponding flux values for each wavelength
for i in range(5):
    # Plot the synthetic spectrum for the current PC
    axs[i].plot(common_wavelength_grids, pca.components_[i])
    if i < len(axs) - 1:
        
        axs[i].set_xticklabels([]) 

    if i == len(axs) - 1:
        axs[i].set_xlabel('Wavelength / Å', fontsize =12)
 
fig.text(0.0001, 0.52, 'Flux', va='center', rotation='vertical', fontsize = 12)

# Set labels for each subfigure indicating the corresponding principal component
for i, ax in enumerate(axs):
    ax.text(0.95, 0.8, f'PC{i+1}', transform=ax.transAxes, ha='left', fontsize =16)

plt.tight_layout()
plt.savefig('figures/eigenspectra')
plt.show()


# ### Plot data as function of angle in PC plot

# In[363]:


data_type_count = {}

# Iterate over the data_type list and count occurrences of each data type
for data_type in data_type:
    data_type_count[data_type] = data_type_count.get(data_type, 0) + 1

for data_type, count in data_type_count.items():
    print(f"{data_type}: {count}")


# In[107]:


type_color_map = {
    'whitedwarf': 'red',
    'F': 'blue',
    'G': 'green',
    'K': 'orange',
    'SF': 'purple',
    'Pass': 'yellow',
    'SN': 'cyan'
}


# In[108]:


# Initialize a dictionary to store the count of each data type
data_type_count = {}

# Iterate over the data_type list and count occurrences of each data type
for dtype in data_type:
    data_type_count[dtype] = data_type_count.get(dtype, 0) + 1

legend_handles = []
legend_labels = []

for dtype, color in type_color_map.items():
    count = data_type_count.get(dtype, 0)
    legend_handles.append(plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10))
    legend_labels.append(f"{dtype} ({count})")
    
    indices = [i for i, dt in enumerate(data_type) if dt == dtype]
    plt.scatter(principal_components[indices, 0], principal_components[indices, 1], c=color, s=0.4)

plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.grid(True)
plt.legend(handles=legend_handles, labels=legend_labels, loc='best',fontsize=8)
plt.savefig('figures/pca_types')
plt.show()


# In[218]:


data_type_count = {}

# Iterate over the data_type list and count occurrences of each data type
for dtype in data_type:
    data_type_count[dtype] = data_type_count.get(dtype, 0) + 1

legend_handles = []
legend_labels = []

for dtype, color in type_color_map.items():
    count = data_type_count.get(dtype, 0)
    legend_handles.append(plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10))
    legend_labels.append(f"{dtype} ({count})")

    indices = [i for i, dt in enumerate(data_type) if dt == dtype]
    plt.scatter(principal_components[indices, 0], principal_components[indices, 1], c=color, s=0.2)

plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.grid(True)
plt.legend(handles=legend_handles, labels=legend_labels, loc='best',fontsize=8)
plt.xlim(-25,20)
plt.ylim(-20,20)
plt.savefig('figures/pca_types_zoom')
plt.show()


# In[256]:


print(min(principal_components[:, 0]))


# ## Elbow / Dunn / Silhouette - angle

# In[ ]:


# Dunn index
def dunn_index(data, labels, centroids):
    min_intercluster_distance = np.inf
    max_intracluster_distance = -np.inf

    # Max intracluster distance
    for i in range(len(centroids)):
        intracluster_distance = np.max(np.linalg.norm(data[labels == i] - centroids[i], axis=1))
        if intracluster_distance > max_intracluster_distance:
            max_intracluster_distance = intracluster_distance

    # Min intercluster distance
    for i in range(len(centroids)):
        for j in range(i + 1, len(centroids)):
            intercluster_distance = np.linalg.norm(centroids[i] - centroids[j])
            if intercluster_distance < min_intercluster_distance:
                min_intercluster_distance = intercluster_distance

    dunn_index = min_intercluster_distance / max_intracluster_distance
    return dunn_index


# ELBOW METHOD AGAIN + DUNN INDEX:

# range of clusters
min_clusters = 1
max_clusters = 20

wcss = []
dunn_indices = []

# Iterate over the range of clusters
for n_clusters in range(min_clusters, max_clusters + 1):
    kmeans = KMeans(n_clusters=n_clusters, init='k-means++', random_state=42)
    kmeans.fit(theta)  
    wcss.append(kmeans.inertia_)  
    
    centroids = kmeans.cluster_centers_
    labels = kmeans.labels_
    dunn_index_value = dunn_index(theta, labels, centroids)
    dunn_indices.append(dunn_index_value) 

# Elbow curve
fig, ax1 = plt.subplots()
color = 'tab:red'
ax1.set_xlabel('Number of Clusters')
ax1.set_ylabel('WCSS', color=color)
ax1.plot(range(min_clusters, max_clusters + 1), wcss, marker='o', linestyle='-', color=color)
ax1.tick_params(axis='y', labelcolor=color)
#plt.title('Elbow Method')
plt.xticks(range(min_clusters, max_clusters + 1))
plt.grid(True)
plt.savefig('figures/Elbow_C_angle')

# Dunn index
fig, ax2 = plt.subplots()
color = 'tab:blue'
ax2.set_xlabel('Number of Clusters')
ax2.set_ylabel('Dunn Index', color=color)
ax2.plot(range(min_clusters, max_clusters + 1), dunn_indices, marker='o', linestyle='-', color=color)
ax2.tick_params(axis='y', labelcolor=color)
#plt.title('Dunn Index')
plt.xticks(range(min_clusters, max_clusters + 1))
plt.grid(True)
plt.savefig('figures/Dunn_C_angle')
plt.show()


# In[ ]:


range_n_clusters = [3,4,5,6,7]
X = theta

# Create a single figure for all silhouette plots
#fig, axs = plt.subplots(len(range_n_clusters), 2, figsize=(18, len(range_n_clusters) * 7))

for n_clusters in range_n_clusters:
    # Create a subplot with 1 row and 2 columns
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_size_inches(18, 7)

    # The 1st subplot is the silhouette plot
    # The silhouette coefficient can range from -1, 1 but in this example all
    # lie within [-0.1, 1]
    ax1.set_xlim([-0.1, 1])
    # The (n_clusters+1)*10 is for inserting blank space between silhouette
    # plots of individual clusters, to demarcate them clearly.
    ax1.set_ylim([0, len(X) + (n_clusters + 1) * 10])

    # Initialize the clusterer with n_clusters value and a random generator
    # seed of 10 for reproducibility.
    clusterer = KMeans(n_clusters=n_clusters, random_state=10)
    cluster_labels = clusterer.fit_predict(X)

    # The silhouette_score gives the average value for all the samples.
    # This gives a perspective into the density and separation of the formed
    # clusters
    silhouette_avg = silhouette_score(X, cluster_labels)
    print(
        "For n_clusters =",
        n_clusters,
        "The average silhouette_score is :",
        silhouette_avg,
    )

    # Compute the silhouette scores for each sample
    sample_silhouette_values = silhouette_samples(X, cluster_labels)

    y_lower = 10
    for i in range(n_clusters):
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = sample_silhouette_values[cluster_labels == i]

        ith_cluster_silhouette_values.sort()

        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i

        color = cm.nipy_spectral(float(i) / n_clusters)
        ax1.fill_betweenx(
            np.arange(y_lower, y_upper),
            0,
            ith_cluster_silhouette_values,
            facecolor=color,
            edgecolor=color,
            alpha=0.7,
        )

        # Label the silhouette plots with their cluster numbers at the middle
        ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples

    #ax1.set_title("The silhouette plot for the various clusters.")
    ax1.set_xlabel("Silhouette Coefficient")
    ax1.set_ylabel("Cluster label")

    # The vertical line for average silhouette score of all the values
    ax1.axvline(x=silhouette_avg, color="red", linestyle="--")

    ax1.set_yticks([])  # Clear the yaxis labels / ticks
    ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

    # 2nd Plot showing the actual clusters formed
    colors = cm.nipy_spectral(cluster_labels.astype(float) / n_clusters)
    ax2.scatter(
        principal_components[:, 0], principal_components[:, 1], marker=".", s=30, lw=0, alpha=0.7, c=colors, edgecolor="k"
    )

    ax2.set_xlabel("PC1")
    ax2.set_ylabel("PC2")

    plt.suptitle(
        "k = %d"
        % n_clusters,
        fontsize=14,
        fontweight="bold",
    )

    plt.savefig("figures/silhouette_angle_k{}.png".format(n_clusters))
    plt.close()
    
plt.show()


# In[ ]:


range_n_clusters = [2,8,9,10,11,12]
X = theta

# Create a single figure for all silhouette plots
#fig, axs = plt.subplots(len(range_n_clusters), 2, figsize=(18, len(range_n_clusters) * 7))

for n_clusters in range_n_clusters:
    # Create a subplot with 1 row and 2 columns
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_size_inches(18, 7)

    # The 1st subplot is the silhouette plot
    # The silhouette coefficient can range from -1, 1 but in this example all
    # lie within [-0.1, 1]
    ax1.set_xlim([-0.1, 1])
    # The (n_clusters+1)*10 is for inserting blank space between silhouette
    # plots of individual clusters, to demarcate them clearly.
    ax1.set_ylim([0, len(X) + (n_clusters + 1) * 10])

    # Initialize the clusterer with n_clusters value and a random generator
    # seed of 10 for reproducibility.
    clusterer = KMeans(n_clusters=n_clusters, random_state=10)
    cluster_labels = clusterer.fit_predict(X)

    # The silhouette_score gives the average value for all the samples.
    # This gives a perspective into the density and separation of the formed
    # clusters
    silhouette_avg = silhouette_score(X, cluster_labels)
    print(
        "For n_clusters =",
        n_clusters,
        "The average silhouette_score is :",
        silhouette_avg,
    )

    # Compute the silhouette scores for each sample
    sample_silhouette_values = silhouette_samples(X, cluster_labels)

    y_lower = 10
    for i in range(n_clusters):
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = sample_silhouette_values[cluster_labels == i]

        ith_cluster_silhouette_values.sort()

        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i

        color = cm.nipy_spectral(float(i) / n_clusters)
        ax1.fill_betweenx(
            np.arange(y_lower, y_upper),
            0,
            ith_cluster_silhouette_values,
            facecolor=color,
            edgecolor=color,
            alpha=0.7,
        )

        # Label the silhouette plots with their cluster numbers at the middle
        ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples

    #ax1.set_title("The silhouette plot for the various clusters.")
    ax1.set_xlabel("Silhouette Coefficient")
    ax1.set_ylabel("Cluster label")

    # The vertical line for average silhouette score of all the values
    ax1.axvline(x=silhouette_avg, color="red", linestyle="--")

    ax1.set_yticks([])  # Clear the yaxis labels / ticks
    ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

    # 2nd Plot showing the actual clusters formed
    colors = cm.nipy_spectral(cluster_labels.astype(float) / n_clusters)
    ax2.scatter(
        principal_components[:, 0], principal_components[:, 1], marker=".", s=30, lw=0, alpha=0.7, c=colors, edgecolor="k"
    )

    ax2.set_xlabel("PC1")
    ax2.set_ylabel("PC2")

    plt.suptitle(
        "k = %d"
        % n_clusters,
        fontsize=14,
        fontweight="bold",
    )

    plt.savefig("figures/silhouette_angle_k{}.png".format(n_clusters))
    plt.close()
    
plt.show()


# ## Clustering for angle plot

# In[109]:


# Calculate the angle using arctan(PC2/PC1)
theta = np.arctan2(principal_components[:, 1], (principal_components[:, 0])+24.233875)

# Reshape theta to a column vector
theta = theta.reshape(-1, 1)

# Perform k-means clustering using theta 
kmeans = KMeans(n_clusters=7, random_state=0) 
cluster_labels = kmeans.fit_predict(theta)

# Count the occurrences of each cluster label
cluster_counts = [np.sum(cluster_labels == i) for i in range(7)]  

# Plot PC1 vs PC2, color-coding the data points based on the clusters
plt.figure(figsize=(8, 6))
for i in range(7): 
    plt.scatter(principal_components[cluster_labels == i, 0], principal_components[cluster_labels == i, 1], s=0.5, label=f'Cluster {i+1} ({cluster_counts[i]})')

plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.grid(True)

legend = plt.legend(loc='upper left',fontsize = 10)
for handle in legend.legendHandles:
    handle.set_sizes([80])  

plt.savefig('figures/pc_cluster')
plt.show()


# In[110]:


# Plot PC1 vs PC2
plt.figure(figsize=(8, 6))
for i in range(7):  
    plt.scatter(principal_components[cluster_labels == i, 0], principal_components[cluster_labels == i, 1], s=0.2, label=f'Cluster {i+1} ({cluster_counts[i]})')

plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.grid(True)

legend = plt.legend(loc='best',fontsize = 9)
for handle in legend.legendHandles:
    handle.set_sizes([80])  

plt.xlim(-25,20)
plt.ylim(-20,20)

plt.savefig('figures/pc_cluster_zoom')
plt.show()


# In[111]:


# Count the number of spectra in each cluster
cluster_counts = np.bincount(kmeans.labels_)

# Convert the counts to a dictionary
cluster_counts_dict = {f"Cluster {i + 1}": count for i, count in enumerate(cluster_counts)}

# Convert the dictionary into a DataFrame
df = pd.DataFrame.from_dict(cluster_counts_dict, orient='index', columns=['Number of Spectra'])

df


# In[112]:


print(len(kmeans.labels_))


# In[113]:


order = ['whitedwarf','F', 'G', 'K', 'SF', 'Pass', 'SN']


# In[114]:


# Count the number of each type of object in each cluster
cluster_percentage_counts_dict = {}

for cluster_id in np.unique(kmeans.labels_):
    cluster_indices = np.where(kmeans.labels_ == cluster_id)[0]
    cluster_origins = [data_type[i] for i in cluster_indices]
    counts = {obj: cluster_origins.count(obj) for obj in order}
    
    total_count = sum(counts.values())
    percentages_counts = {obj: (f"{(counts[obj] / total_count) * 100:.2f}% ({counts[obj]})") for obj in order}
    
    cluster_percentage_counts_dict[f"Cluster {cluster_id + 1}"] = percentages_counts

# Convert the dictionary into a DataFrame
df_percentage_counts = pd.DataFrame.from_dict(cluster_percentage_counts_dict, orient='index')

# Transpose the DataFrame
df_percentage_counts = df_percentage_counts.transpose()

df_percentage_counts


# In[115]:


print(df_percentage_counts.to_latex())


# In[116]:


# Count the total number of occurrences of each object type across all clusters
total_counts = {}
for obj in order:
    total_counts[obj] = sum(int(row.split('(')[1].split(')')[0]) for row in df_percentage_counts.loc[obj])

# Calculate the percentage of each type in each cluster based on the total count
cluster_percentage_counts_dict = {}
for cluster_id in np.unique(kmeans.labels_):
    percentages_counts = {}
    for obj in order:
        count = df_percentage_counts.at[obj, f'Cluster {cluster_id + 1}'].split('(')[1].split(')')[0]
        percentage = (int(count) / total_counts[obj]) * 100 if total_counts[obj] != 0 else 0
        percentages_counts[obj] = f"{percentage:.2f}% ({count})"
    
    cluster_percentage_counts_dict[f"Cluster {cluster_id + 1}"] = percentages_counts

# Convert the dictionary into a DataFrame
df_percentage = pd.DataFrame.from_dict(cluster_percentage_counts_dict, orient='index')

df_percentage = df_percentage.transpose()

df_percentage


# In[117]:


# Create DataFrame
data = {
    'Cluster 1': {'whitedwarf': '0.08', 'F': '3.78', 'G': '6.13', 'K': '31.87', 'SF': '7.03', 'Pass': '38.70', 'SN': '4.28'},
    'Cluster 2': {'whitedwarf': '23.58', 'F': '25.61', 'G': '0.08', 'K': '0.00', 'SF': '4.49', 'Pass': '0.01', 'SN': '30.51'},
    'Cluster 3': {'whitedwarf': '0.94', 'F': '2.11', 'G': '27.28', 'K': '5.60', 'SF': '38.24', 'Pass': '60.40', 'SN': '20.55'},
    'Cluster 4': {'whitedwarf': '4.12', 'F': '26.69', 'G': '24.91', 'K': '0.00', 'SF': '50.18', 'Pass': '0.77', 'SN': '18.78'},
    'Cluster 5': {'whitedwarf': '0.00', 'F': '38.81', 'G': '19.66', 'K': '26.30', 'SF': '0.07', 'Pass': '0.10', 'SN': '0.15'},
    'Cluster 6': {'whitedwarf': '71.27', 'F': '0.00', 'G': '0.02', 'K': '0.00', 'SF': '0.00', 'Pass': '0.01', 'SN': '25.72'},
    'Cluster 7': {'whitedwarf': '0.00', 'F': '3.00', 'G': '21.92', 'K': '36.23', 'SF': '0.00', 'Pass': '0.00', 'SN': '0.00'}
}

df = pd.DataFrame(data)

# Function to extract percentage from string
def float_conv(s):
    return float(s)

# Apply the function to convert strings to numeric values
df = df.applymap(float_conv)

#Create bar chart:

fig, ax = plt.subplots(figsize=(12, 8))
df.plot(kind='bar', stacked=True, ax=ax)

for col in df.columns:
    max_val = df[col].max()
    max_idx = df[col].idxmax()
    idx = df.index.get_loc(max_idx)
    cluster_num = col.split()[1]
    ax.text(idx, max_val, f"Cluster {cluster_num}", ha='center', va='bottom', color='black')

ax.set_ylabel('%')
#ax.set_title('Percentage Distribution of Object Types across Clusters')

plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()


# In[118]:


# Plotting
import seaborn as sns
# Create DataFrame
data = {
    'Cluster 1': {'whitedwarf': '0.08', 'F': '3.78', 'G': '6.13', 'K': '31.87', 'SF': '7.03', 'Pass': '38.70', 'SN': '4.28'},
    'Cluster 2': {'whitedwarf': '23.58', 'F': '25.61', 'G': '0.08', 'K': '0.00', 'SF': '4.49', 'Pass': '0.01', 'SN': '30.51'},
    'Cluster 3': {'whitedwarf': '0.94', 'F': '2.11', 'G': '27.28', 'K': '5.60', 'SF': '38.24', 'Pass': '60.40', 'SN': '20.55'},
    'Cluster 4': {'whitedwarf': '4.12', 'F': '26.69', 'G': '24.91', 'K': '0.00', 'SF': '50.18', 'Pass': '0.77', 'SN': '18.78'},
    'Cluster 5': {'whitedwarf': '0.00', 'F': '38.81', 'G': '19.66', 'K': '26.30', 'SF': '0.07', 'Pass': '0.10', 'SN': '0.15'},
    'Cluster 6': {'whitedwarf': '71.27', 'F': '0.00', 'G': '0.02', 'K': '0.00', 'SF': '0.00', 'Pass': '0.01', 'SN': '25.72'},
    'Cluster 7': {'whitedwarf': '0.00', 'F': '3.00', 'G': '21.92', 'K': '36.23', 'SF': '0.00', 'Pass': '0.00', 'SN': '0.00'}
}

df = pd.DataFrame(data)

# Function to extract percentage from string
def float_conv(s):
    return float(s)/100

# Apply the function to convert strings to numeric values
df = df.applymap(float_conv)
# Create heatmap
fig, ax = plt.subplots(figsize=(8, 8))
sns.heatmap(df.T, annot=True, cmap='Blues', fmt='.2f', linewidths=.5, ax=ax, cbar=False)

#ax.set_title('Percentage Distribution of Object Types across Clusters')
#ax.set_xlabel('Object Types')
#ax.set_ylabel('Clusters')
plt.savefig('figures/distributionheatmap')
plt.show()


# In[119]:


print(df_percentage_counts.to_latex())


# ## Galactic / Extragalactic - angle

# In[353]:


# Define the groups
star_groups = {'F', 'G', 'K', 'whitedwarf'}
galaxy_groups = {'Pass', 'SN', 'SF'}
n_clusters = 7
# Calculate total count of objects in each group within each cluster
total_counts_per_cluster = {}
for cluster_id in np.unique(kmeans.labels_):
    cluster_indices = np.where(kmeans.labels_ == cluster_id)[0]
    cluster_origins = [data_type[i] for i in cluster_indices]
    
    total_count_stars = sum(1 for obj in cluster_origins if obj in star_groups)
    total_count_galaxies = sum(1 for obj in cluster_origins if obj in galaxy_groups)
    
    total_counts_per_cluster[cluster_id] = {'Galactic': total_count_stars, 'Extragalactic': total_count_galaxies}

# List of rows for the table
table_data = []
for group in ['Galactic', 'Extragalactic']:
    row = [group]
    for cluster_id in range(n_clusters):
        total_count_in_group = total_counts_per_cluster[cluster_id][group]
        total_count_in_cluster = sum(total_counts_per_cluster[cluster_id].values())
        percentage = (total_count_in_group / total_count_in_cluster) * 100 if total_count_in_cluster != 0 else 0
        row.append(f"{percentage:.2f}% ({total_count_in_group})")
    table_data.append(row)

headers = ["Group"] + [f"Cluster {cluster_id + 1}" for cluster_id in range(n_clusters)]

df_group_counts = pd.DataFrame(table_data, columns=headers)

df_group_counts


# In[357]:


print(df_group_counts.to_latex())


# In[355]:


# Define the groups
star_groups = {'F', 'G', 'K', 'whitedwarf'}
galaxy_groups = {'Pass', 'SN', 'SF'}

# Calculate total count of objects in each group within each cluster
total_counts_per_cluster = {}
for cluster_id in np.unique(kmeans.labels_):
    cluster_indices = np.where(kmeans.labels_ == cluster_id)[0]
    cluster_origins = [data_type[i] for i in cluster_indices]
    
    total_count_stars = sum(1 for obj in cluster_origins if obj in star_groups)
    total_count_galaxies = sum(1 for obj in cluster_origins if obj in galaxy_groups)
    
    total_counts_per_cluster[cluster_id] = {'Galactic': total_count_stars, 'Extragalactic': total_count_galaxies}

# List of rows for the table
table_data = []
for group in ['Galactic', 'Extragalactic']:
    row = [group]
    for cluster_id in range(n_clusters):
        total_count_in_group = total_counts_per_cluster[cluster_id][group]
        total_count_in_cluster = sum(total_counts_per_cluster[cluster_id].values())
        percentage = (total_count_in_group / total_count_in_cluster) * 100 if total_count_in_cluster != 0 else 0
        row.append(f"{percentage:.2f}%")
    table_data.append(row)

headers = ["Group"] + [f"Cluster {cluster_id + 1}" for cluster_id in range(n_clusters)]

df_group_counts1 = pd.DataFrame(table_data, columns=headers)

#Creat bar chart:

# Convert percentage strings to float values
for col in df_group_counts1.columns[1:]:
    df_group_counts1[col] = df_group_counts1[col].str.rstrip('%').astype(float)
                                
# Set the index to 'Group' column for easier plotting
df_group_counts1.set_index('Group', inplace=True)

# Plot the bar chart
ax = df_group_counts1.T.plot(kind='bar', figsize=(10, 6))
#plt.title('Percentage of Groups in Each Cluster')
for idx, col in enumerate(df_group_counts1.columns):
    for i in ax.patches[idx::len(df_group_counts1.columns)]:
        ax.text(i.get_x() + i.get_width() / 2, i.get_height() + 0.5,
                f'{i.get_height()}', ha='center', va='bottom')
        
plt.xlabel('Clusters')
plt.ylabel('Percentage / %')
plt.xticks(rotation=0)
plt.legend(title='Object Type')
plt.tight_layout()
plt.savefig('figures/clusters_bar_angle')
plt.show()


# ## Elbow/ Dunn/ Silhouette - angle

# In[328]:


# Dunn index
def dunn_index(data, labels, centroids):
    min_intercluster_distance = np.inf
    max_intracluster_distance = -np.inf

    # Max intracluster distance
    for i in range(len(centroids)):
        intracluster_distance = np.max(np.linalg.norm(data[labels == i] - centroids[i], axis=1))
        if intracluster_distance > max_intracluster_distance:
            max_intracluster_distance = intracluster_distance

    # Min intercluster distance
    for i in range(len(centroids)):
        for j in range(i + 1, len(centroids)):
            intercluster_distance = np.linalg.norm(centroids[i] - centroids[j])
            if intercluster_distance < min_intercluster_distance:
                min_intercluster_distance = intercluster_distance

    dunn_index = min_intercluster_distance / max_intracluster_distance
    return dunn_index


# ELBOW METHOD AGAIN + DUNN INDEX:

# range of clusters
min_clusters = 1
max_clusters = 20

wcss = []
dunn_indices = []

# Iterate over the range of clusters
for n_clusters in range(min_clusters, max_clusters + 1):
    kmeans = KMeans(n_clusters=n_clusters, init='k-means++', random_state=42)
    kmeans.fit(theta)  
    wcss.append(kmeans.inertia_)  
    
    centroids = kmeans.cluster_centers_
    labels = kmeans.labels_
    dunn_index_value = dunn_index(theta, labels, centroids)
    dunn_indices.append(dunn_index_value) 

# Elbow curve
fig, ax1 = plt.subplots()
color = 'tab:red'
ax1.set_xlabel('Number of Clusters')
ax1.set_ylabel('WCSS', color=color)
ax1.plot(range(min_clusters, max_clusters + 1), wcss, marker='o', linestyle='-', color=color)
ax1.tick_params(axis='y', labelcolor=color)
#plt.title('Elbow Method')
plt.xticks(range(min_clusters, max_clusters + 1))
plt.grid(True)
plt.savefig('figures/Elbow_C_angle')

# Dunn index
fig, ax2 = plt.subplots()
color = 'tab:blue'
ax2.set_xlabel('Number of Clusters')
ax2.set_ylabel('Dunn Index', color=color)
ax2.plot(range(min_clusters, max_clusters + 1), dunn_indices, marker='o', linestyle='-', color=color)
ax2.tick_params(axis='y', labelcolor=color)
#plt.title('Dunn Index')
plt.xticks(range(min_clusters, max_clusters + 1))
plt.grid(True)
plt.savefig('figures/Dunn_C_angle')
plt.show()


# In[120]:


range_n_clusters = [3,4,5,6,7]
X = theta

# Create a single figure for all silhouette plots
#fig, axs = plt.subplots(len(range_n_clusters), 2, figsize=(18, len(range_n_clusters) * 7))

for n_clusters in range_n_clusters:
    # Create a subplot with 1 row and 2 columns
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_size_inches(18, 7)

    # The 1st subplot is the silhouette plot
    # The silhouette coefficient can range from -1, 1 but in this example all
    # lie within [-0.1, 1]
    ax1.set_xlim([-0.1, 1])
    # The (n_clusters+1)*10 is for inserting blank space between silhouette
    # plots of individual clusters, to demarcate them clearly.
    ax1.set_ylim([0, len(X) + (n_clusters + 1) * 10])

    # Initialize the clusterer with n_clusters value and a random generator
    # seed of 10 for reproducibility.
    clusterer = KMeans(n_clusters=n_clusters, random_state=10)
    cluster_labels = clusterer.fit_predict(X)

    # The silhouette_score gives the average value for all the samples.
    # This gives a perspective into the density and separation of the formed
    # clusters
    silhouette_avg = silhouette_score(X, cluster_labels)
    print(
        "For n_clusters =",
        n_clusters,
        "The average silhouette_score is :",
        silhouette_avg,
    )

    # Compute the silhouette scores for each sample
    sample_silhouette_values = silhouette_samples(X, cluster_labels)

    y_lower = 10
    for i in range(n_clusters):
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = sample_silhouette_values[cluster_labels == i]

        ith_cluster_silhouette_values.sort()

        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i

        color = cm.nipy_spectral(float(i) / n_clusters)
        ax1.fill_betweenx(
            np.arange(y_lower, y_upper),
            0,
            ith_cluster_silhouette_values,
            facecolor=color,
            edgecolor=color,
            alpha=0.7,
        )

        # Label the silhouette plots with their cluster numbers at the middle
        ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples

    #ax1.set_title("The silhouette plot for the various clusters.")
    ax1.set_xlabel("Silhouette Coefficient")
    ax1.set_ylabel("Cluster label")

    # The vertical line for average silhouette score of all the values
    ax1.axvline(x=silhouette_avg, color="red", linestyle="--")

    ax1.set_yticks([])  # Clear the yaxis labels / ticks
    ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

    # 2nd Plot showing the actual clusters formed
    colors = cm.nipy_spectral(cluster_labels.astype(float) / n_clusters)
    ax2.scatter(
        principal_components[:, 0], principal_components[:, 1], marker=".", s=30, lw=0, alpha=0.7, c=colors, edgecolor="k"
    )

    ax2.set_xlabel("PC1")
    ax2.set_ylabel("PC2")

    plt.suptitle(
        "k = %d"
        % n_clusters,
        fontsize=14,
        fontweight="bold",
    )

    plt.savefig("figures/silhouette_angle_k{}.png".format(n_clusters))
    plt.close()
    
plt.show()


# In[122]:


range_n_clusters = [2,8,9,10,11,12]
X = theta

# Create a single figure for all silhouette plots
#fig, axs = plt.subplots(len(range_n_clusters), 2, figsize=(18, len(range_n_clusters) * 7))

for n_clusters in range_n_clusters:
    # Create a subplot with 1 row and 2 columns
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_size_inches(18, 7)

    # The 1st subplot is the silhouette plot
    # The silhouette coefficient can range from -1, 1 but in this example all
    # lie within [-0.1, 1]
    ax1.set_xlim([-0.1, 1])
    # The (n_clusters+1)*10 is for inserting blank space between silhouette
    # plots of individual clusters, to demarcate them clearly.
    ax1.set_ylim([0, len(X) + (n_clusters + 1) * 10])

    # Initialize the clusterer with n_clusters value and a random generator
    # seed of 10 for reproducibility.
    clusterer = KMeans(n_clusters=n_clusters, random_state=10)
    cluster_labels = clusterer.fit_predict(X)

    # The silhouette_score gives the average value for all the samples.
    # This gives a perspective into the density and separation of the formed
    # clusters
    silhouette_avg = silhouette_score(X, cluster_labels)
    print(
        "For n_clusters =",
        n_clusters,
        "The average silhouette_score is :",
        silhouette_avg,
    )

    # Compute the silhouette scores for each sample
    sample_silhouette_values = silhouette_samples(X, cluster_labels)

    y_lower = 10
    for i in range(n_clusters):
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = sample_silhouette_values[cluster_labels == i]

        ith_cluster_silhouette_values.sort()

        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i

        color = cm.nipy_spectral(float(i) / n_clusters)
        ax1.fill_betweenx(
            np.arange(y_lower, y_upper),
            0,
            ith_cluster_silhouette_values,
            facecolor=color,
            edgecolor=color,
            alpha=0.7,
        )

        # Label the silhouette plots with their cluster numbers at the middle
        ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples

    #ax1.set_title("The silhouette plot for the various clusters.")
    ax1.set_xlabel("Silhouette Coefficient")
    ax1.set_ylabel("Cluster label")

    # The vertical line for average silhouette score of all the values
    ax1.axvline(x=silhouette_avg, color="red", linestyle="--")

    ax1.set_yticks([])  # Clear the yaxis labels / ticks
    ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

    # 2nd Plot showing the actual clusters formed
    colors = cm.nipy_spectral(cluster_labels.astype(float) / n_clusters)
    ax2.scatter(
        principal_components[:, 0], principal_components[:, 1], marker=".", s=30, lw=0, alpha=0.7, c=colors, edgecolor="k"
    )

    ax2.set_xlabel("PC1")
    ax2.set_ylabel("PC2")

    plt.suptitle(
        "k = %d"
        % n_clusters,
        fontsize=14,
        fontweight="bold",
    )

    plt.savefig("figures/silhouette_angle_k{}.png".format(n_clusters))
    plt.close()
    
plt.show()


# # Using PC1 and PC2 as classification parameters

# ## Elbow Method

# In[17]:


min_clusters = 1
max_clusters = 20  

wcss = []

# Iterate over the range of clusters
for n_clusters in range(min_clusters, max_clusters + 1):
    kmeans = KMeans(n_clusters=n_clusters, init='k-means++', random_state=42)
    kmeans.fit(principal_components[:, :2]) 
    wcss.append(kmeans.inertia_) 

# Plot the elbow curve
plt.plot(range(min_clusters, max_clusters + 1), wcss, marker='o', linestyle='-')
plt.xlabel('Number of Clusters')
plt.ylabel('Within-Cluster Sum of Squares (WCSS)')
plt.title('Elbow Method')
plt.xticks(range(min_clusters, max_clusters + 1))
plt.grid(True)
plt.show()


# ## Dunn Index

# In[18]:


# Dunn index
def dunn_index(data, labels, centroids):
    min_intercluster_distance = np.inf
    max_intracluster_distance = -np.inf

    # Max intracluster distance
    for i in range(len(centroids)):
        intracluster_distance = np.max(np.linalg.norm(data[labels == i] - centroids[i], axis=1))
        if intracluster_distance > max_intracluster_distance:
            max_intracluster_distance = intracluster_distance

    # Min intercluster distance
    for i in range(len(centroids)):
        for j in range(i + 1, len(centroids)):
            intercluster_distance = np.linalg.norm(centroids[i] - centroids[j])
            if intercluster_distance < min_intercluster_distance:
                min_intercluster_distance = intercluster_distance

    dunn_index = min_intercluster_distance / max_intracluster_distance
    return dunn_index


# ELBOW METHOD AGAIN + DUNN INDEX:

# range of clusters
min_clusters = 1
max_clusters = 20

wcss = []
dunn_indices = []

# Iterate over the range of clusters
for n_clusters in range(min_clusters, max_clusters + 1):
    kmeans = KMeans(n_clusters=n_clusters, init='k-means++', random_state=42)
    kmeans.fit(principal_components[:, :2])  
    wcss.append(kmeans.inertia_)  
    
    centroids = kmeans.cluster_centers_
    labels = kmeans.labels_
    dunn_index_value = dunn_index(principal_components[:, :2], labels, centroids)
    dunn_indices.append(dunn_index_value) 

# Elbow curve
fig, ax1 = plt.subplots()
color = 'tab:red'
ax1.set_xlabel('Number of Clusters')
ax1.set_ylabel('WCSS', color=color)
ax1.plot(range(min_clusters, max_clusters + 1), wcss, marker='o', linestyle='-', color=color)
ax1.tick_params(axis='y', labelcolor=color)
#plt.title('Elbow Method')
plt.xticks(range(min_clusters, max_clusters + 1))
plt.grid(True)
plt.savefig('figures/Elbow_C')

# Dunn index
fig, ax2 = plt.subplots()
color = 'tab:blue'
ax2.set_xlabel('Number of Clusters')
ax2.set_ylabel('Dunn Index', color=color)
ax2.plot(range(min_clusters, max_clusters + 1), dunn_indices, marker='o', linestyle='-', color=color)
ax2.tick_params(axis='y', labelcolor=color)
#plt.title('Dunn Index')
plt.xticks(range(min_clusters, max_clusters + 1))
plt.grid(True)
plt.savefig('figures/Dunn_C')
plt.show()


# ## Silhouette Score

# In[19]:


range_n_clusters = [3,4,5,6]
X = principal_components[:, :2]

# Create a single figure for all silhouette plots
#fig, axs = plt.subplots(len(range_n_clusters), 2, figsize=(18, len(range_n_clusters) * 7))

for n_clusters in range_n_clusters:
    # Create a subplot with 1 row and 2 columns
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_size_inches(18, 7)

    # The 1st subplot is the silhouette plot
    # The silhouette coefficient can range from -1, 1 but in this example all
    # lie within [-0.1, 1]
    ax1.set_xlim([-0.1, 1])
    # The (n_clusters+1)*10 is for inserting blank space between silhouette
    # plots of individual clusters, to demarcate them clearly.
    ax1.set_ylim([0, len(X) + (n_clusters + 1) * 10])

    # Initialize the clusterer with n_clusters value and a random generator
    # seed of 10 for reproducibility.
    clusterer = KMeans(n_clusters=n_clusters, random_state=10)
    cluster_labels = clusterer.fit_predict(X)

    # The silhouette_score gives the average value for all the samples.
    # This gives a perspective into the density and separation of the formed
    # clusters
    silhouette_avg = silhouette_score(X, cluster_labels)
    print(
        "For n_clusters =",
        n_clusters,
        "The average silhouette_score is :",
        silhouette_avg,
    )

    # Compute the silhouette scores for each sample
    sample_silhouette_values = silhouette_samples(X, cluster_labels)

    y_lower = 10
    for i in range(n_clusters):
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = sample_silhouette_values[cluster_labels == i]

        ith_cluster_silhouette_values.sort()

        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i

        color = cm.nipy_spectral(float(i) / n_clusters)
        ax1.fill_betweenx(
            np.arange(y_lower, y_upper),
            0,
            ith_cluster_silhouette_values,
            facecolor=color,
            edgecolor=color,
            alpha=0.7,
        )

        # Label the silhouette plots with their cluster numbers at the middle
        ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples

    #ax1.set_title("The silhouette plot for the various clusters.")
    ax1.set_xlabel("Silhouette Coefficient")
    ax1.set_ylabel("Cluster label")

    # The vertical line for average silhouette score of all the values
    ax1.axvline(x=silhouette_avg, color="red", linestyle="--")

    ax1.set_yticks([])  # Clear the yaxis labels / ticks
    ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

    # 2nd Plot showing the actual clusters formed
    colors = cm.nipy_spectral(cluster_labels.astype(float) / n_clusters)
    ax2.scatter(
        X[:, 0], X[:, 1], marker=".", s=30, lw=0, alpha=0.7, c=colors, edgecolor="k"
    )

    # Labeling the clusters
    centers = clusterer.cluster_centers_
    # Draw white circles at cluster centers
    ax2.scatter(
        centers[:, 0],
        centers[:, 1],
        marker="o",
        c="white",
        alpha=1,
        s=200,
        edgecolor="k",
    )

    for i, c in enumerate(centers):
        ax2.scatter(c[0], c[1], marker="$%d$" % i, alpha=1, s=50, edgecolor="k")

    #ax2.set_title("The visualization of the clustered data.")
    ax2.set_xlabel("PC1")
    ax2.set_ylabel("PC2")

    plt.suptitle(
        "k = %d"
        % n_clusters,
        fontsize=14,
        fontweight="bold",
    )

    plt.savefig("figures/silhouette_k{}.png".format(n_clusters))
    plt.close()
    
plt.show()


# In[20]:


range_n_clusters = [3, 4, 5, 6]
silhouette_scores = []

for n_clusters in range_n_clusters:
    clusterer = KMeans(n_clusters=n_clusters, random_state=10)
    cluster_labels = clusterer.fit_predict(X)
    silhouette_avg = silhouette_score(X, cluster_labels)
    silhouette_scores.append(round(silhouette_avg, 2))
    
    
# Create a DataFrame to store silhouette scores
df_silhouette = pd.DataFrame({
    'Cluster': range_n_clusters,
    'Silhouette Score': silhouette_scores
})

df_silhouette


# In[21]:


print(df_silhouette.to_latex())


# ## Clusters

# In[52]:


# Reorder data_type
order = ['whitedwarf','F', 'G', 'K', 'SF', 'Pass', 'SN']


# In[53]:


# Perform k-means clustering on the PCA-transformed data
n_clusters = 4  # Number of clusters
kmeans = KMeans(n_clusters=n_clusters)
kmeans.fit(principal_components[:, :2])

# Count the number of data points in each cluster
cluster_counts = np.bincount(kmeans.labels_)

# Convert the counts to a dictionary
cluster_counts_dict = {f"Cluster {i + 1}": count for i, count in enumerate(cluster_counts)}

# Convert the dictionary into a DataFrame
df = pd.DataFrame.from_dict(cluster_counts_dict, orient='index', columns=['Number of Spectra'])

df


# In[54]:


print(df.to_latex())


# In[55]:


# Count the number of each type of object in each cluster
for cluster_id in np.unique(kmeans.labels_):
    cluster_indices = np.where(kmeans.labels_ == cluster_id)[0]
    cluster_origins = [data_type[i] for i in cluster_indices]
    counts = {obj: cluster_origins.count(obj) for obj in order}
    cluster_counts_dict[f"Cluster {cluster_id + 1}"] = counts

# Convert the dictionary into a DataFrame
df = pd.DataFrame.from_dict(cluster_counts_dict, orient='index')

df


# In[56]:


# Count the number of each type of object in each cluster
cluster_percentage_counts_dict = {}

for cluster_id in np.unique(kmeans.labels_):
    cluster_indices = np.where(kmeans.labels_ == cluster_id)[0]
    cluster_origins = [data_type[i] for i in cluster_indices]
    counts = {obj: cluster_origins.count(obj) for obj in order}
    
    total_count = sum(counts.values())
    percentages_counts = {obj: (f"{(counts[obj] / total_count) * 100:.2f}% ({counts[obj]})") for obj in order}
    
    cluster_percentage_counts_dict[f"Cluster {cluster_id + 1}"] = percentages_counts

# Convert the dictionary into a DataFrame
df_percentage_counts = pd.DataFrame.from_dict(cluster_percentage_counts_dict, orient='index')

df_percentage_counts = df_percentage_counts.transpose()

df_percentage_counts


# In[57]:


print(df_percentage_counts.to_latex())


# In[58]:


# Define the groups
star_groups = {'F', 'G', 'K', 'whitedwarf'}
galaxy_groups = {'Pass', 'SN', 'SF'}

# Calculate total count of objects in each group within each cluster
total_counts_per_cluster = {}
for cluster_id in np.unique(kmeans.labels_):
    cluster_indices = np.where(kmeans.labels_ == cluster_id)[0]
    cluster_origins = [data_type[i] for i in cluster_indices]
    
    total_count_stars = sum(1 for obj in cluster_origins if obj in star_groups)
    total_count_galaxies = sum(1 for obj in cluster_origins if obj in galaxy_groups)
    
    total_counts_per_cluster[cluster_id] = {'Galactic': total_count_stars, 'Extragalactic': total_count_galaxies}

# List of rows for the table
table_data = []
for group in ['Galactic', 'Extragalactic']:
    row = [group]
    for cluster_id in range(n_clusters):
        total_count_in_group = total_counts_per_cluster[cluster_id][group]
        total_count_in_cluster = sum(total_counts_per_cluster[cluster_id].values())
        percentage = (total_count_in_group / total_count_in_cluster) * 100 if total_count_in_cluster != 0 else 0
        row.append(f"{percentage:.2f}% ({total_count_in_group})")
    table_data.append(row)

headers = ["Group"] + [f"Cluster {cluster_id + 1}" for cluster_id in range(n_clusters)]

df_group_counts = pd.DataFrame(table_data, columns=headers)

df_group_counts


# In[59]:


print(df_group_counts.to_latex())


# In[60]:


# Define the groups
star_groups = {'F', 'G', 'K', 'whitedwarf'}
galaxy_groups = {'Pass', 'SN', 'SF'}

# Calculate total count of objects in each group within each cluster
total_counts_per_cluster = {}
for cluster_id in np.unique(kmeans.labels_):
    cluster_indices = np.where(kmeans.labels_ == cluster_id)[0]
    cluster_origins = [data_type[i] for i in cluster_indices]
    
    total_count_stars = sum(1 for obj in cluster_origins if obj in star_groups)
    total_count_galaxies = sum(1 for obj in cluster_origins if obj in galaxy_groups)
    
    total_counts_per_cluster[cluster_id] = {'Galactic': total_count_stars, 'Extragalactic': total_count_galaxies}

# List of rows for the table
table_data = []
for group in ['Galactic', 'Extragalactic']:
    row = [group]
    for cluster_id in range(n_clusters):
        total_count_in_group = total_counts_per_cluster[cluster_id][group]
        total_count_in_cluster = sum(total_counts_per_cluster[cluster_id].values())
        percentage = (total_count_in_group / total_count_in_cluster) * 100 if total_count_in_cluster != 0 else 0
        row.append(f"{percentage:.2f}%")
    table_data.append(row)

headers = ["Group"] + [f"Cluster {cluster_id + 1}" for cluster_id in range(n_clusters)]

df_group_counts1 = pd.DataFrame(table_data, columns=headers)



# Convert percentage strings to float values
for col in df_group_counts1.columns[1:]:
    df_group_counts1[col] = df_group_counts1[col].str.rstrip('%').astype(float)
                                
# Set the index to 'Group' column for easier plotting
df_group_counts1.set_index('Group', inplace=True)

# Plot the bar chart
ax = df_group_counts1.T.plot(kind='bar', figsize=(10, 6))
#plt.title('Percentage of Groups in Each Cluster')
for idx, col in enumerate(df_group_counts1.columns):
    for i in ax.patches[idx::len(df_group_counts1.columns)]:
        ax.text(i.get_x() + i.get_width() / 2, i.get_height() + 0.5,
                f'{i.get_height()}', ha='center', va='bottom')
        
plt.xlabel('Clusters')
plt.ylabel('Percentage / %')
plt.xticks(rotation=0)
plt.legend(title='Object Type')
plt.tight_layout()
plt.savefig('figures/clusters_bar')
plt.show()


# In[61]:


# Plot 
plt.figure(figsize=(8, 6))
scatter = plt.scatter(principal_components[:, 0], principal_components[:, 1], c=kmeans.labels_, cmap='Set2', alpha=0.8, s=4)

# Cluster centres
for i, center in enumerate(kmeans.cluster_centers_):
    plt.scatter(center[0], center[1], marker='o', s=200, facecolors='none', edgecolors='black')
    plt.text(center[0], center[1], str(i + 1), fontsize=12, color='black', ha='center', va='center')

plt.xlabel('PC1')
plt.ylabel('PC2')

# Create custom legend labels with cluster labels and data point counts
legend_labels = [f"Cluster {i+1} ({count})" for i, count in enumerate(cluster_counts)]

# Plot legend with custom labels
plt.legend(handles=scatter.legend_elements()[0], labels=legend_labels)
plt.grid(True)
plt.savefig('figures/clusters')
plt.show()


# ## Subcluster

# In[62]:


# Initial K-means clustering
#n_clusters = 7  
#kmeans = KMeans(n_clusters=n_clusters)
#kmeans.fit(principal_components[:, :2])

# Identify the large cluster
large_cluster_idx = np.argmax(np.bincount(kmeans.labels_))

# Extract data points belonging to the large cluster
large_cluster_mask = kmeans.labels_ == large_cluster_idx
large_cluster_data = principal_components[large_cluster_mask]

# Elbow method
wcss = []
max_subclusters = 20  
for i in range(1, max_subclusters + 1):
    subkmeans = KMeans(n_clusters=i)
    subkmeans.fit(large_cluster_data)
    wcss.append(subkmeans.inertia_)

# Plot Elbow
fig, ax1 = plt.subplots()
color = 'tab:red'
ax1.set_xlabel('Number of Subclusters')
ax1.set_ylabel('WCSS', color=color)
ax1.plot(range(1, max_subclusters + 1), wcss, marker='o', linestyle='-', color=color)
ax1.tick_params(axis='y', labelcolor=color)
#plt.title('Elbow Method')
plt.grid(True)
plt.xticks(range(1, max_subclusters + 1))
plt.savefig('figures/Elbow_SC')
# Dunn index 
dunn_indices = []
for i in range(1, max_subclusters + 1):
    subkmeans = KMeans(n_clusters=i)
    subkmeans.fit(large_cluster_data)
    centroids = subkmeans.cluster_centers_
    labels = subkmeans.labels_
    dunn_index_value = dunn_index(large_cluster_data, labels, centroids)
    dunn_indices.append(dunn_index_value)

# Plot Dunn index
fig, ax2 = plt.subplots()
color = 'tab:blue'
ax2.set_xlabel('Number of Subclusters')
ax2.set_ylabel('Dunn Index', color=color)
ax2.plot(range(1, max_subclusters + 1), dunn_indices, marker='o', linestyle='-', color=color)
ax2.tick_params(axis='y', labelcolor=color)
#plt.title('Dunn Index')
plt.xticks(range(1, max_subclusters + 1))
plt.grid(True)
plt.savefig('figures/Dunn_SC')
plt.show()


# In[33]:


range_n_subclusters = [3,4,5,6,7,8]
X = large_cluster_data

for n_clusters in range_n_subclusters:
    # Create a subplot with 1 row and 2 columns
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_size_inches(18, 7)

    # The 1st subplot is the silhouette plot
    # The silhouette coefficient can range from -1, 1 but in this example all
    # lie within [-0.1, 1]
    ax1.set_xlim([-0.1, 1])
    # The (n_clusters+1)*10 is for inserting blank space between silhouette
    # plots of individual clusters, to demarcate them clearly.
    ax1.set_ylim([0, len(X) + (n_clusters + 1) * 10])

    # Initialize the clusterer with n_clusters value and a random generator
    # seed of 10 for reproducibility.
    clusterer = KMeans(n_clusters=n_clusters, random_state=10)
    cluster_labels = clusterer.fit_predict(X)

    # The silhouette_score gives the average value for all the samples.
    # This gives a perspective into the density and separation of the formed
    # clusters
    silhouette_avg = silhouette_score(X, cluster_labels)
    print(
        "For n_clusters =",
        n_clusters,
        "The average silhouette_score is :",
        silhouette_avg,
    )

    # Compute the silhouette scores for each sample
    sample_silhouette_values = silhouette_samples(X, cluster_labels)

    y_lower = 10
    for i in range(n_clusters):
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = sample_silhouette_values[cluster_labels == i]

        ith_cluster_silhouette_values.sort()

        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i

        color = cm.nipy_spectral(float(i) / n_clusters)
        ax1.fill_betweenx(
            np.arange(y_lower, y_upper),
            0,
            ith_cluster_silhouette_values,
            facecolor=color,
            edgecolor=color,
            alpha=0.7,
        )

        # Label the silhouette plots with their cluster numbers at the middle
        ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples

    #ax1.set_title("The silhouette plot for the various clusters.")
    ax1.set_xlabel("Silhouette Coefficient")
    ax1.set_ylabel("Cluster label")

    # The vertical line for average silhouette score of all the values
    ax1.axvline(x=silhouette_avg, color="red", linestyle="--")

    ax1.set_yticks([])  # Clear the yaxis labels / ticks
    ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

    # 2nd Plot showing the actual clusters formed
    colors = cm.nipy_spectral(cluster_labels.astype(float) / n_clusters)
    ax2.scatter(
        X[:, 0], X[:, 1], marker=".", s=30, lw=0, alpha=0.7, c=colors, edgecolor="k"
    )

    # Labeling the clusters
    centers = clusterer.cluster_centers_
    # Draw white circles at cluster centers
    ax2.scatter(
        centers[:, 0],
        centers[:, 1],
        marker="o",
        c="white",
        alpha=1,
        s=200,
        edgecolor="k",
    )

    for i, c in enumerate(centers):
        ax2.scatter(c[0], c[1], marker="$%d$" % i, alpha=1, s=50, edgecolor="k")

    #ax2.set_title("The visualization of the clustered data.")
    ax2.set_xlabel("PC1")
    ax2.set_ylabel("PC2")

    plt.suptitle(
        r"$k_{s} = %d$"
        % n_clusters,
        fontsize=14,
        fontweight="bold",
    )

    plt.savefig("figures/silhouette_sub_k{}.png".format(n_clusters))
    plt.close()
    
#plt.savefig("figures/silhouette_sc")
plt.show()


# In[34]:


range_n_subclusters = [3,4,5,6,7,8]
silhouette_scores = []

for n_clusters in range_n_subclusters:
    clusterer = KMeans(n_clusters=n_clusters, random_state=10)
    cluster_labels = clusterer.fit_predict(X)
    silhouette_avg = silhouette_score(X, cluster_labels)
    silhouette_scores.append(round(silhouette_avg, 2))
    
    
# Create a DataFrame to store silhouette scores
df_silhouette2 = pd.DataFrame({
    'SubCluster': range_n_subclusters,
    'Silhouette Score': silhouette_scores
})

df_silhouette2


# In[35]:


print(df_silhouette2.to_latex())


# ## k_s = 4

# In[63]:


# Initial K-means clustering
#n_clusters = 7  
#kmeans = KMeans(n_clusters=n_clusters)
#kmeans.fit(principal_components[:, :2])

# Identify the large cluster
#large_cluster_idx = np.argmax(np.bincount(kmeans.labels_))

# Extract data points belonging to the large cluster
#large_cluster_mask = kmeans.labels_ == large_cluster_idx
#large_cluster_data = principal_components[large_cluster_mask]

optimal_subclusters = 4 

# K-means on large cluster
subkmeans = KMeans(n_clusters=optimal_subclusters)
subkmeans.fit(large_cluster_data)

# Count the number of data points in each cluster
subcluster_counts = np.bincount(subkmeans.labels_)

# Convert the counts to a dictionary
subcluster_counts_dict = {f"Subcluster {i + 1}": count for i, count in enumerate(subcluster_counts)}

# Convert the dictionary into a DataFrame
df = pd.DataFrame.from_dict(subcluster_counts_dict, orient='index', columns=[f"Number of Spectra: Cluster {large_cluster_idx+1} ({np.sum(large_cluster_mask)})"])

df


# In[64]:


print(df.to_latex())


# In[65]:


# Count the number of each type of object in each cluster
subcluster_percentage_counts_dict = {}

for subcluster_id in np.unique(subkmeans.labels_):
    subcluster_indices = np.where(subkmeans.labels_ == subcluster_id)[0]
    subcluster_origins = [data_type[i] for i in subcluster_indices]
    counts = {obj: subcluster_origins.count(obj) for obj in order}
    
    total_count = sum(counts.values())
    percentages_counts = {obj: (f"{(counts[obj] / total_count) * 100:.2f}% ({counts[obj]})") for obj in order}
    
    subcluster_percentage_counts_dict[f"Cluster {subcluster_id + 1}"] = percentages_counts

# Convert the dictionary into a DataFrame
df_percentage = pd.DataFrame.from_dict(subcluster_percentage_counts_dict, orient='index')

# Transpose the DataFrame
df_percentage = df_percentage.transpose()

df_percentage


# In[39]:


print(df_percentage.to_latex())


# In[40]:


# Calculate total count of objects in each group within each cluster
total_counts_per_subcluster = {}
for subcluster_id in np.unique(subkmeans.labels_):
    subcluster_indices = np.where(subkmeans.labels_ == subcluster_id)[0]
    subcluster_origins = [data_type[i] for i in subcluster_indices]
    
    total_count_stars = sum(1 for obj in subcluster_origins if obj in star_groups)
    total_count_galaxies = sum(1 for obj in subcluster_origins if obj in galaxy_groups)
    
    total_counts_per_subcluster[subcluster_id] = {'Galactic': total_count_stars, 'Extragalactic': total_count_galaxies}

# List of rows for the table
table_data = []
for group in ['Galactic', 'Extragalactic']:
    row = [group]
    for subcluster_id in range(optimal_subclusters):
        total_count_in_group = total_counts_per_subcluster[subcluster_id][group]
        total_count_in_subcluster = sum(total_counts_per_subcluster[subcluster_id].values())
        percentage = (total_count_in_group / total_count_in_subcluster) * 100 if total_count_in_subcluster != 0 else 0
        row.append(f"{percentage:.2f}% ({total_count_in_group})")
    table_data.append(row)

headers = ["Group"] + [f"Cluster {subcluster_id + 1}" for subcluster_id in range(optimal_subclusters)]

df_group = pd.DataFrame(table_data, columns=headers)

df_group


# In[41]:


print(df_group.to_latex())


# In[42]:


# Calculate total count of objects in each group within each cluster
total_counts_per_subcluster = {}
for subcluster_id in np.unique(subkmeans.labels_):
    subcluster_indices = np.where(subkmeans.labels_ == subcluster_id)[0]
    subcluster_origins = [data_type[i] for i in subcluster_indices]
    
    total_count_stars = sum(1 for obj in subcluster_origins if obj in star_groups)
    total_count_galaxies = sum(1 for obj in subcluster_origins if obj in galaxy_groups)
    
    total_counts_per_subcluster[subcluster_id] = {'Galactic': total_count_stars, 'Extragalactic': total_count_galaxies}

# List of rows for the table
table_data = []
for group in ['Galactic', 'Extragalactic']:
    row = [group]
    for subcluster_id in range(optimal_subclusters):
        total_count_in_group = total_counts_per_subcluster[subcluster_id][group]
        total_count_in_subcluster = sum(total_counts_per_subcluster[subcluster_id].values())
        percentage = (total_count_in_group / total_count_in_subcluster) * 100 if total_count_in_subcluster != 0 else 0
        row.append(f"{percentage:.2f}%")
    table_data.append(row)

headers = ["Group"] + [f"Cluster {subcluster_id + 1}" for subcluster_id in range(optimal_subclusters)]

df_group1 = pd.DataFrame(table_data, columns=headers)

# Convert percentage strings to float values
for col in df_group1.columns[1:]:
    df_group1[col] = df_group1[col].str.rstrip('%').astype(float)
                                
# Set the index to 'Group' column for easier plotting
df_group1.set_index('Group', inplace=True)

# Plot the bar chart
ax = df_group1.T.plot(kind='bar', figsize=(10, 6))
#plt.title('Percentage of Groups in Each Cluster')
# Iterate over each column and annotate the bars with percentages
for idx, col in enumerate(df_group1.columns):
    for i in ax.patches[idx::len(df_group1.columns)]:
        ax.text(i.get_x() + i.get_width() / 2, i.get_height() + 0.5,
                f'{i.get_height()}', ha='center', va='bottom')
        
plt.xlabel('Clusters')
plt.ylabel('Percentage / %')
plt.xticks(rotation=0)
plt.legend(title='Object Type')
plt.tight_layout()
plt.savefig('figures/subclusters_bar_4')
plt.show()


# In[43]:


# Plot the subclusters
plt.figure(figsize=(8, 6))
for cluster_num in range(optimal_subclusters):
    cluster_data = large_cluster_data[subkmeans.labels_ == cluster_num]
    plt.scatter(cluster_data[:, 0], cluster_data[:, 1], label=f'SC{cluster_num+1} ({len(cluster_data)})',s=0.5)

# Cluster centers
for i, center in enumerate(subkmeans.cluster_centers_):
    plt.scatter(center[0], center[1], marker='o', s=200, facecolors='none', edgecolors='black')
    plt.text(center[0], center[1], str(i+1), fontsize=12, color='black', ha='center', va='center')
#plt.scatter(subkmeans.cluster_centers_[:, 0], subkmeans.cluster_centers_[:, 1], marker='x', color='black', label='Centroids')
 #plt.text(center[0], center[1], str(i+1), fontsize=12, color='black', ha='center', va='center')
plt.title(f'Cluster {large_cluster_idx+1}')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.legend()
plt.savefig('figures/Subclusters_4')
plt.show()


# In[78]:


order = ['whitedwarf','F', 'G', 'K', 'SF', 'Pass', 'SN']

n_clusters = 4  # Number of clusters
kmeans = KMeans(n_clusters=n_clusters)
kmeans.fit(principal_components[:, :2])

# Count the number of spectra in each cluster
cluster_counts = np.bincount(kmeans.labels_)

# Convert the counts to a dictionary
cluster_counts_dict = {f"Cluster {i + 1}": count for i, count in enumerate(cluster_counts)}

# Count the number of each type of object in each cluster
cluster_percentage_counts_dict = {}

for cluster_id in np.unique(kmeans.labels_):
    cluster_indices = np.where(kmeans.labels_ == cluster_id)[0]
    cluster_origins = [data_type[i] for i in cluster_indices]
    counts = {obj: cluster_origins.count(obj) for obj in order}
    
    total_count = sum(counts.values())
    percentages_counts = {obj: (f"{(counts[obj] / total_count) * 100:.2f}% ({counts[obj]})") for obj in order}
    
    cluster_percentage_counts_dict[f"Cluster {cluster_id + 1}"] = percentages_counts

# Identify the large cluster
large_cluster_idx = np.argmax(np.bincount(kmeans.labels_))

# Extract spectra belonging to the large cluster
large_cluster_mask = kmeans.labels_ == large_cluster_idx
large_cluster_data = principal_components[:, :2][large_cluster_mask]

optimal_subclusters = 4 

# K-means on large cluster
subkmeans = KMeans(n_clusters=optimal_subclusters)
subkmeans.fit(large_cluster_data)

# Count the number of spectra in each cluster
subcluster_counts = np.bincount(subkmeans.labels_)

# Convert the counts to a dictionary
subcluster_counts_dict = {f"Subcluster {i + 1}": count for i, count in enumerate(subcluster_counts)}

# Count the number of each type of object in each cluster
subcluster_percentage_counts_dict = {}

for subcluster_id in np.unique(subkmeans.labels_):
    subcluster_indices = np.where(subkmeans.labels_ == subcluster_id)[0]
    subcluster_origins = [data_type[i] for i in subcluster_indices]
    counts = {obj: subcluster_origins.count(obj) for obj in order}
    
    total_count = sum(counts.values())
    percentages_counts = {obj: (f"{(counts[obj] / total_count) * 100:.2f}% ({counts[obj]})") for obj in order}
    
    subcluster_percentage_counts_dict[f"Subcluster {subcluster_id + 1}"] = percentages_counts
    
# Convert the dictionary into a DataFrame
df_percentage = pd.DataFrame.from_dict(subcluster_percentage_counts_dict, orient='index')

# Transpose the DataFrame
df_percentage = df_percentage.transpose()

df_percentage


# ## k_s = 7

# In[44]:


# Initial K-means clustering
#n_clusters = 7  
#kmeans = KMeans(n_clusters=n_clusters)
#kmeans.fit(principal_components[:, :2])

# Identify the large cluster
large_cluster_idx = np.argmax(np.bincount(kmeans.labels_))

# Extract data points belonging to the large cluster
large_cluster_mask = kmeans.labels_ == large_cluster_idx
large_cluster_data = principal_components[large_cluster_mask]

optimal_subclusters = 7 

# K-means on large cluster
subkmeans = KMeans(n_clusters=optimal_subclusters)
subkmeans.fit(large_cluster_data)

# Count the number of data points in each cluster
subcluster_counts = np.bincount(subkmeans.labels_)

# Convert the counts to a dictionary
subcluster_counts_dict = {f"Subcluster {i + 1}": count for i, count in enumerate(subcluster_counts)}

# Convert the dictionary into a DataFrame
df = pd.DataFrame.from_dict(subcluster_counts_dict, orient='index', columns=[f"Number of Spectra: Cluster {large_cluster_idx+1} ({np.sum(large_cluster_mask)})"])

df


# In[45]:


print(df.to_latex())


# In[46]:


# Count the number of each type of object in each cluster
subcluster_percentage_counts_dict = {}

for subcluster_id in np.unique(subkmeans.labels_):
    subcluster_indices = np.where(subkmeans.labels_ == subcluster_id)[0]
    subcluster_origins = [data_type[i] for i in subcluster_indices]
    counts = {obj: subcluster_origins.count(obj) for obj in order}
    
    total_count = sum(counts.values())
    percentages_counts = {obj: (f"{(counts[obj] / total_count) * 100:.2f}% ({counts[obj]})") for obj in order}
    
    subcluster_percentage_counts_dict[f"Subcluster {subcluster_id + 1}"] = percentages_counts

# Convert the dictionary into a DataFrame
df_percentage = pd.DataFrame.from_dict(subcluster_percentage_counts_dict, orient='index')

# Transpose the DataFrame
df_percentage = df_percentage.transpose()

df_percentage


# In[47]:


print(df_percentage.to_latex())


# In[48]:


# Calculate total count of objects in each group within each cluster
total_counts_per_subcluster = {}
for subcluster_id in np.unique(subkmeans.labels_):
    subcluster_indices = np.where(subkmeans.labels_ == subcluster_id)[0]
    subcluster_origins = [data_type[i] for i in subcluster_indices]
    
    total_count_stars = sum(1 for obj in subcluster_origins if obj in star_groups)
    total_count_galaxies = sum(1 for obj in subcluster_origins if obj in galaxy_groups)
    
    total_counts_per_subcluster[subcluster_id] = {'Galactic': total_count_stars, 'Extragalactic': total_count_galaxies}

# List of rows for the table
table_data = []
for group in ['Galactic', 'Extragalactic']:
    row = [group]
    for subcluster_id in range(optimal_subclusters):
        total_count_in_group = total_counts_per_subcluster[subcluster_id][group]
        total_count_in_subcluster = sum(total_counts_per_subcluster[subcluster_id].values())
        percentage = (total_count_in_group / total_count_in_subcluster) * 100 if total_count_in_subcluster != 0 else 0
        row.append(f"{percentage:.2f}% ({total_count_in_group})")
    table_data.append(row)

headers = ["Group"] + [f"Subcluster {subcluster_id + 1}" for subcluster_id in range(optimal_subclusters)]

df_group = pd.DataFrame(table_data, columns=headers)

df_group


# In[49]:


print(df_group.to_latex())


# In[50]:


# Calculate total count of objects in each group within each cluster
total_counts_per_subcluster = {}
for subcluster_id in np.unique(subkmeans.labels_):
    subcluster_indices = np.where(subkmeans.labels_ == subcluster_id)[0]
    subcluster_origins = [data_type[i] for i in subcluster_indices]
    
    total_count_stars = sum(1 for obj in subcluster_origins if obj in star_groups)
    total_count_galaxies = sum(1 for obj in subcluster_origins if obj in galaxy_groups)
    
    total_counts_per_subcluster[subcluster_id] = {'Galactic': total_count_stars, 'Extragalactic': total_count_galaxies}

# List of rows for the table
table_data = []
for group in ['Galactic', 'Extragalactic']:
    row = [group]
    for subcluster_id in range(optimal_subclusters):
        total_count_in_group = total_counts_per_subcluster[subcluster_id][group]
        total_count_in_subcluster = sum(total_counts_per_subcluster[subcluster_id].values())
        percentage = (total_count_in_group / total_count_in_subcluster) * 100 if total_count_in_subcluster != 0 else 0
        row.append(f"{percentage:.2f}%")
    table_data.append(row)

headers = ["Group"] + [f"Subcluster {subcluster_id + 1}" for subcluster_id in range(optimal_subclusters)]

df_group1 = pd.DataFrame(table_data, columns=headers)

# Convert percentage strings to float values
for col in df_group1.columns[1:]:
    df_group1[col] = df_group1[col].str.rstrip('%').astype(float)
                                
# Set the index to 'Group' column for easier plotting
df_group1.set_index('Group', inplace=True)

# Plot the bar chart
ax = df_group1.T.plot(kind='bar', figsize=(10, 6))
#plt.title('Percentage of Groups in Each Cluster')
# Iterate over each column and annotate the bars with percentages
for idx, col in enumerate(df_group1.columns):
    for i in ax.patches[idx::len(df_group1.columns)]:
        ax.text(i.get_x() + i.get_width() / 2, i.get_height() + 0.5,
                f'{i.get_height()}', ha='center', va='bottom')
        
plt.xlabel('Subclusters')
plt.ylabel('Percentage / %')
plt.xticks(rotation=0)
plt.legend(title='Object Type')
plt.tight_layout()
plt.savefig('figures/subclusters_bar_7')
plt.show()


# In[51]:


# Plot the subclusters
plt.figure(figsize=(8, 6))
for cluster_num in range(optimal_subclusters):
    cluster_data = large_cluster_data[subkmeans.labels_ == cluster_num]
    plt.scatter(cluster_data[:, 0], cluster_data[:, 1], label=f'SC{cluster_num+1} ({len(cluster_data)})',s=0.5)

# Cluster centers
for i, center in enumerate(subkmeans.cluster_centers_):
    plt.scatter(center[0], center[1], marker='o', s=200, facecolors='none', edgecolors='black')
    plt.text(center[0], center[1], str(i+1), fontsize=12, color='black', ha='center', va='center')
#plt.scatter(subkmeans.cluster_centers_[:, 0], subkmeans.cluster_centers_[:, 1], marker='x', color='black', label='Centroids')
 #plt.text(center[0], center[1], str(i+1), fontsize=12, color='black', ha='center', va='center')
plt.title(f'Cluster {large_cluster_idx+1}')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.legend()
plt.savefig('figures/Subclusters_7')
plt.show()


# In[ ]:




