#!/usr/bin/env python
# coding: utf-8

# In[219]:


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


# In[ ]:


import tarfile
from astropy.io import fits

# Function to check header of the first FITS file in a TAR file
def check_header(tar_filename):
    with tarfile.open(tar_filename, 'r') as tar:
        # Get the first FITS file within the TAR file
        first_fits_member = next((member for member in tar.getmembers() if member.name.endswith('.fits')), None)
        if first_fits_member:
            # Extract the first FITS file
            tar.extract(first_fits_member)
            # Open the FITS file
            with fits.open(first_fits_member.name) as hdul:
                # Print the header
                print(hdul[1].header)
                # Close the FITS file
                hdul.close()
            # Remove the extracted FITS file
            os.remove(first_fits_member.name)

# Usage example
check_header('/its/home/ecb41/4MOST/files2/deltaCep/AT_deltaCep_S09_List1_singlespec.tar.gz')


# In[97]:


all_fluxes = []
data_type = []
root_folder = '/its/home/ecb41/4MOST/files2/Galaxies'


# In[ ]:


import csv
import re

# Open the text file in read mode
with open('/its/mnt/lustre/scratch/astro/loveday/4most/OpR2.5/opr25_metadata_05042022.txt', 'r') as file:
    # Create a CSV reader object
    reader = csv.reader(file, delimiter=' ')  # Assuming space-separated values
    
    # Initialize an empty dictionary to store the results
    results = {}
    
    # Iterate over each row in the file
    for row in reader:
        # Each row is a list of values
        
        # Check if the row contains "LJ1" and ends with ".fits"
        if "LJ1" in row[0] and row[0].endswith(".fits"):
            # Extract the redshift value from the other .fits file in the row
            other_fits = [f for f in row[1:] if f.endswith(".fits") and "LJ1" not in f][0]
            redshift_match = re.search(r'redshift([\d.]+)\.fits', other_fits)
            if redshift_match:
                redshift_value = redshift_match.group(1)
                
                # Extract the relevant values
                lj1_file = row[0]
                
                # Store the LJ1 filename and its corresponding redshift value in the dictionary
                results[lj1_file] = redshift_value

# Paginate and print the results dictionary
items_per_page = 1000  # Adjust as needed
total_items = len(results)
pages = (total_items + items_per_page - 1) // items_per_page

for page in range(pages):
    print(f"Page {page + 1}:")
    start_index = page * items_per_page
    end_index = min((page + 1) * items_per_page, total_items)
    
    for lj1_file, redshift_value in list(results.items())[start_index:end_index]:
        print(f"{lj1_file} : {redshift_value}")

    # Prompt the user to view the next page
    if page < pages - 1:
        input("Press Enter to view the next page...")


# In[ ]:


# Dictionary to store LJ1 .fits filenames and their corresponding redshift values
lj1_redshifts = {}

# Open the text file and extract LJ1 .fits filenames and their corresponding redshift values
with open('/its/mnt/lustre/scratch/astro/loveday/4most/OpR2.5/opr25_metadata_05042022.txt', 'r') as file:
    reader = csv.reader(file, delimiter=' ')
    for row in reader:
        if "LJ1" in row[0] and row[0].endswith(".fits"):
            other_fits = [f for f in row[1:] if f.endswith(".fits") and "LJ1" not in f][0]
            redshift_match = re.search(r'redshift([\d.]+)\.fits', other_fits)
            if redshift_match:
                redshift_value = float(redshift_match.group(1))  # Convert to float
                lj1_redshifts[row[0]] = redshift_value

print("lj1_redshifts:")
for lj1_file, redshift_value in list(lj1_redshifts.items())[78850:]:
    print(f"{lj1_file}: {redshift_value}")
print(len(lj1_redshifts.items()))


# In[5]:


lj1_redshifts = {}

# Open the text file and extract LJ1 .fits filenames and their corresponding redshift values
with open('/its/mnt/lustre/scratch/astro/loveday/4most/OpR2.5/opr25_metadata_05042022.txt', 'r') as file:
    reader = csv.reader(file, delimiter=' ')
    for row in reader:
        if "LJ1" in row[0] and row[0].endswith(".fits"):
            other_fits = [f for f in row[1:] if f.endswith(".fits")][0]
          #  print("Other fits:", other_fits)  # Debug statement
            # Search for 'redshift' pattern
            redshift_match = re.search(r'redshift([\d.]+)\.fits', other_fits)
            if redshift_match:
                redshift_value = float(redshift_match.group(1))
                lj1_redshifts[row[0]] = redshift_value
                continue  # Move to the next row

            # Search for 'z' followed by one or two digits, optionally followed by a decimal point and two digits pattern
           # z_match = re.search(r'z([\d]{1,2}(?:\.\d{2})?)', other_fits)
            #if z_match:
             #   redshift_value = float(z_match.group(1))  # Convert to decimal
              #  lj1_redshifts[row[0]] = redshift_value
               # continue  # Move to the next row
            # Search for 'z' followed by one or two digits, optionally followed by 'pt' and two digits pattern
            z_match = re.search(r'z([\d]{1,2}(?:pt\d{2})?)', other_fits)
            if z_match:
                # Replace 'pt' with '.'
                redshift_str = z_match.group(1).replace('pt', '.')
                redshift_value = float(redshift_str)  # Convert to decimal
                lj1_redshifts[row[0]] = redshift_value
                continue  # Move to the next row

            # Search for 'z-' followed by one or more digits, optionally followed by a decimal point and two digits pattern
            z_minus_match = re.search(r'z-([\d]+\.\d{2})', other_fits)
            if z_minus_match:
                redshift_value = float(z_minus_match.group(1))
                lj1_redshifts[row[0]] = redshift_value
                continue  # Move to the next row

print("lj1_redshifts:")
#for lj1_file, redshift_value in lj1_redshifts.items():
 #   print(f"{lj1_file}: {redshift_value}")
print(len(lj1_redshifts.items()))


# In[ ]:


threshold = 0.25
count_below_threshold = sum(1 for redshift in lj1_redshifts.values() if redshift == threshold)
print(f"Number of keys with redshift below {threshold}: {count_below_threshold}")


# In[ ]:


# Counter to track the number of filenames containing 'z-' in the same row as another FITS file with 'LJ1' in its name
count_z_minus_with_LJ1 = 0

# Open the text file and count the filenames containing 'z-' in the same row as another FITS file with 'LJ1' in its name
with open('/its/mnt/lustre/scratch/astro/loveday/4most/OpR2.5/opr25_metadata_05042022.txt', 'r') as file:
    reader = csv.reader(file, delimiter=' ')
    for row in reader:
        # Check if there's a filename in the row containing 'LJ1'
        lj1_filenames = [filename for filename in row if filename.endswith('.fits') and 'LJ1' in filename]
        if lj1_filenames:
            # Check if there's also a filename containing 'z-' in the same row
            if any(filename.endswith('.fits') and 'z-' in filename for filename in row):
                count_z_minus_with_LJ1 += 1

# Print the total count of filenames containing 'z-' in the same row as another FITS file with 'LJ1' in its name
print("Total filenames containing 'z-' in the same row as another FITS file with 'LJ1' in its name:", count_z_minus_with_LJ1)


# In[ ]:


# Counter to track the number of extracted FITS files
count = 0
all_fluxes = []
all_redshifts = []

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
        
        # Open the tar file
        with tarfile.open(tar_filepath, 'r:gz') as tar:
            # Extract all FITS files from the tar file
            fits_filenames = [member for member in tar.getmembers() if member.name.endswith('.fits') and "LJ1" in member.name]
            
            # Total number of FITS files in the tar file
            total_fits_files = len(fits_filenames)
            
            # Counter to track the number of processed FITS files in the current tar file
            processed_count = 0
            
            # Extract each FITS file
            for fits_member in fits_filenames:
                tar.extract(fits_member)
                
                # Extract flux
                flux = extract_flux(fits_member.name)
                all_fluxes.append(flux)
                
                print(f"FITS file: {fits_member.name}")  # Print the name of the extracted FITS file
    
                # Match LJ1 filename with its corresponding redshift value
                if fits_member.name in lj1_redshifts:
                    redshift_value = lj1_redshifts[fits_member.name]
                    all_redshifts.append(redshift_value)
                    print(f"Redshift: {redshift_value}")
                else:
                    all_redshifts.append(None)  # Handle cases where redshift value is not found
                
                count += 1
                processed_count += 1
                
                print(f"Extracted data from {fits_member.name}. Total extracted: {count}")
                
                # Check if all FITS files in the tar file have been processed
                if processed_count >= total_fits_files:
                    break

                # Delete the extracted FITS file -> save space
                os.remove(fits_member.name)
                
            # Check if all FITS files in the tar file have been processed
            if processed_count >= total_fits_files:
                break


# In[ ]:


# Counter to track the number of extracted FITS files
count = 0
all_fluxes = []
all_redshifts = []

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
        
        # Open the tar file
        with tarfile.open(tar_filepath, 'r:gz') as tar:
            # Extract all FITS files from the tar file
            fits_filenames = [member for member in tar.getmembers() if member.name.endswith('.fits') and "LJ1" in member.name]
            
            # Total number of FITS files in the tar file
            total_fits_files = len(fits_filenames)
            
            # Counter to track the number of processed FITS files in the current tar file
            processed_count = 0
            
            # Extract each FITS file
            for fits_member in fits_filenames:
                processed_count += 1
                # Check if the FITS file matches with LJ1 filenames in lj1_redshifts
                if fits_member.name in lj1_redshifts:
                    tar.extract(fits_member)
                    
                    # Extract flux
                    flux = extract_flux(fits_member.name)
                    all_fluxes.append(flux)
                    
                    print(f"FITS file: {fits_member.name}")  # Print the name of the extracted FITS file
        
                    # Match LJ1 filename with its corresponding redshift value
                    redshift_value = lj1_redshifts[fits_member.name]
                    all_redshifts.append(redshift_value)
                    print(f"Redshift: {redshift_value}")
                    
                    count += 1
                    
                    
                    print(f"Extracted data from {fits_member.name}. Total extracted: {count}")
                    
                    # Delete the extracted FITS file -> save space
                    os.remove(fits_member.name)
                    
                # Check if all FITS files in the tar file have been processed
                if processed_count >= total_fits_files:
                    break

                
            


# In[ ]:


import numpy as np  # Import NumPy for numerical operations

# Counter to track the number of extracted FITS files
count = 0
all_fluxes = []
all_redshifts = []
all_wavelengths = []
rest_frame_fluxes = []
rest_frame_wavelengths = []

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
        
        # Open the tar file
        with tarfile.open(tar_filepath, 'r:gz') as tar:
            # Extract all FITS files from the tar file
            fits_filenames = [member for member in tar.getmembers() if member.name.endswith('.fits') and "LJ1" in member.name]
            
            # Total number of FITS files in the tar file
            total_fits_files = len(fits_filenames)
            
            # Counter to track the number of processed FITS files in the current tar file
            processed_count = 0
            
            # Extract each FITS file
            for fits_member in fits_filenames:
                processed_count += 1
                # Check if the FITS file matches with LJ1 filenames in lj1_redshifts
                if fits_member.name in lj1_redshifts:
                    tar.extract(fits_member)
                    
                    # Extract flux and wavelength data
                    flux, wavelengths = extract_flux_wave(fits_member.name)
                    all_fluxes.append(flux)
                    all_wavelengths.append(wavelengths)
                    
                    #print(f"FITS file: {fits_member.name}")  # Print the name of the extracted FITS file
        
                    # Match LJ1 filename with its corresponding redshift value
                    redshift_value = lj1_redshifts[fits_member.name]
                    all_redshifts.append(redshift_value)
                    #print(f"Redshift: {redshift_value}")
                    
                    # Correct flux and wavelength for redshift
                    #rest_frame_flux = flux / (1 + redshift_value)
                    rest_frame_wavelength = wavelengths / (1 + redshift_value)
                    
                    rest_frame_fluxes.append(rest_frame_flux)
                    rest_frame_wavelengths.append(rest_frame_wavelength)
                    
                    data_type.append(folder_name)
                    
                   # print(f"Rest-frame Flux: {rest_frame_flux}")
                   # print(f"Rest-frame Wavelength: {rest_frame_wavelength}")
                    
                    count += 1
                    
                    print(f"Extracted data from {fits_member.name}. Total extracted: {count}")
                    
                    # Delete the extracted FITS file -> save space
                    os.remove(fits_member.name)
                    
                # Check if all FITS files in the tar file have been processed
                if processed_count >= total_fits_files:
                    break


# In[ ]:


rest_frame_fluxes = np.asarray(rest_frame_fluxes)
print(rest_frame_fluxes.shape)
rest_frame_wavelengths = np.asarray(rest_frame_wavelengths)
print(rest_frame_wavelengths.shape)


# In[ ]:


min = ((wavelengths[0] for wavelengths in rest_frame_wavelengths))
print(min)


# In[ ]:


# Step 1: Determine the minimum and maximum wavelength values
min_wavelength = 3000
max_wavelength = 9000

for wavelengths in rest_frame_wavelengths:
    if np.min(wavelengths) < min_wavelength:
        min_wavelength = np.min(wavelengths)
    if np.max(wavelengths) > max_wavelength:
        max_wavelength = np.max(wavelengths)

# Step 2: Define the common wavelength range
desired_num_points = 2000
common_wavelength_range = np.linspace(min_wavelength, max_wavelength, num=desired_num_points)

# Step 3: Interpolate each flux spectrum onto the common wavelength range
truncated_fluxes = []
for wavelengths, fluxes in zip(rest_frame_wavelengths, rest_frame_fluxes):
    interpolated_flux = np.interp(common_wavelength_range, wavelengths, fluxes)
    truncated_fluxes.append(interpolated_flux)


# In[ ]:


# Find the overall minimum and maximum wavelengths
min_wavelength = (min(wavelengths) for wavelengths in rest_frame_wavelengths)
max_wavelength = (max(wavelengths) for wavelengths in rest_frame_wavelengths)
print(max(all_redshifts))
print(min_wavelength,max_wavelength)

# Truncate the fluxes to the common wavelength range
truncated_fluxes = []
for wavelengths, fluxes in zip(rest_frame_wavelengths, rest_frame_fluxes):
    truncated_flux = []
    for wavelength, flux in zip(wavelengths, fluxes):
        if min_wavelength <= wavelength <= max_wavelength:
            truncated_flux.append(flux)
    truncated_fluxes.append(truncated_flux)


# In[ ]:


# Flatten the list of wavelengths for all spectra
all_wavelengths_flat = np.concatenate(rest_frame_wavelength)

# Find the overall minimum and maximum wavelengths
min_wavelength = np.min(all_wavelengths_flat)
max_wavelength = np.max(all_wavelengths_flat)

# Define the common wavelength range
common_wavelength_range = (min_wavelength, max_wavelength)
print(common_wavelength_range)


# In[98]:


from specutils import Spectrum1D
from astropy import units as u

count = 0
all_spectra = []

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
                    
                    if redshift_value <= 0.25 :
                    
                        tar.extract(fits_member)
                    
                        # Extract flux and wavelength data
                        flux, wavelengths = extract_flux_wave(fits_member.name)
                        
                        # Convert flux to float32
                        flux = flux.astype(np.float32)
                        
                        # Correct flux and wavelength for redshift
                        #rest_frame_flux = flux / (1 + redshift_value)
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


# In[99]:


root_folder2 = '/its/home/ecb41/4MOST/files2/Stars'


# In[100]:


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


# In[ ]:


flux_data = np.array([spectrum.flux.value for spectrum in all_spectra])


# In[ ]:


# Find the overall minimum and maximum wavelengths
min_wavelength = (min(spectrum.spectral_axis) for spectrum in all_spectra)
max_wavelength = (max(spectrum.spectral_axis) for spectrum in all_spectra)

mi = max(min_wavelength)
ma = min(max_wavelength)

print(mi,ma)


# In[10]:


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


# In[ ]:


#from specutils.manipulation import extract_region

# Find the overall minimum and maximum wavelengths
#min_wavelength = min(min(spectrum.spectral_axis) for spectrum in all_spectra)
#max_wavelength = max(max(spectrum.spectral_axis) for spectrum in all_spectra)

#min_wavelength = 3800 
#max_wavelength = 8000 

# Define the common wavelength range
#common_wavelength_range = (min_wavelength, max_wavelength) * u.AA

# Truncate the spectra to the common wavelength range
#truncated_spectra = []
#for spectrum in all_spectra:
 #   truncated_spectrum = extract_region(spectrum, common_wavelength_range)
  #  truncated_spectra.append(truncated_spectrum)


# In[ ]:


#truncated_fluxes = np.asarray(truncated_fluxes)
#print(truncated_fluxes.shape)


# In[11]:


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


# In[21]:


print(flux_data[29])


# In[ ]:


# Function to plot a spectrum
def plot_spectrum(spectrum, ax=None, **kwargs):
    if ax is None:
        ax = plt.gca()
    ax.plot(spectrum.spectral_axis, spectrum.flux, **kwargs)
    ax.set_xlabel('Wavelength (Angstrom)')
    ax.set_ylabel('Flux')
    ax.set_title('Resampled Spectrum')
    ax.grid(True)

# Plot some resampled spectra
num_spectra_to_plot = 5
fig, axs = plt.subplots(num_spectra_to_plot, figsize=(10, 8), sharex=True, sharey=True)

for i in range(num_spectra_to_plot):
    plot_spectrum(resampled_spectra[i], ax=axs[i])

plt.tight_layout()
plt.show()


# In[12]:


# Normalize flux data
scaler = StandardScaler()
normalized_flux = scaler.fit_transform(flux_data)


# ## PCA

# In[146]:


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


# In[104]:


# Variance in each PC
variance_df = pd.DataFrame({'Principal Component': range(1,11),
                            'Explained Variance Ratio': variance[:10]})

variance_df


# In[105]:


print(variance_df.to_latex())


# In[216]:


# Compute the cumulative explained variance ratio
Cumalative_variance = pca.explained_variance_ratio_.cumsum()

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


# ## Elbow Method

# In[106]:


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

# In[145]:


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

# In[231]:


range_n_clusters = [6,7,8,9]
X = principal_components[:, :2]

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

plt.savefig("figures/silhouette")
plt.show()


# In[226]:


range_n_clusters = [6, 7, 8, 9]
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


# In[227]:


print(df_silhouette.to_latex())


# ## Clusters

# In[151]:


# Reorder data_type
order = ['whitedwarf','F', 'G', 'K', 'SF', 'Pass', 'SN']


# In[152]:


# Perform k-means clustering on the PCA-transformed data
n_clusters = 7  # Number of clusters
kmeans = KMeans(n_clusters=n_clusters)
kmeans.fit(principal_components[:, :2])

# Count the number of data points in each cluster
cluster_counts = np.bincount(kmeans.labels_)

# Convert the counts to a dictionary
cluster_counts_dict = {f"Cluster {i + 1}": count for i, count in enumerate(cluster_counts)}

# Convert the dictionary into a DataFrame
df = pd.DataFrame.from_dict(cluster_counts_dict, orient='index', columns=['Number of Spectra'])

df


# In[153]:


print(df.to_latex())


# In[154]:


# Count the number of each type of object in each cluster
#for cluster_id in np.unique(kmeans.labels_):
  #  cluster_indices = np.where(kmeans.labels_ == cluster_id)[0]
 #   cluster_origins = [data_type[i] for i in cluster_indices]
#    cluster_counts_dict[f"Cluster {cluster_id + 1}"] = {obj: cluster_origins.count(obj) for obj in set(cluster_origins)}


# Count the number of each type of object in each cluster
for cluster_id in np.unique(kmeans.labels_):
    cluster_indices = np.where(kmeans.labels_ == cluster_id)[0]
    cluster_origins = [data_type[i] for i in cluster_indices]
    counts = {obj: cluster_origins.count(obj) for obj in order}
    cluster_counts_dict[f"Cluster {cluster_id + 1}"] = counts

# Convert the dictionary into a DataFrame
df = pd.DataFrame.from_dict(cluster_counts_dict, orient='index')

df


# In[155]:


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


# In[156]:


print(df_percentage_counts.to_latex())


# In[157]:


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


# In[158]:


print(df_group_counts.to_latex())


# In[217]:


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


# In[203]:


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

# In[161]:


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
plt.savefig('figures/Dunn_SC')
plt.show()


# In[232]:


range_n_subclusters = [5,6,7,8]
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
        "k = %d"
        % n_clusters,
        fontsize=14,
        fontweight="bold",
    )

plt.savefig("figures/silhouette_sc")
plt.show()


# In[229]:


range_n_subclusters = [5,6,7,8]
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


# In[230]:


print(df_silhouette2.to_latex())


# In[187]:


# Initial K-means clustering
#n_clusters = 7  
#kmeans = KMeans(n_clusters=n_clusters)
#kmeans.fit(principal_components[:, :2])

# Identify the large cluster
large_cluster_idx = np.argmax(np.bincount(kmeans.labels_))

# Extract data points belonging to the large cluster
large_cluster_mask = kmeans.labels_ == large_cluster_idx
large_cluster_data = principal_components[large_cluster_mask]

optimal_subclusters = 8 

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


# In[188]:


print(df.to_latex())


# In[191]:


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


# In[192]:


print(df_percentage.to_latex())


# In[194]:


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


# In[195]:


print(df_group.to_latex())


# In[218]:


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
plt.savefig('figures/subclusters_bar')
plt.show()


# In[204]:


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
plt.savefig('figures/Subclusters')
plt.show()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[25]:


subclusters = subkmeans.fit_predict(large_cluster_data)
# Count the number of each type of object in each subcluster
subcluster_counts = {}
for subcluster_id in range(optimal_subclusters):
    subcluster_indices = np.where(subclusters == subcluster_id)[0]
    subcluster_origins = [data_type[i] for i in np.where(large_cluster_mask)[0][subcluster_indices]]
    subcluster_counts[subcluster_id] = {obj: subcluster_origins.count(obj) for obj in set(subcluster_origins)}

for subcluster_id, counts in subcluster_counts.items():
    print(f"Subcluster {subcluster_id}: {counts}")


# In[26]:


# Calculate total count of each object across all clusters
total_counts = {}
for subcluster_counts_dict in subcluster_counts.values():
    for obj, count in subcluster_counts_dict.items():
        total_counts[obj] = total_counts.get(obj, 0) + count

# list of rows for the table
table_data = []
for obj in total_counts.keys():
    row = [obj]
    for subcluster_id in range(optimal_subclusters):
        count = subcluster_counts[subcluster_id].get(obj, 0)
        percentage = (count / total_counts[obj]) * 100 if total_counts[obj] != 0 else 0
        row.append(f"{percentage:.2f}%")
    table_data.append(row)


headers = ["Object"] + [f"Subcluster {subcluster_id}" for subcluster_id in range(optimal_subclusters)]

df = pd.DataFrame(table_data, columns=headers)

df


# In[27]:


# Calculate total count of each object within each cluster
total_counts_per_cluster = {}
for subcluster_id, subcluster_counts_dict in subcluster_counts.items():
    total_count = sum(subcluster_counts_dict.values())
    total_counts_per_cluster[subcluster_id] = total_count

# List of rows for the table
table_data = []
for obj in total_counts.keys():
    row = [obj]
    for subcluster_id in range(optimal_subclusters):
        count = subcluster_counts[subcluster_id].get(obj, 0)
        total_count_in_cluster = total_counts_per_cluster.get(subcluster_id, 0)
        percentage = (count / total_count_in_cluster) * 100 if total_count_in_cluster != 0 else 0
        row.append(f"{percentage:.2f}%")
    table_data.append(row)

headers = ["Object"] + [f"Subcluster {subcluster_id}" for subcluster_id in range(optimal_subclusters)]

df = pd.DataFrame(table_data, columns=headers)

df


# In[34]:


# Define the groups
star_groups = {'F', 'G', 'K', 'whitedwarf'}
galaxy_groups = {'Pass', 'SN', 'SF'}

# Calculate total count of objects in each group within each cluster
total_counts_per_cluster = {}
for subcluster_id, subcluster_counts_dict in subcluster_counts.items():
    total_count_stars = sum(count for obj, count in subcluster_counts_dict.items() if obj in star_groups)
    total_count_galaxies = sum(count for obj, count in subcluster_counts_dict.items() if obj in galaxy_groups)
    total_counts_per_cluster[subcluster_id] = {'Galactic': total_count_stars, 'Extragalactic': total_count_galaxies}

# List of rows for the table
table_data = []
for group in ['Galactic', 'Extragalactic']:
    row = [group]
    for subcluster_id in range(optimal_subclusters):
        total_count_in_group = total_counts_per_cluster[subcluster_id][group]
        total_count_in_cluster = sum(total_counts_per_cluster[subcluster_id].values())
        percentage = (total_count_in_group / total_count_in_cluster) * 100 if total_count_in_cluster != 0 else 0
        row.append(f"{percentage:.2f}%")
    table_data.append(row)

headers = ["Group"] + [f"Subcluster {subcluster_id}" for subcluster_id in range(optimal_subclusters)]

df = pd.DataFrame(table_data, columns=headers)

df


# In[ ]:


print(df.to_latex())


# In[38]:


from sklearn.metrics import confusion_matrix

# Define the true labels and predicted labels based on the percentage of objects in each cluster
true_labels = []  # True labels based on the actual object types
predicted_labels = []  # Predicted labels based on the cluster assignments

# Assign true labels based on the object types
for obj in total_counts.keys():
    if obj in star_groups:
        true_labels.append('Galactic')
    elif obj in galaxy_groups:
        true_labels.append('Extragalactic')
print(total_counts.keys())
print(true_labels)
# Assign predicted labels based on the cluster assignments
for subcluster_id in range(optimal_subclusters):
    total_count_stars = sum(count for obj, count in subcluster_counts[subcluster_id].items() if obj in star_groups)
    total_count_galaxies = sum(count for obj, count in subcluster_counts[subcluster_id].items() if obj in galaxy_groups)
    if total_count_stars > total_count_galaxies:
        predicted_labels.append('Galactic')
    else:
        predicted_labels.append('Extragalactic')
print(predicted_labels)
# Create the confusion matrix
conf_matrix = confusion_matrix(true_labels, predicted_labels, labels=['Galactic', 'Extragalactic'])

print("Confusion Matrix:")
print(conf_matrix)


# In[36]:


from sklearn.metrics import ConfusionMatrixDisplay

# Display the confusion matrix
disp = ConfusionMatrixDisplay(confusion_matrix=conf_matrix, display_labels=['Star', 'Galaxy'])
disp.plot()
plt.title('Confusion Matrix')
plt.show()


# In[43]:


import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Define the original groups
original_groups = ['G', 'SF', 'K', 'Pass', 'F', 'whitedwarf', 'SN']

provided_counts = {0: {'G': 1419, 'SF': 1027, 'K': 801, 'Pass': 2754, 'F': 955, 'whitedwarf': 3, 'SN': 14},
1: {'G': 20, 'SF': 86, 'F': 66, 'whitedwarf': 2862, 'SN': 81},
2: {'G': 444, 'SF': 108, 'K': 47, 'Pass': 397, 'F': 305, 'whitedwarf': 48},
3: {'G': 384, 'SF': 2460, 'K': 26, 'Pass': 391, 'F': 230, 'whitedwarf': 1918, 'SN': 4137},
4: {'F': 20, 'SF': 2, 'whitedwarf': 1132, 'SN': 1},
5: {'G': 420, 'K': 377, 'F': 247},
6: {'G': 723, 'SF': 7070, 'K': 451, 'Pass': 3430, 'F': 521, 'whitedwarf': 249, 'SN': 4288}}    


# In[45]:


# Initialize a dictionary to store counts for each subcluster
subcluster_counts = {subcluster_id: {group: 0 for group in original_groups} for subcluster_id in range(optimal_subclusters)}

# Populate the counts dictionary
for subcluster_id, counts_dict in subcluster_counts.items():
    provided_counts_subcluster = provided_counts.get(subcluster_id, {})
    for group, count in counts_dict.items():
        subcluster_counts[subcluster_id][group] = provided_counts_subcluster.get(group, 0)

# Convert counts to percentages
total_counts = {group: sum(counts[group] for counts in subcluster_counts.values()) for group in original_groups}
percentage_counts = {subcluster_id: {group: (counts[group] / total_counts[group]) * 100 for group in original_groups} for subcluster_id, counts in subcluster_counts.items()}

# Create a confusion matrix plot
fig, ax = plt.subplots(figsize=(10, 8))
sns.heatmap(pd.DataFrame(percentage_counts).transpose(), annot=True, cmap="YlGnBu", fmt=".2f", ax=ax)
plt.title('Percentage of Original Groups in Each Subcluster')
plt.xlabel('Original Groups')
plt.ylabel('Subclusters')
plt.xticks(rotation=45)
plt.yticks(rotation=0)
plt.tight_layout()
plt.show()


# In[46]:


# Define the groups
star_groups = {'F', 'G', 'K', 'whitedwarf'}
galaxy_groups = {'Pass', 'SN', 'SF'}

# Initialize a dictionary to store counts for each subcluster and group
provided_counts = {subcluster_id: {'Galactic': 0, 'Extragalactic': 0} for subcluster_id in range(optimal_subclusters)}

# Populate the counts dictionary with the provided counts
for subcluster_id, subcluster_counts_dict in subcluster_counts.items():
    provided_counts[subcluster_id]['Galactic'] = sum(count for obj, count in subcluster_counts_dict.items() if obj in star_groups)
    provided_counts[subcluster_id]['Extragalactic'] = sum(count for obj, count in subcluster_counts_dict.items() if obj in galaxy_groups)

# Define the original groups
original_groups = ['Galactic', 'Extragalactic']

# Initialize a dictionary to store counts for each subcluster
subcluster_counts = {subcluster_id: {group: 0 for group in original_groups} for subcluster_id in range(optimal_subclusters)}

# Populate the counts dictionary
for subcluster_id, counts_dict in subcluster_counts.items():
    provided_counts_subcluster = provided_counts.get(subcluster_id, {})
    for group, count in counts_dict.items():
        subcluster_counts[subcluster_id][group] = provided_counts_subcluster.get(group, 0)

# Convert counts to percentages
total_counts = {group: sum(counts[group] for counts in subcluster_counts.values()) for group in original_groups}
percentage_counts = {subcluster_id: {group: (counts[group] / total_counts[group]) * 100 for group in original_groups} for subcluster_id, counts in subcluster_counts.items()}

# Create a confusion matrix plot
fig, ax = plt.subplots(figsize=(8, 6))
sns.heatmap(pd.DataFrame(percentage_counts).transpose(), annot=True, cmap="YlGnBu", fmt=".2f", ax=ax)
plt.title('Percentage of Stars and Galaxies in Each Subcluster')
plt.xlabel('Original Groups')
plt.ylabel('Subclusters')
plt.xticks(rotation=45)
plt.yticks(rotation=0)
plt.tight_layout()
plt.show()


# In[47]:


# Convert percentages in the DataFrame to numerical values
for col in df.columns[1:]:
    df[col] = df[col].str.rstrip('%').astype(float)

# Set the index to 'Group'
df.set_index('Group', inplace=True)

# Create a heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(df, annot=True, cmap='YlGnBu', fmt=".2f")
plt.title('Percentage of Objects in Each Subcluster by Group')
plt.xlabel('Subclusters')
plt.ylabel('Group')
plt.xticks(rotation=45)
plt.yticks(rotation=0)
plt.tight_layout()
plt.show()


# In[53]:


df = pd.DataFrame({
    'Group': ['Galactic', 'Extragalactic'],
    'Subcluster 0': ['45.58%', '54.42%'],
    'Subcluster 1': ['94.64%', '5.36%'],
    'Subcluster 2': ['62.56%', '37.44%'],
    'Subcluster 3': ['26.80%', '73.20%'],
    'Subcluster 4': ['99.74%', '0.26%'],
    'Subcluster 5': ['100%', '0.00%'],
    'Subcluster 6': ['11.62%', '88.38%']
})
# Convert percentage strings to float values
for col in df.columns[1:]:
    df[col] = df[col].str.rstrip('%').astype(float)

# Set the index to 'Group' column for easier plotting
df.set_index('Group', inplace=True)

# Plot the bar chart
ax = df.T.plot(kind='bar', figsize=(10, 6))
plt.title('Percentage of Groups in Each Subcluster')
plt.xlabel('Subclusters')
plt.ylabel('Percentage')
plt.xticks(rotation=0)
plt.legend(title='Group')
plt.tight_layout()
plt.savefig('figures/Gal_ex_bar')
plt.show()


# In[56]:


# Calculate total count of each object within each cluster
total_counts_per_cluster = {}
for subcluster_id, subcluster_counts_dict in subcluster_counts.items():
    total_count = sum(subcluster_counts_dict.values())
    total_counts_per_cluster[subcluster_id] = total_count

# List of rows for the table
table_data = []
for obj in total_counts.keys():
    row = [obj]
    for subcluster_id in range(optimal_subclusters):
        count = subcluster_counts[subcluster_id].get(obj, 0)
        total_count_in_cluster = total_counts_per_cluster.get(subcluster_id, 0)
        percentage = (count / total_count_in_cluster) * 100 if total_count_in_cluster != 0 else 0
        row.append(f"{percentage:.2f}%")
    table_data.append(row)

headers = ["Object"] + [f"Subcluster {subcluster_id}" for subcluster_id in range(optimal_subclusters)]

df = pd.DataFrame(table_data, columns=headers)  
# Convert percentage strings to float values
for col in df.columns[1:]:
    df[col] = df[col].str.rstrip('%').astype(float)

# Set the index to 'Group' column for easier plotting
df.set_index('Object', inplace=True)

# Plot the bar chart
ax = df.T.plot(kind='bar', figsize=(10, 6))
plt.title('Percentage of Groups in Each Subcluster')
plt.xlabel('Subclusters')
plt.ylabel('Percentage')
plt.xticks(rotation=0)
plt.legend(title='Group')
plt.tight_layout()
plt.show()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


from specutils.manipulation import FluxConservingResampler

disp_grid = np.arange(3800, 8000, 3) * u.AA
fluxcon = FluxConservingResampler()
spec = fluxcon(all_spectra, disp_grid)


# In[ ]:


print(len(all_fluxes))
print(len(all_redshifts))
print(all_redshifts[ :1000])


# In[ ]:


# Dictionary to store LJ1 .fits filenames and their corresponding redshift values
lj1_redshifts = {}

# Open the text file and extract LJ1 .fits filenames and their corresponding redshift values
with open('/its/mnt/lustre/scratch/astro/loveday/4most/OpR2.5/opr25_metadata_05042022.txt', 'r') as file:
    reader = csv.reader(file, delimiter=' ')
    for row in reader:
        if "LJ1" in row[0] and row[0].endswith(".fits"):
            other_fits = [f for f in row[1:] if f.endswith(".fits") and "LJ1" not in f][0]
            redshift_match = re.search(r'redshift([\d.]+)\.fits', other_fits)
            if redshift_match:
                redshift_value = float(redshift_match.group(1))  # Convert to float
                lj1_redshifts[row[0]] = redshift_value

print("lj1_redshifts:")
for lj1_file, redshift_value in list(lj1_redshifts.items())[:10]:
    print(f"{lj1_file}: {redshift_value}")

# Counter to track the number of extracted FITS files
count = 0
all_fluxes = []
all_redshifts = []

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
        
        # Open the tar file
        with tarfile.open(tar_filepath, 'r:gz') as tar:
            # Extract all FITS files from the tar file
            fits_filenames = [member for member in tar.getmembers() if member.name.endswith('.fits') and "LJ1" in member.name]
            
            # Total number of FITS files in the tar file
            total_fits_files = len(fits_filenames)
            
            # Counter to track the number of processed FITS files in the current tar file
            processed_count = 0
            
            # Extract each FITS file
            for fits_member in fits_filenames:
                tar.extract(fits_member)
                
                # Extract flux
                flux = extract_flux(fits_member.name)
                all_fluxes.append(flux)
                
                print(f"FITS file: {fits_member.name}")  # Print the name of the extracted FITS file
    
                # Match LJ1 filename with its corresponding redshift value
                if fits_member.name in lj1_redshifts:
                    redshift_value = lj1_redshifts[fits_member.name]
                    all_redshifts.append(redshift_value)
                    print(f"Redshift: {redshift_value}")
                else:
                    print(f"No redshift found for {fits_member.name}")  # Debugging message
                    all_redshifts.append(None)  # Handle cases where redshift value is not found
                
                count += 1
                processed_count += 1
                
                print(f"Extracted data from {fits_member.name}. Total extracted: {count}")
                
                # Check if all FITS files in the tar file have been processed
                if processed_count >= total_fits_files:
                    break

                # Delete the extracted FITS file -> save space
                os.remove(fits_member.name)
                
            # Check if all FITS files in the tar file have been processed
            if processed_count >= total_fits_files:
                break


# In[ ]:


# Counter to track the number of matches
match_count = 0

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
        
        # Open the tar file
        with tarfile.open(tar_filepath, 'r:gz') as tar:
            # Extract all FITS files from the tar file
            fits_filenames = [member.name for member in tar.getmembers() if member.name.endswith('.fits') and "LJ1" in member.name]
            
            # Check for matches
            for fits_member_name in fits_filenames:
                if fits_member_name in lj1_redshifts:
                    match_count += 1
        print(match_count)
# Print the total number of matches
print("Total number of matches:", match_count)


# In[ ]:


# Dictionary to store LJ1 .fits filenames and their corresponding redshift values
lj1_redshifts = {}

# Open the text file from the previous code snippet and extract the LJ1 .fits filenames and redshift values
with open('/its/mnt/lustre/scratch/astro/loveday/4most/OpR2.5/opr25_metadata_05042022.txt', 'r') as file:
    reader = csv.reader(file, delimiter=' ')
    for row in reader:
        if "LJ1" in row[0] and row[0].endswith(".fits"):
            other_fits = [f for f in row[1:] if f.endswith(".fits") and "LJ1" not in f][0]
            redshift_match = re.search(r'redshift([\d.]+)\.fits', other_fits)
            if redshift_match:
                redshift_value = redshift_match.group(1)
                lj1_redshifts[row[0]] = redshift_value
print(len(lj1_redshifts))


# In[ ]:


# Dictionary to store LJ1 .fits filenames and their corresponding redshift values
lj1_redshifts = {}

# Open the text file from the previous code snippet and extract the LJ1 .fits filenames and redshift values
with open('/its/mnt/lustre/scratch/astro/loveday/4most/OpR2.5/opr25_metadata_05042022.txt', 'r') as file:
    reader = csv.reader(file, delimiter=' ')
    for row in reader:
        if "LJ1" in row[0] and row[0].endswith(".fits"):
            other_fits = [f for f in row[1:] if f.endswith(".fits") and "LJ1" not in f][0]
            redshift_match = re.search(r'redshift([\d.]+)\.fits', other_fits)
            if redshift_match:
                redshift_value = redshift_match.group(1)
                lj1_redshifts[row[0]] = redshift_value
print(lj1_redshifts[:10])
# Counter to track the number of extracted FITS files
count = 0
all_fluxes = []
all_redshifts = []

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
        
        # Open the tar file
        with tarfile.open(tar_filepath, 'r:gz') as tar:
            # Extract all FITS files from the tar file
            fits_filenames = [member for member in tar.getmembers() if member.name.endswith('.fits') and "LJ1" in member.name]
            
            # Extract each FITS file
            for fits_member in fits_filenames:
                tar.extract(fits_member)
                
                # Extract flux
                flux = extract_flux(fits_member.name)
                all_fluxes.append(flux)
                
                # Match LJ1 filename with its corresponding redshift value
                if fits_member.name in lj1_redshifts:
                    redshift_value = lj1_redshifts[fits_member.name]
                    all_redshifts.append(redshift_value)
                else:
                    all_redshifts.append(None)  # Handle cases where redshift value is not found
                
                count += 1
                
                print(f"Extracted data from {fits_member.name}. Total extracted: {count}")
                
                # Check if the count exceeds 50,000 and break the loop if so
                if count >= 50000:
                    break

                # Delete the extracted FITS file -> save space
                os.remove(fits_member.name)
                
            # Check if the count exceeds 50,000 and break the outer loop if so
            if count >= 50000:
                break
                
    # Check if the count exceeds 50,000 and break the outermost loop if so
    if count >= 50000:
        break


# In[ ]:


print(len(all_fluxes))
print(len(all_redshifts))
print(all_redshifts[ :1000])


# In[ ]:


import numpy as np

# Convert all_fluxes and all_redshifts to NumPy arrays for efficient element-wise operations
all_fluxes = np.array(all_fluxes)
all_redshifts = np.array(all_redshifts)

# De-redshift the fluxes using NumPy
all_de_redshifted_fluxes = all_fluxes / (1 + all_redshifts[:, np.newaxis])


# In[ ]:


# Dictionary to store LJ1 .fits filenames and their corresponding redshift values
lj1_redshifts = {}

# Open the text file and extract LJ1 .fits filenames and their corresponding redshift values
with open('/its/mnt/lustre/scratch/astro/loveday/4most/OpR2.5/opr25_metadata_05042022.txt', 'r') as file:
    reader = csv.reader(file, delimiter=' ')
    for row in reader:
        if "LJ1" in row[0] and row[0].endswith(".fits"):
            other_fits = [f for f in row[1:] if f.endswith(".fits") and "LJ1" not in f][0]
            redshift_match = re.search(r'redshift([\d.]+)\.fits', other_fits)
            if redshift_match:
                redshift_value = float(redshift_match.group(1))  # Convert to float
                lj1_redshifts[row[0]] = redshift_value

print("lj1_redshifts:")
#print(lj1_redshifts)
for lj1_file, redshift_value in list(lj1_redshifts.items())[:10]:
    print(f"{lj1_file}: {redshift_value}")
# Counter to track the number of extracted FITS files
count = 0
all_fluxes = []
all_redshifts = []

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
        
        # Open the tar file
        with tarfile.open(tar_filepath, 'r:gz') as tar:
            # Extract all FITS files from the tar file
            fits_filenames = [member for member in tar.getmembers() if member.name.endswith('.fits') and "LJ1" in member.name]
            
            # Extract each FITS file
            for fits_member in fits_filenames:
                tar.extract(fits_member)
                
                # Extract flux
                flux = extract_flux(fits_member.name)
                all_fluxes.append(flux)
                
                print(f"FITS file: {fits_member.name}")  # Print the name of the extracted FITS file
    
                # Match LJ1 filename with its corresponding redshift value
                if fits_member.name in lj1_redshifts:
                    redshift_value = lj1_redshifts[fits_member.name]
                    all_redshifts.append(redshift_value)
                    print(f"Redshift: {redshift_value}")
                else:
                    all_redshifts.append(None)  # Handle cases where redshift value is not found
                
                count += 1
                
                print(f"Extracted data from {fits_member.name}. Total extracted: {count}")
                
                # Check if the count exceeds 50,000 and break the loop if so
                if count >= 50000:
                    break

                # Delete the extracted FITS file -> save space
                os.remove(fits_member.name)
                
            # Check if the count exceeds 50,000 and break the outer loop if so
            if count >= 50000:
                break
                
    # Check if the count exceeds 50,000 and break the outermost loop if so
    if count >= 50000:
        break


# In[ ]:


print(len(all_fluxes))
print(len(all_redshifts))
print(all_redshifts[ :1000])


# In[ ]:


# Dictionary to store LJ1 .fits filenames and their corresponding redshift values
lj1_redshifts = {}

# Open the text file and extract LJ1 .fits filenames and their corresponding redshift values
with open('/its/mnt/lustre/scratch/astro/loveday/4most/OpR2.5/opr25_metadata_05042022.txt', 'r') as file:
    reader = csv.reader(file, delimiter=' ')
    for row in reader:
        if "LJ1" in row[0] and row[0].endswith(".fits"):
            other_fits = [f for f in row[1:] if f.endswith(".fits") and "LJ1" not in f][0]
            print(other_fits)  # Debug print statement
            redshift_match = re.search(r'redshift([\d.]+)\.fits', other_fits)
            if redshift_match:
                redshift_value = float(redshift_match.group(1))
                lj1_redshifts[row[0]] = redshift_value


# In[ ]:


import glob
# Define the filename to check
filename_to_check = "20220904/singlespec/qmost_20581800-4119150_20220904_100114_LJ1.fits"

# Check if the file exists in the text file
file_found_in_text_file = False
with open('/its/mnt/lustre/scratch/astro/loveday/4most/OpR2.5/opr25_metadata_05042022.txt', 'r') as file:
    reader = csv.reader(file, delimiter=' ')
    for row in reader:
        if filename_to_check in row:
            file_found_in_text_file = True
            break

# Check if the file exists in any of the tar files
file_found_in_tar_files = False
root_folder = '/its/home/ecb41/4MOST/files2'  # Update with the path to the root folder
tar_files = glob.glob(os.path.join(root_folder, "**/*.tar.gz"), recursive=True)
for tar_filepath in tar_files:
    # Open the tar file
    with tarfile.open(tar_filepath, 'r:gz') as tar:
        # Check if the file exists in the tar file
        if filename_to_check in tar.getnames():
            file_found_in_tar_files = True
            break

# Print the results
print(f"File {filename_to_check} found in text file: {file_found_in_text_file}")
print(f"File {filename_to_check} found in tar files: {file_found_in_tar_files}")


# In[ ]:





# In[ ]:


# Open the text file and extract LJ1 .fits filenames and their corresponding redshift values
with open('/its/mnt/lustre/scratch/astro/loveday/4most/OpR2.5/opr25_metadata_05042022.txt', 'r') as file:
    reader = csv.reader(file, delimiter=' ')
    for row in reader:
        if "LJ1" in row[0] and row[0].endswith(".fits"):
            other_fits = [f for f in row[1:] if f.endswith(".fits") and "LJ1" not in f][0]
            print(f"Other FITS filename: {other_fits}")  # Print the other FITS filename for debugging
            redshift_match = re.search(r'redshift(\d+\.\d+)\.fits', other_fits)
            if redshift_match:
                redshift_value = redshift_match.group(1)
                lj1_redshifts[row[0]] = redshift_value
               # print(f"Redshift value: {redshift_value}")  # Print the extracted redshift value for debugging
        count += 1
        if count >= 10:  # Stop after printing the first 10 other_fits filenames
            break
    
print(lj1_redshifts[:10])  # Print a portion of lj1_redshifts for debugging


# In[ ]:




