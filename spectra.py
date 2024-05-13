#!/usr/bin/env python
# coding: utf-8

# In[4]:


import os
import glob
import tarfile
from astropy.io import fits
import matplotlib.pyplot as plt

name = ['F', 'G', 'K', 'Pass', 'SF', 'SN', 'Pass','whitedwarf']

# Function to plot flux vs. wavelength
def plot_flux_vs_wavelength(ax, wavelength, flux, label):
    ax.plot(wavelength, flux, label=label)
    ax.set_xlabel('Wavelength / Å')
    ax.set_ylabel('Flux / $\mathrm{erg \, s^{-1} \, cm^{-2} \, \AA^{-1}}$')  
    ax.legend(loc='upper right')
    
  #  ax.set_title('Flux vs. Wavelength')

# Directory containing the subdirectories
subdirectories_dir = '/its/home/ecb41/4MOST/files2/Galaxies/*'

# Number of subplots per row and column
num_rows = 7
num_cols = 1

# Create subplots
fig, axs = plt.subplots(num_rows, num_cols, figsize=(15, 25))
axs = axs.ravel()  # Flatten the array of subplots for easier indexing

# Counter for subplot index
subplot_index = 0

    # Iterate over subdirectories
for subdirectory in glob.glob(subdirectories_dir):
    # Check if the subdirectory name is in the 'name' list and if it's within the Galaxies directory
    if os.path.basename(subdirectory.rstrip('/')) in name and '/its/home/ecb41/4MOST/files2/Galaxies/' in subdirectory:
   
    # Find all tar files in the subdirectory
        tar_files = glob.glob(os.path.join(subdirectory, '*.tar.gz'))
    
    # Iterate over tar files
        for tar_file_path in tar_files:
            # Open the tar file
            with tarfile.open(tar_file_path, 'r:gz') as tar:
                # Extract the first FITS file with 'LJ1' in the member name
                fits_member = next((member for member in tar.getmembers() if member.name.endswith('.fits') and 'LJ1' in member.name), None)
                if fits_member:
                    #fits_member = fits_members[0]  # Take the first FITS file
                    fits_file = tar.extractfile(fits_member)
                    with fits.open(fits_file) as hdul:
                        flux = hdul[1].data['FLUX'].flatten()
                        wavelength = hdul[1].data['WAVE'].flatten()
                        folder_name = os.path.basename(subdirectory.rstrip('/'))
                        plot_flux_vs_wavelength(axs[subplot_index], wavelength, flux, label=folder_name)
                        if subplot_index == 0:
                            axs[subplot_index].set_ylim(-0.3e-15, 0.3e-15)
                        if subplot_index == 2:
                            axs[subplot_index].set_ylim(-0.25e-15, 0.45e-15)
                        axs[subplot_index].legend(fontsize='large')   # Display legend on each subplot
                        subplot_index += 1
                        if subplot_index >= num_rows * num_cols:
                            break
        if subplot_index >= num_rows * num_cols:
            break  # Break outer loop if all subplots are plotted

            
# Directory containing the subdirectories
subdirectories_dir2 = '/its/home/ecb41/4MOST/files2/Stars/*'
    # Iterate over subdirectories
for subdirectory in glob.glob(subdirectories_dir2):
    # Check if the subdirectory name is in the 'name' list and if it's within the Galaxies directory
    if os.path.basename(subdirectory.rstrip('/')) in name and '/its/home/ecb41/4MOST/files2/Stars/' in subdirectory:
   
    # Find all tar files in the subdirectory
        tar_files = glob.glob(os.path.join(subdirectory, '*.tar.gz'))
    
    # Iterate over tar files
        for tar_file_path in tar_files:
            # Open the tar file
            with tarfile.open(tar_file_path, 'r:gz') as tar:
                # Extract the first FITS file with 'LJ1' in the member name
                fits_member = next((member for member in tar.getmembers() if member.name.endswith('.fits') and 'LJ1' in member.name), None)
                if fits_member:
                    #fits_member = fits_members[0]  # Take the first FITS file
                    fits_file = tar.extractfile(fits_member)
                    with fits.open(fits_file) as hdul:
                        flux = hdul[1].data['FLUX'].flatten()
                        wavelength = hdul[1].data['WAVE'].flatten()
                        folder_name = os.path.basename(subdirectory.rstrip('/'))
                        plot_flux_vs_wavelength(axs[subplot_index], wavelength, flux, label=folder_name)
                        axs[subplot_index].legend(fontsize='large')   # Display legend on each subplot
                        subplot_index += 1
                        if subplot_index >= num_rows * num_cols:
                            break
        if subplot_index >= num_rows * num_cols:
            break  # Break outer loop if all subplots are plotted

# Adjust layout
plt.tight_layout()
plt.savefig('figures/examplespectra')
plt.show()


# In[9]:


# Counter for subplot index
subplot_index = 0

# Define the path to the tar file
tar_file_path = '/its/home/ecb41/4MOST/files2/Stars/F/FGKTeff_6500-6700_singlespec.tar.gz'

# Open the tar file
with tarfile.open(tar_file_path, 'r:gz') as tar:
    # Extract the first FITS file with 'LJ1' in the member name
    fits_member = next((member for member in tar.getmembers() if member.name.endswith('.fits') and 'LJ1' in member.name), None)
    if fits_member:
        # Extract filename
        fits_filename = os.path.basename(fits_member.name)
        
        # Print the .fits filename
        print(f"FITS Filename for 'F' category: {fits_filename}")


# In[ ]:




