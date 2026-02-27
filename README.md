# Template Cleaner

Template(-based E-to-B leakage) cleaner for CMB B-mode analysis.

## Overview

This package implements the E-mode recycling method for cleaning B-mode maps from template-based E-to-B leakage, commonly encountered in CMB polarization analysis. The method works by identifying the B-mode leakage pattern induced by a mask from T and E modes, and then subtracting the leakage contribution from the observed B-mode map.

## Installation

### Requirements
- Python 3.7+
- numpy
- healpy

### From source

```bash
pip install .
```

Or, for development installation:

```bash
pip install -e .
```

## Quick Start

The simplest way to clean your B-mode map is using the main wrapper function:

```python
from template_cleaner import get_cleanedBmap
import numpy as np

# Your input data
IQU_map = np.array([I, Q, U])  # 3xNpix array
mask = np.array([...])          # 1D array of length Npix

# Clean the B-mode map
cleaned_B = get_cleanedBmap(IQU_map, mask)
```

If you need the fitted leakage parameters:

```python
cleaned_B, beta_0, beta_1 = get_cleanedBmap(IQU_map, mask, return_fit=True)
```

## get_cleanedBmap — Recommended Function

### Description

`get_cleanedBmap()` is the recommended high-level interface for B-mode cleaning. It handles all steps automatically:

1. Initializes the E-mode recycling workspace
2. Computes the T and E mode template
3. Fits the linear relationship between template and observed B-modes
4. Produces the cleaned B-mode map

### Signature

```python
get_cleanedBmap(map_IQU, mask_bin, lmax_sht=None, beta_0=None, beta_1=None, return_fit=False)
```

### Parameters

- **map_IQU** (ndarray): Input IQU sky map as a 3×Npix array containing Intensity, Q, and U Stokes parameters
- **mask_bin** (ndarray): Binary or weighted mask as a 1D array of length Npix
- **lmax_sht** (int, optional): Maximum multipole moment. If None, defaults to 3×nside−1
- **beta_0** (float, optional): Fixed intercept for linear fit. If None, will be computed
- **beta_1** (float, optional): Fixed slope for linear fit. If None, will be computed
- **return_fit** (bool, optional): If True, returns fitted parameters along with cleaned map. Default is False

### Returns

- **B_f** (ndarray): Cleaned B-mode map
- **beta_0, beta_1** (float, float): Fitted linear parameters (only if return_fit=True)

### Example

```python
from template_cleaner import get_cleanedBmap

# Basic usage
cleaned_B = get_cleanedBmap(IQU_map, mask)

# With fitted parameters
cleaned_B, intercept, slope = get_cleanedBmap(IQU_map, mask, return_fit=True)
print(f"Linear fit: B_clean = B_observed - {intercept} - {slope} * B_template")

# With fixed parameters
cleaned_B = get_cleanedBmap(IQU_map, mask, beta_0=0.0, beta_1=1.0)
```

## API Reference

### Emode_recycler

Main class for E-mode recycling workspace setup and B-mode cleaning.

#### Constructor

```python
Emode_recycler(IQU_in, mask_in, lmax_in=None)
```

**Parameters:**
- **IQU_in** (ndarray): Input IQU sky map as a 3×Npix array
- **mask_in** (ndarray): Binary or weighted mask as a 1D array of length Npix
- **lmax_in** (int, optional): Maximum multipole moment. If None, defaults to 3×nside−1

**Attributes:**
- `IQU`: Input IQU map
- `nside`: HEALPix nside parameter
- `lmax`: Maximum multipole moment
- `msk`: Mask array
- `B_c`: Measured B-mode map (computed during initialization)
- `B_t`: Template B-mode map (computed by `compute_template()`)
- `B_f`: Cleaned B-mode map (computed by `clean_Bmap()`)

#### Methods

##### compute_template()

Compute the T and E mode template and corresponding B-mode leakage map.

This method extracts the T and E mode spherical harmonic coefficients from the full TEB alms, transforms them back to map space, and then re-analyzes to identify the B-mode leakage pattern induced by the mask.

```python
recycler.compute_template()
```

##### clean_Bmap(beta_0=None, beta_1=None, return_fit=False)

Clean the B-mode map by subtracting template-based E-to-B leakage.

Fits a linear model between the template B-mode and the measured B-mode, then subtracts the fitted leakage component.

**Parameters:**
- **beta_0** (float, optional): Fixed intercept parameter. If None, will be computed
- **beta_1** (float, optional): Fixed slope parameter. If None, will be computed
- **return_fit** (bool, optional): If True, returns fitted parameters. Default is False

**Returns:**
- tuple or None: (beta_0, beta_1) if return_fit=True, else None. Cleaned map stored in self.B_f

**Example:**
```python
recycler = Emode_recycler(IQU_map, mask)
recycler.compute_template()
beta_0, beta_1 = recycler.clean_Bmap(return_fit=True)
cleaned_B = recycler.B_f
```

### templateclean_Blm

Clean B-mode spherical harmonic coefficients using the E-mode template method.

```python
B_clean = templateclean_Blm(alm_TEB, nside, lmax, mask_bin)
```

**Parameters:**
- **alm_TEB** (ndarray): Input TEB spherical harmonic coefficients as a 3×(lmax+1)×(lmax+1) array
- **nside** (int): HEALPix nside parameter for the map
- **lmax** (int): Maximum multipole moment
- **mask_bin** (ndarray): Binary or weighted mask as a 1D array

**Returns:**
- **B_clean** (ndarray): Cleaned B-mode spherical harmonic coefficients

### get_residual

Compute residual B-mode map by comparing cleaned and observed B-modes.

```python
B_residual = get_residual(recycler, IQU_full, ret_full=False)
```

**Parameters:**
- **recycler** (Emode_recycler): An Emode_recycler instance with computed cleaning
- **IQU_full** (ndarray): Full IQU sky map as a 3×Npix array
- **ret_full** (bool, optional): If True, returns both residual and observed B-modes. Default is False

**Returns:**
- **B_r** (ndarray): Residual B-mode map (cleaned minus observed)
- **B_o** (ndarray): Observed B-mode map (only if ret_full=True)

### get_leakage

Compute E-to-B leakage map from full-sky IQU map using mask.

```python
leakage_map = get_leakage(IQU_full, msk_in)
```

**Parameters:**
- **IQU_full** (ndarray): Full IQU sky map as a 3×Npix array
- **msk_in** (ndarray): Binary or weighted mask as a 1D array

**Returns:**
- **L** (ndarray): E-to-B leakage B-mode map

## Advanced Usage

### Working with the class directly

For more control over the process, you can use the `Emode_recycler` class directly:

```python
from template_cleaner import Emode_recycler

# Initialize
recycler = Emode_recycler(IQU_map, mask, lmax_in=256)

# Compute template (automatically called by clean_Bmap if needed)
recycler.compute_template()

# Access intermediate results
print(recycler.nside)      # HEALPix resolution
print(recycler.B_c.shape)  # Measured B-mode shape
print(recycler.B_t.shape)  # Template B-mode shape

# Clean with fitted parameters
beta_0, beta_1 = recycler.clean_Bmap(return_fit=True)
cleaned_B = recycler.B_f
```

### Using fixed fit parameters

If you have pre-computed or known fit parameters:

```python
recycler = Emode_recycler(IQU_map, mask)
# Use fixed parameters
recycler.clean_Bmap(beta_0=0.5, beta_1=0.95)
cleaned_B = recycler.B_f
```

## Citation

This code implements Liu, H. et al (2019)   
ADS: https://ui.adsabs.harvard.edu/abs/2019JCAP...04..046L/abstract


## License

This program is free software licensed under the GNU General Public License v3.
See LICENSE file for details.

## Contact

For more information, visit: https://github.com/1cosmologist/template_cleaner

Contact: Shamik Ghosh (shamik@ustc.edu.cn)
