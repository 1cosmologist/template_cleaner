#######################################################################
#
# Template(-based E-to-B leakage) cleaner for CMB B-mode analysis.
# Copyright (C) 2021  Shamik Ghosh
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# For more information about CMBframe please visit 
# <https://github.com/1cosmologist/template_cleaner> or contact Shamik Ghosh 
# at shamik@ustc.edu.cn
#
#########################################################################
import numpy as np
import healpy as hp 
import os 

skytools_datapath = os.environ['SKYTOOLS_DATA']

class Emode_recycler:
    '''
        E-mode recycling method workspace setup.
    '''

    def __init__(self, IQU_in, mask_in, lmax_in=None):
        """
        Initialize the E-mode recycling workspace.
        
        Parameters
        ----------
        IQU_in : ndarray
            Input IQU sky map as a 3xNpix array containing Intensity, Q, and U Stokes parameters.
        mask_in : ndarray
            Binary or weighted mask as a 1D array of length Npix.
        lmax_in : int, optional
            Maximum multipole moment for spherical harmonic transform. If None, defaults to 3*nside-1.
        """
        self.IQU = IQU_in
        
        self.nside = hp.npix2nside(len(self.IQU[0,:]))
        if lmax_in == None :
            self.lmax = 3*self.nside - 1
        else: 
            self.lmax = lmax_in 

        self.msk = np.array(np.copy(mask_in))
        self.__msk_arr = np.array([self.msk, self.msk, self.msk])

        self.__alm = hp.map2alm(np.copy(self.IQU)*self.__msk_arr, lmax=self.lmax, use_weights=True, datapath=skytools_datapath)
        self.B_c = hp.alm2map(np.copy(self.__alm[2]), self.nside, lmax=self.lmax, pol=False) * self.msk

    def compute_template(self):
        """
        Compute the T and E mode template and corresponding B-mode leakage map.
        
        This method extracts the T and E mode spherical harmonic coefficients from the 
        full TEB alms, transforms them back to map space, and then re-analyzes to identify
        the B-mode leakage pattern induced by the mask.
        """
        TE_lm = np.zeros_like(self.__alm)
        TE_lm[0] = np.copy(self.__alm[0])
        TE_lm[1] = np.copy(self.__alm[1])

        self.IQU_TE = hp.alm2map(TE_lm, self.nside, lmax=self.lmax, pol=True)

        alm_tilde = hp.map2alm(np.copy(self.IQU_TE)*self.__msk_arr, lmax=self.lmax, use_weights=True, datapath=skytools_datapath)

        self.B_t = hp.alm2map(alm_tilde[2], self.nside, lmax=self.lmax, pol=False) * self.msk

    def __lin_fit(self, x, y, intercept_in=None, slope_in=None):
        """
        Perform linear regression fit to data.
        
        Parameters
        ----------
        x : ndarray
            Independent variable data.
        y : ndarray
            Dependent variable data.
        intercept_in : float, optional
            Fixed intercept value. If provided, slope is computed keeping intercept fixed.
        slope_in : float, optional
            Fixed slope value. If provided, intercept is computed keeping slope fixed.
            
        Returns
        -------
        intercept : float
            Y-intercept of the fitted line.
        slope : float
            Slope of the fitted line.
        """
        # Slope parameter:
        if slope_in == None:
            slope = np.cov(x,y)[0,1] / np.var(x, dtype=np.float64)
        else:
            slope = slope_in

        # Intercept paramter:
        if intercept_in == None:
            intercept = np.mean(y, dtype=np.float64) - slope * np.mean(x, dtype=np.float64)
        else:
            intercept = intercept_in
        return intercept, slope

    def clean_Bmap(self, beta_0=None, beta_1=None, return_fit=False):
        """
        Clean the B-mode map by subtracting template-based E-to-B leakage.
        
        Fits a linear model between the template B-mode (B_t) and the measured B-mode (B_c),
        then subtracts the fitted leakage component to produce a cleaned B-mode map.
        
        Parameters
        ----------
        beta_0 : float, optional
            Fixed intercept parameter for the linear fit. If None, will be computed.
        beta_1 : float, optional
            Fixed slope parameter for the linear fit. If None, will be computed.
        return_fit : bool, optional
            If True, returns the fitted parameters (beta_0, beta_1). Default is False.
            
        Returns
        -------
        tuple or None
            If return_fit is True, returns (beta_0, beta_1). Otherwise, returns None.
            The cleaned B-mode map is stored in self.B_f.
        """
        if not hasattr(self, 'B_t'):
            self.compute_template()

        beta_0, beta_1 = self.__lin_fit(np.copy(self.B_t[np.where(self.msk > 0.9)]), np.copy(self.B_c[np.where(self.msk > 0.9)]), intercept_in=beta_0, slope_in=beta_1)

        # print(beta_0,beta_1)

        self.B_f = (self.B_c - beta_0 - (beta_1 * self.B_t))*self.msk 

        if return_fit:
            return beta_0, beta_1

def get_cleanedBmap(map_IQU, mask_bin, lmax_sht=None, beta_0=None, beta_1=None, return_fit=False):
    """
    Compute cleaned B-mode map using the E-mode recycling method.
    
    Wrapper function that instantiates an Emode_recycler object, computes the template,
    and performs B-mode cleaning in a single call.
    
    Parameters
    ----------
    map_IQU : ndarray
        Input IQU sky map as a 3xNpix array.
    mask_bin : ndarray
        Binary or weighted mask as a 1D array of length Npix.
    lmax_sht : int, optional
        Maximum multipole moment. If None, defaults to 3*nside-1.
    beta_0 : float, optional
        Fixed intercept for linear fit. If None, will be computed.
    beta_1 : float, optional
        Fixed slope for linear fit. If None, will be computed.
    return_fit : bool, optional
        If True, returns fitted parameters along with cleaned map. Default is False.
        
    Returns
    -------
    B_f : ndarray
        Cleaned B-mode map.
    beta_0, beta_1 : float, float
        Fitted linear parameters. Only returned if return_fit is True.
    """
    cleaner = Emode_recycler(map_IQU, mask_bin, lmax_in=lmax_sht)
    cleaner.compute_template()

    if return_fit:
        beta_0, beta_1 = cleaner.clean_Bmap(beta_0=beta_0, beta_1=beta_1, return_fit=True)
        return cleaner.B_f, beta_0, beta_1
    else: 
        cleaner.clean_Bmap(beta_0=beta_0, beta_1=beta_1)
        return cleaner.B_f

def __lin_fit(x, y, intercept_in=None, slope_in=None):
    """
    Perform linear regression fit to data (module-level function).
    
    Parameters
    ----------
    x : ndarray
        Independent variable data.
    y : ndarray
        Dependent variable data.
    intercept_in : float, optional
        Fixed intercept value. If provided, slope is computed keeping intercept fixed.
    slope_in : float, optional
        Fixed slope value. If provided, intercept is computed keeping slope fixed.
        
    Returns
    -------
    intercept : float
        Y-intercept of the fitted line.
    slope : float
        Slope of the fitted line.
    """
    # Slope parameter:
    if slope_in == None:
        slope = np.cov(x,y)[0,1] / np.var(x)
    else:
        slope = slope_in

    # Intercept paramter:
    if intercept_in == None:
        intercept = np.mean(y) - slope * np.mean(x)
    else:
        intercept = intercept_in
    return intercept, slope


def templateclean_Blm(alm_TEB, nside, lmax, mask_bin):
    """
    Clean B-mode spherical harmonic coefficients using E-mode template method.
    
    Extracts the T and E modes from the input alm, transforms to map space, applies mask,
    and re-analyzes to identify E-to-B leakage. Performs linear fit to subtract leakage
    from the original B-mode alm.
    
    Parameters
    ----------
    alm_TEB : ndarray
        Input TEB spherical harmonic coefficients as a 3x(lmax+1)x(lmax+1) array.
    nside : int
        HEALPix nside parameter for the map.
    lmax : int
        Maximum multipole moment.
    mask_bin : ndarray
        Binary or weighted mask as a 1D array.
        
    Returns
    -------
    B_clean : ndarray
        Cleaned B-mode spherical harmonic coefficients.
    """
    alm_TE = np.copy(alm_TEB)
    alm_TE[2] = 0.+0.j
    alm_TE = np.ascontiguousarray(alm_TE)
    TE_iqu = hp.alm2map(alm_TE, nside, lmax=lmax)
    leakage_alm = hp.map2alm(TE_iqu*mask_bin, lmax=lmax, use_weights=True, datapath=skytools_datapath)[2]

    del alm_TE, TE_iqu

    beta0, beta1 = __lin_fit(leakage_alm, alm_TEB[2])

    return alm_TEB[2] - beta0 - (beta1 * leakage_alm)

def get_residual(recyler, IQU_full, ret_full=False):
    """
    Compute residual B-mode map by comparing cleaned and observed B-modes.
    
    Uses an Emode_recycler instance that has already been cleaned to compute the difference
    between the cleaned B-mode (B_f) and the observed B-mode from the full IQU map.
    
    Parameters
    ----------
    recyler : Emode_recycler
        An Emode_recycler instance with computed cleaning.
    IQU_full : ndarray
        Full IQU sky map as a 3xNpix array.
    ret_full : bool, optional
        If True, returns both residual and observed B-modes. Default is False.
        
    Returns
    -------
    B_r : ndarray
        Residual B-mode map (cleaned minus observed).
    B_o : ndarray
        Observed B-mode map. Only returned if ret_full is True.
    """
    if not hasattr(recyler, 'B_f'):
        recyler.clean_Bmap()

    alm_o = hp.map2alm(np.copy(IQU_full), lmax=recyler.lmax, use_weights=True, datapath=skytools_datapath)
    
    B_o = hp.alm2map(np.copy(alm_o[2]), recyler.nside, lmax=recyler.lmax, pol=False) * recyler.msk

    B_r = np.copy(recyler.B_f) - B_o
    
    if ret_full :
        return B_r, B_o
    else :
        return B_r

def get_leakage(IQU_full, msk_in):
    """
    Compute E-to-B leakage map from full-sky IQU map using mask.
    
    Extracts T and E modes from the input IQU map, applies the mask, and re-analyzes
    to obtain the B-mode leakage pattern induced by the mask.
    
    Parameters
    ----------
    IQU_full : ndarray
        Full IQU sky map as a 3xNpix array.
    msk_in : ndarray
        Binary or weighted mask as a 1D array.
        
    Returns
    -------
    L : ndarray
        E-to-B leakage B-mode map.
    """
    nside = hp.npix2nside(len(IQU_full[0,:]))
    lmax = 3*nside - 1

    alm_o = hp.map2alm(np.copy(IQU_full), lmax=lmax, use_pixel_weights=True, datapath='/home/shamik/DATA/HPX_pix_wgts')
    full_TE_lm = np.zeros_like(alm_o)
    full_TE_lm[0] = np.copy(alm_o[0])
    full_TE_lm[1] = np.copy(alm_o[1])

    IQU_TE = hp.alm2map(full_TE_lm, nside, lmax=lmax, pol=True)

    alm_m = hp.map2alm(np.copy(IQU_TE)*[msk_in, msk_in, msk_in], lmax=lmax, use_pixel_weights=True, datapath='/home/shamik/DATA/HPX_pix_wgts')

    L = hp.alm2map(np.copy(alm_m[2]), nside, lmax=lmax, pol=False) * msk_in

    return L 



        