import numpy as np
import ipdb
import pickle
import matplotlib.pylab as pl
pl.ion()
datadir = 'data/'
nside = 8192


def measure_ksz_rms_at_sptsz_like_clusters():
    measure_ksz_rms_for_catalog(load_sptsz_like_catalog())

    
def measure_ksz_rms_at_ssdf_like_clusters():
    measure_ksz_rms_for_catalog(load_ssdf_like_catalog())

    
def measure_ksz_rms_for_catalog(d):
    from astropy.io import fits
    tmp = fits.open(datadir+'219_ksz_healpix.fits')[1].data['signal']
    ksz_uk = tmp*jy_per_steradian_to_k_cmb(219)*1e6
    add_healpix_coordinates(d)
    rms = ksz_uk[d['ind_hpix']].std()
    print 'RMS is %0.1f uK-CMB'%(rms)


def add_healpix_coordinates(d):
    import healpy as hp
    phi = d['ra']*np.pi/180.
    theta = (90.-d['dec'])*np.pi/180.
    d['ind_hpix'] = hp.ang2pix(nside, theta, phi, nest=False)

    
def load_sptsz_like_catalog():
    return load_halo_catalog(m500c_min=2e14, z_min=0.1)

def load_ssdf_like_catalog():
    return load_halo_catalog(m500c_min=1.8e14, z_min=1.3)


def load_halo_catalog(m500c_min=-1, m500c_max=9e99, z_min=-1, z_max=100):    
    # see http://lambda.gsfc.nasa.gov/toolbox/tb_sim_readme_halo.cfm for definitions
    # and note that all masses are w.r.t critical densities and have no h's in them.
    import os.path
    if not(os.path.isfile(datadir+'halo_nbody.npy')):
        make_halo_nbody_npy()
    tmp = np.load(datadir+'halo_nbody.npy')
    wh_keep = np.where( (tmp[:,14]>m500c_min) & (tmp[:,14]<m500c_max) & (tmp[:,0]>z_min) & (tmp[:,0]<z_max) )[0]
    tmp = tmp[wh_keep, :]
    keys = ['z','ra','dec',
            'xpos','ypos','zpos',
            'xvel','yvel','zvel',
            'mfof','mvir','rvir',
            'm200c','r200c',
            'm500c','r500c',
            'm1500c','r1500c',
            'm2500c','r2500c']
    d = {}
    for i in range(len(keys)):
        d[keys[i]] = tmp[:, i]
    return d


def make_halo_nbody_npy():
    tmp = np.loadtxt(datadir+'halo_nbody.ascii')
    np.save(datadir+'halo_nbody', tmp)

    
def jy_per_steradian_to_k_cmb(nu):
    tcmb = 2.72548
    return tcmb/{
        30:7.364967e7,
        90:5.526540e8,
        148:1.072480e9,
        219:1.318837e9,
        277:1.182877e9,
        350:8.247628e8}[nu]

    
