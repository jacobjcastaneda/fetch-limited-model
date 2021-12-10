"""
@Author: Jacob Castaneda
@Date: 8-16-21
"""

import numpy as np
import numpy.ma as ma
from stompy.grid.unstructured_grid import UnTRIM08Grid
from scipy.interpolate import NearestNDInterpolator, LinearNDInterpolator
import matplotlib.pyplot as plt
from shapely import geometry
import pdb
import copy
import time
import pickle

def array_NDInterp(xf, yf, xy=None, Z=None, type='Linear', interpolator=None):
    """
    xy: ndarray of x,y coordinates
    Z: scalar value at x,y
    xf, yf: array of new xf, yf values to be interpolated on (not a meshgrid)
    type: specify if Nearest Neighbor or Linear interpolation desired
    interpolator: if you define an interpolator above the scope of this function
    then pass it in to avoid re-doing work
    -----
    RETURN: interpolated values in ndarray
    """
    if interpolator:
        interp = interpolator
    else:
        if type != 'Linear':
            interp = NearestNDInterpolator(xy, Z)
        else:
            interp = LinearNDInterpolator(xy, Z)
            
    res = np.zeros(xf.shape[0])
    for i in range(xf.shape[0]):
        res[i] = interp(xf[i], yf[i])
    return res
    
def find_polygon(grd, poly, xl, yl, L, xc, yc):
    lnLst = [geometry.LineString(zip(xl, yl))]
    lnlngth = len([x for x in lnLst[-1].coords])
    while lnlngth > 3:
        mdpt = lnlngth // 2
        lnLst.append(geometry.LineString(lnLst[-1].coords[:mdpt]))
        if poly.contains(lnLst[-1]):
            break
        lnlngth = len([x for x in lnLst[-1].coords])
    result = np.asarray(lnLst[-1].coords)
    return result[:, 0], result[:, 1], result[-1, 0], result[-1, 1]
        
def distance(x1, y1, x2, y2):
    return np.sqrt( (x2-x1)**2 + (y2-y1)**2 )

def fetch_line(xc, yc, L, thetaw, res):
    xe = xc + L*np.cos(thetaw)
    ye = yc + L*np.sin(thetaw)
    return np.linspace(xc, xe, res), np.linspace(yc, ye, res), xe, ye

def get_close_to_bdry(poly, ln, xc, yc, xe, ye, thetaw, res):
    count = 0
    if poly.contains(ln):
        while poly.contains(ln):
            L = distance(xc, yc, xe, ye)*1.02
            xl, yl, xe, ye = fetch_line(xc, yc, L, thetaw, res)
            ln = geometry.LineString(zip(xl, yl))
            #if L < 25:
             #   xl = np.array([xc])
              #  yl = np.array([yc])
               # xe, ye = xc, yc
            count = count + 1
            if count > 1000:
                with open('log.logs', 'a') as dst:
                    dst.write('xc: {}, yc: {}, xe: {}, ye: {}\n'.format(xc, yc, xe, ye))
                break
    else:
        while not poly.contains(ln):
            L = distance(xc, yc, xe, ye)*.98
            xl, yl, xe, ye = fetch_line(xc, yc, L, thetaw, res)
            ln = geometry.LineString(zip(xl, yl))
            count = count + 1
            if count > 1000:
                with open('log.logs', 'a') as dst:
                    dst.write('xc: {}, yc: {}, xe: {}, ye: {}\n'.format(xc, yc, xe, ye))
                    
                break
                
    return xl, yl, xe, ye


def fetch(cell, grd, thetaw, res, d=None, ctrs=None, nodes=None, poly=None, 
          interp=None, calc_geo=False, SHOW=False, VERBOSE=False):
    """
    grd: pass in stompy UnstructuredGrid Object
    d: bathymetry corresponding to grd cell-centers (made for SWAN)
    thetaw: dir- of wind source 0 rad in x+ dir- (+) CCW (Cartesian Convention, takes dir- of source!!)
    res: resolution for the 1D grid to determine fetch
    
    ----
    RETURN: numpy array of the fetch for cell centers
    """
    if calc_geo:
        ctrs = grd.cells_center()

        nodes = grd.nodes['x'][grd.boundary_cycle()]
        poly = geometry.Polygon(nodes)

    L = np.max([np.max(ctrs[:, 0]) - np.min(ctrs[:, 0]),
                np.max(ctrs[:, 1]) - np.min(ctrs[:, 1])])
    if interp is None:
        if d is None:
            assert(False, "Please enter array of depths for cell centers")
        interp = NearestNDInterpolator(ctrs, d)
    
    if VERBOSE:
        print("Begin computing Fetch...")
        start = time.time()
   
    xc, yc = ctrs[cell, :]
        
    xl, yl, xe, ye = fetch_line(xc, yc, L, thetaw, res)       # find first fetch_line. Big Guess
    xl, yl, xe, ye = find_polygon(grd, poly, xl, yl, L, xc, yc)  # reduce line to one within Poly... Halve it
    ln = geometry.LineString(zip(xl, yl))
    xl, yl, xe, ye = get_close_to_bdry(poly, ln, xc, yc, xe, ye, thetaw, res)
       
    dl = array_NDInterp(xl, yl, interpolator=interp) 
    F = distance(xc, yc, xe, ye)
    Df = np.mean(dl)
    if VERBOSE:
        print(time.time() - start)
    if SHOW:
        fig, ax = plt.subplots()
        grd.plot_edges(ax=ax)
        ax.scatter([xc,xe], [yc,ye], c='r', zorder=10)
        ax.plot(xl, yl, 'k--')
        plt.show()
       
    return F, Df

def read_node_depths(fp):
    """
    fp: pointer to file with depth information... assumes format of depths
    passed to SWAN for UNSTRUC grid

    ---
    RETURN: numpy array of the depths at the grid nodes
    """
    with open(fp, 'rt') as src:
        lines = src.readlines()

    for idx, elem in enumerate(lines):
        tmp = elem.rstrip()
        lines[idx] = float(tmp)
    return np.asarray(lines)


def visualize_grid(grd, vals, xlim, ylim):
    fig, ax = plt.subplots()
    grd.plot_cells(centers=True, values=vals)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

def get_mask(xy, xlim, ylim, lighten=(None, 1)):
    """
    xy: ndarray of the x,y coordinates for points passed in
    xlim: ndarray limits on x to filter out ~AOI
    ylim: ndarray limites on y to filter out ~AOI
    lighten: tuple(boolean, iterations) whether to lighten the point 
    density and by how many iterations. Each iteration halves point density
    -------
    RETURN: A mask of True/False for the points passed in (built to handle 
    cell-centers of Unstructured Grid)... False indicates point is undesired, 
    and True indicates point is in AOI
    NOTE: this is opposite of mask for numpy.ma, should therefore set mask=~mask
    for numpy masked array
    """
    TWO = 2
    x1, x2 = xlim
    y1, y2 = ylim
    res = []
    for i in range(xy.shape[0]):
        if (xy[i, 0] > x1 and xy[i, 0] < x2) and (xy[i, 1] > y1 and xy[i, 1] < y2):
            res.append(True)
        else:
            res.append(False)
    # if lighten != None, then make more sparse by factor given
    if lighten[0]:
        for i in range(lighten[1]):
            cnt = 0
            for j in range(len(res)):
                if res[j]:
                    if cnt%TWO == 0:
                        res[j] = False
                    cnt = cnt + 1
    return np.asarray(res)
       
    

"""
n=1000 works fairly well for resolution here where the max Length is ~O(100,000)
... result is about +/- 100 m on a grid with O(50m-100m) resolution... seems fair 
for fetch computation...
"""
if __name__ == "__main__":
    start = time.time()
    print("Running in Main...")
    # compute for wind source from West, with n=1000 grid resolution
    THETAW = 144*np.pi/180 # mean wind dir- July 2018 
    RES = 100
    
    fpgrd = 'input/sfei_v20_net_mod.grd'
    fpbathy = 'input/sfei_v20_net_mod_depth'

    grd = UnTRIM08Grid(fpgrd)
    ctrs = grd.cells_center()
    
    d = array_NDInterp(ctrs[:, 0], ctrs[:, 1], xy=grd.nodes['x'],
                       Z=read_node_depths(fpbathy), type='Nearest')
    
    #pdb.set_trace()
    interp = NearestNDInterpolator(ctrs, d)
    # P2:573131.50, 4160346.95 // P3: 567858.71, 4160104.62
    
    cell = grd.point_to_cell(np.array([573131.50, 4160346.95]))
    F, Df = fetch(cell, grd, THETAW, RES, interp=interp)


    