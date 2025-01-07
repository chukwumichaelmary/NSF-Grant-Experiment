import struct
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection 
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from PyQt5.QtCore import pyqtRemoveInputHook, QCoreApplication
from adjustText import adjust_text
import shapefile
import warnings
import tkinter as tk
from tkinter import simpledialog
import os


# The output file will be named transects_new.asc.

run_folder = 'historical_NP32'
run_name   = 'historical_NP32'

#-----------------------------------------------

nTR=0

tkROOT = tk.Tk()
tkROOT.withdraw()

warnings.filterwarnings("ignore", 'Creating legend with loc="best" can be slow with large amounts of data.')

plt.ion()
plt.close('all')
pyqtRemoveInputHook()

# transects_new.asc is the working "new" transects file; it can be replaced or appended to.
fname='transect_inputs/{}/transects_new.asc'.format(run_folder)
if os.path.exists(fname):
    ans=input('delete transects_new.asc? (type "n" to keep and append new transects; anything else will delete): ')
    if ans!='n':
        print('deleted.')
        os.remove(fname)
        with open(fname,'w') as f1:
            f1.write('0\n')

print('generating base plot, be patient...')    


with open('transect_inputs/{}/{}.hyd'.format(run_folder,run_name)) as f1:
    lines=f1.readlines()
for line in lines:
    if line.startswith('number-water-quality-segments-per-layer'):
        NOSEGL=int(line.split()[1])
    elif line.startswith('number-horizontal-exchanges'):
        NOQ1=int(line.split()[1])
print('NOSEGL = {}'.format(NOSEGL))
print('NOQ1 = {}'.format(NOQ1))



class Click():
    def __init__(self, ax, func):
        self.ax=ax
        self.func=func
        self.press=False
        self.move = False
        self.c1=self.ax.figure.canvas.mpl_connect('button_press_event', self.onpress)
        self.c2=self.ax.figure.canvas.mpl_connect('button_release_event', self.onrelease)
        self.c3=self.ax.figure.canvas.mpl_connect('motion_notify_event', self.onmove)
        
    def onclick(self,event):
        if event.inaxes == self.ax:
            self.func(event, self.ax)
    def onpress(self,event):
        self.press=True
    def onmove(self,event):
        if self.press:
            self.move=True
    def onrelease(self,event):
        if self.press and not self.move:
            self.onclick(event)
        self.press=False; self.move=False


f1=open('transect_inputs/{}/{}.poi'.format(run_folder,run_name),'rb')
alldata=list(struct.iter_unpack('i',f1.read(4*NOQ1*4))) #4-byte  integers * NOQ1 rows * 4 columns; each entry is a 1-element tuple
f1.close()

#get indices of open-boundary segments
tmp0=np.array(alldata).reshape((NOQ1,4))[:,:2]  #array of [from,to] cell pairs that define exchange surfaces and positive flow direction
ind=tmp0[:,0]<0   #first column is negative for boundary surfaces (positive flow is always into domain)
tmp0=tmp0[ind,1]
tmp0=[np.mod(val,NOSEGL) for val in tmp0] #collapse all vertical cells into a horizontal-only index
obind=[NOSEGL if val==0 else val for val in tmp0] 

tmp=[np.mod(tup[0],NOSEGL) if tup[0]>0 else -999 for tup in alldata] # -999 for open boundaries
tmp=[NOSEGL if val==0 else val for val in tmp]
poi = np.array(tmp).reshape((NOQ1,4))[:,:2]  #matrix with from,to horizontal cell-index pairs as rows


f2=nc.Dataset('transect_inputs/{}/{}_waqgeom.nc'.format(run_folder,run_name))
tmp=f2.variables['FlowElemContour_x'][:].tolist()
xv=[list(filter(lambda x: x!=-999,l0)) for l0 in tmp]
tmp=f2.variables['FlowElemContour_y'][:].tolist()
yv=[list(filter(lambda x: x!=-999,l0)) for l0 in tmp]
verts=[list(zip(x,y)) for x,y in list(zip(xv,yv))]
f2.close()

#open "old" transect file, if it exists, for plotting
tfile = 'transect_inputs/{}/transects_new.asc'.format(run_folder)
if os.path.exists(tfile):
    with open(tfile,'r') as f1:
        nxs=int(f1.readline().split(';')[0])  # number of transects
        xsnames=[]
        xspts=[]
        for ixs in range(nxs):
            line = f1.readline()
            splt = line.split()
            xsnames.append(splt[0].strip('"'))
            npts = int(splt[2])
            ntot = 0
            pts = []
            while ntot<npts:
                line = f1.readline()
                pts.extend([int(val) for val in line.split()])  #transects are composed of indices into the poi file. negative indices reverse the positive flow direction. the poi file contains pairs of from,to indices (plus 2 other indices in each row that are ignored), where each of these is an index into a cell in the waqgeom file. When these indices in the poi file are greater than NOSEGL, they correspond to vertical layers (N<=NOSEGL is top layer, N+NOSEGL is 2nd layer, etc). If the index in the poi file is negative, it is a sequential index of a boundary, where there is actually no "from" cell, only a "to" cell (2nd index in the pair).That is, all boundaries have a positive flow direction into the model domain.
                ntot=len(pts)
            xspts.append(pts)


#open waste load points file, if it exists, for plotting
wfile = 'transect_inputs/{}/wasteload_pts.txt'.format(run_folder)
if os.path.exists(wfile):
    with open(wfile,'r') as f1:
        npts=int(f1.readline().split(';')[0])  # number of points
        wlpts=[]
        wlnames=[]
        for ipt in range(npts):
            line = f1.readline()
            splt = line.split()
            wlnames.append(line.split()[3].strip("'"))
            pt=int(splt[0])
            val=np.mod(pt,NOSEGL)
            pt=NOSEGL if val==0 else val
            wlpts.append(pt)


#plot model grid
fig,ax0=plt.subplots(figsize=(10,8)) #(14,12))
ax0.axis([484000,665000,4130000,4300000])
ax0.set_aspect('equal')
ax0.axis('off')
coll = PolyCollection(verts,facecolors=(0.8,0.8,0.8), edgecolors='white', linewidths=0.2)
ax0.add_collection(coll)
ax0.autoscale_view()

#plot wasteload points
if os.path.exists(wfile):
    for ipt,wlpt in enumerate(wlpts):
        x,y=list(zip(*verts[wlpt-1]))
        ax0.fill(x,y,'k', edgecolor='k', linewidth=0.2)
        name = wlnames[ipt]
        #ax0.text(x[0],y[0],name)

#plot old transects
if os.path.exists(tfile):
    #texts=[]
    for ixs in range(nxs):
        pts = xspts[ixs]
        name = xsnames[ixs]
        for ipt in pts:
            if ipt<0:
                ito,ifrom = poi[-ipt-1,:]  #ipt<0 means use abs(ipt) as index but reverse positive flow direction (switch from,to)
            else:
                ifrom,ito = poi[ipt-1,:]
            if ifrom!=-999:  #no "from" cell to paint for boundary cells
                x,y=list(zip(*verts[ifrom-1]))
                ax0.fill(x,y,'lightskyblue', edgecolor='w', linewidth=0.2)
            x,y=list(zip(*verts[ito-1]))
            ax0.fill(x,y,'lightsalmon', edgecolor='w', linewidth=0.2)
        ax0.text(x[0],y[0],name)
        #texts.append(ax0.annotate(name, xy=(x[0], y[0]), xytext=(x[0],y[0]+.3)))
    #adjust_text(texts) #, arrowprops=dict(arrowstyle="-", color='k', lw=0.5))


#plot other map layers
#sf = shapefile.Reader('transect_inputs/shpfiles/delta-subregions_2007-01-04_utm10.shp')
#for shape in sf.shapes():
    #points = shape.points
    #ap=plt.Polygon(points, fill=False, edgecolor='g',label='TMDL zones')
    #ax0.add_patch(ap)
sf = shapefile.Reader('transect_inputs/shpfiles/bcs_and_pumps.shp')
for shape in sf.shapes():
    points = shape.points
    ax0.plot(*zip(*points),'m',linewidth=3.0,label='bcs and pumps')

#raise Exception('stop')

#plot boundary cells
for ipt in obind:
    x,y=list(zip(*verts[ipt-1]))
    ax0.fill(x,y,facecolor='none', edgecolor='r', linewidth=0.8, label='boundary cells')
    
handles, labels = ax0.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax0.legend(by_label.values(), by_label.keys())


onFrom = True
fromlist = []
tolist = []

#create main function that allows creation of new transects
def func(event, ax0):
    global onFrom, nTR  #bool that indicates if the selected point is a "from" (True) or "to" (False) point
    global fromlist, tolist  # storage of volume indices
    global selected_iseg_last
    
    x,y = event.xdata, event.ydata
    b = event.button  # = 1 for left, 3 for right
            
    point = Point((x,y))
    print('x,y= {}, {}'.format(x,y))
    
    selected_iseg=-9999
    print('matching segment...')
    for iseg,segcoords in enumerate(verts):
        polygon = Polygon(segcoords)
        if polygon.contains(point):
            selected_iseg=iseg
            break  #uses current iseg
    print('search complete.')
    
    if selected_iseg==-9999:
        if len(fromlist)==len(tolist) and len(fromlist)>0:  
            
            print('Saving new transects file...')
            #find exchange-surface indices for all layers
            indES = np.array([])
            for ifrom,ito in zip(fromlist,tolist):
                if ifrom==-999:
                    #open-boundary segment
                    ipoiob = np.where((poi == (-999,ito+1)).all(axis=1))[0]
                    if len(ipoiob)==0:
                        raise Exception('Error: no match found in poi for {}, {}'.format(ifrom,ito))
                    elif len(ipoiob)==10:
                        indES=np.hstack((indES,ipoiob+1)) 
                    else:
                        raise Exception('Error: number of matches different than multiple of 10 found in poi for {}, {}'.format(ifrom,ito))
                else:
                    ipoipos = np.where((poi == (ifrom+1,ito+1)).all(axis=1))[0]
                    if len(ipoipos)==0:
                        ipoineg = np.where((poi == (ito+1,ifrom+1)).all(axis=1))[0]
                        if len(ipoineg)==0:
                            raise Exception('Error: no match found in poi for {}, {}'.format(ifrom,ito))
                        elif len(ipoineg)==10:
                            indES=np.hstack((indES,-(ipoineg+1))) 
                        else:
                            raise Exception('Error: number of matches different than multiple of 10 found in poi for {}, {}'.format(ifrom,ito))
                    elif len(ipoipos)==10:
                        indES=np.hstack((indES,ipoipos+1))
                    else:
                        raise Exception('Error: number of matches different than multiple of 10 found in poi for {}, {}'.format(ifrom,ito))

            indES = indES.reshape((int(len(indES)/5),5))
            trname = simpledialog.askstring(title="input",prompt="transect name:")
            if os.path.exists(fname):
                f1=open(fname,'a')
            else:
                f1=open(fname,'w')
            np.savetxt(f1,indES,header='"{:s}"          1      {:4d}'.format(trname,np.size(indES)),fmt='   %8d   %8d   %8d   %8d   %8d',comments='')
            f1.close()
            nTR+=1
            print('Saved new transect to transects_new.asc file.')
            fromlist=[]
            tolist=[]
            onFrom=True
            return
        else:
            print('No exchange pairs defined, or a from, to pair is incomplete. Closing.')
            with open(fname,'r') as f1:
                lines = f1.readlines()
                nTR_orig = int(lines[0].split(';')[0])
                nTR = nTR + nTR_orig
                line=[str(nTR)+'\n']
            with open(fname,'w') as f1:
                f1.writelines(line+lines[1:])
            QCoreApplication.quit()
            plt.close()

    
    if b==1:
        xseg,yseg=list(zip(*verts[selected_iseg]))
        if onFrom:
            print('added "from" segment {}'.format(selected_iseg))
            fromlist.append(selected_iseg)
            onFrom=False
            selected_iseg_last = selected_iseg
            ax0.fill(xseg,yseg,facecolor='blue', edgecolor='white', linewidth=0.2)
        else:
            if selected_iseg == selected_iseg_last:
                if selected_iseg+1 in obind:
                    print('boundary surface selected for segment {}.'.format(selected_iseg))
                    fromlist[-1]=-999
                    tolist.append(selected_iseg)
                    onFrom=True
                    ax0.fill(xseg,yseg,facecolor='purple', edgecolor='black', linewidth=0.2)
                else:
                    print('"from" and "to" segments are the same, and the segment is not a boundary segment. ignoring selection.')
                    return
            else:    
                print('added "to" segment {}'.format(selected_iseg))
                tolist.append(selected_iseg)
                onFrom=True
                ax0.fill(xseg,yseg,facecolor='red', edgecolor='black', linewidth=0.2)
    else:
        indfrom = [i for i,val in enumerate(fromlist) if val==selected_iseg]
        indto   = [i for i,val in enumerate(tolist) if val==selected_iseg]
        indseg = list(set(indfrom).union(set(indto)))
        indseg2 = [i for i in indseg if i<len(tolist)]   #Might need to delete a from seg that doesn't yet have a matching to seg
        delete_segs = list(set([fromlist[i] for i in indseg]).union(set([tolist[i] for i in indseg2])))  #delete from AND to values in each pair if pair has one match to selected_iseg
        [fromlist.pop(i) for i,val in enumerate(fromlist) if i in indseg]
        [tolist.pop(i) for i,val in enumerate(tolist) if i in indseg2]
        for iseg in delete_segs:
            print('deleting segment {} and any paired segments'.format(iseg))
            xseg,yseg=list(zip(*verts[iseg]))
            ax0.fill(xseg,yseg,facecolor=(0.8,0.8,0.8), edgecolor='white', linewidth=0.5)
        for isegf,isegt in zip(fromlist,tolist):
            xsegf,ysegf=list(zip(*verts[isegf]))
            xsegt,ysegt=list(zip(*verts[isegt]))
            if isegf==-999:
                ax0.fill(xsegt,ysegt,facecolor='purple', edgecolor='black', linewidth=0.2)
            else:
                ax0.fill(xsegf,ysegf,facecolor='blue', edgecolor='white', linewidth=0.2)
                ax0.fill(xsegt,ysegt,facecolor='red', edgecolor='black', linewidth=0.2)
        onFrom=True

    ax0.figure.canvas.draw()

print('\nLeft-click to define from, to pairs (blue, red). Click only once and be patient.')
print('Select boundary segment as both from and to segment to choose open boundary surface. Cell will turn purple in this case. Directionality will be into model domain.')
print('Right-click to delete associated from/to pair(s).')
print('Left-click outside of model mesh to finish current transect and save it to transect file.')
print('Left-click outside of model mesh without defining new transect to quit.\n')
print('New transect file is transects_new.asc. Rename as needed (e.g., to "transects_loop1.asc").\n')

click = Click(ax0, func)


