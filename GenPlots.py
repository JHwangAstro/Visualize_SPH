import numpy as np
import math
import matplotlib.pyplot as plt

def Double_Plot(x0,x1,y0,y1,mask0,mask1,xlabel,ylabel,fname,log):

  #shell characteristic plots
  plt.plot(np.compress(mask0,x0), np.compress(mask0,y0), color = 'k', linestyle='--')
  plt.plot(np.compress(mask1,x1), np.compress(mask1,y1), color = 'k')
  plt.xlabel(xlabel, fontsize=20)
  plt.ylabel(ylabel, fontsize=20)
  if log == 1: plt.yscale('log')

  plt.savefig(fname + '.ps', facecolor='w', edgecolor='w', format='ps',dpi=200)
  plt.savefig(fname + '.png', facecolor='w', edgecolor='w', format='png',dpi=200)
  plt.clf()




def Particle_Plot(x,y,x1,y1,x2,y2,xlabel,ylabel,fname):

  #shell characteristic plots
  c = ['r','g']
  max_r = 0
  for i in range(2):
    plt.plot(x[i,:],y[i,:],color=c[i],marker='.',markersize=.5,linestyle='')
    max_r = max([max((x[i,:]**2+y[i,:]**2)**0.5),max_r])
  plt.plot(x1,y1,color='b',marker='.',markersize=.5,linestyle='')
  plt.plot(x2,y2,color='y',marker='.',markersize=.5,linestyle='')
  plt.xlabel(xlabel,fontsize=20)
  plt.ylabel(ylabel,fontsize=20)
  plt.axis([-max_r,max_r,-max_r,max_r])

  plt.savefig(fname + '.ps', facecolor='w', edgecolor='w', format='ps',dpi=200)
  plt.savefig(fname + '.png', facecolor='w', edgecolor='w', format='png',dpi=200)
  plt.clf()




def Zoomed_Particle_Plot(x,y,x1,y1,x2,y2,xlabel,ylabel,fname):

  #shell characteristic plots
  c = ['r','g']
  max_r = 0
  for i in range(len(x[:,0])):
    plt.plot(x[i,:],y[i,:],color=c[i],marker='.',markersize=.5,linestyle='')
  plt.plot(x1,y1,color='b',marker='.',markersize=.5,linestyle='')
  plt.plot(x2,y2,color='y',marker='.',markersize=.5,linestyle='')
  plt.xlabel(xlabel,fontsize=20)
  plt.ylabel(ylabel,fontsize=20)

  for i in range(len(x[:,0])):
    dx = max(x[i,:])-min(x[i,:])
    dy = max(x[i,:])-min(x[i,:])
    plt.axis([min(x[i,:])-dx,max(x[i,:])+dx,min(y[i,:])-dy,max(y[i,:])+dy])
    plt.savefig(fname + str(i) + '.ps',  facecolor='w', edgecolor='w', format='ps', dpi=200)
    plt.savefig(fname + str(i) + '.png', facecolor='w', edgecolor='w', format='png',dpi=200)
  plt.clf()















