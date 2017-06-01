import numpy as np
import math
import os
import shutil
import subprocess


class NBody:

    def __init__(self, snapshot):


        g   = 6.67390*math.pow(10,-8)
        ms = 1.98855*math.pow(10.,33.)
        rs  = 6.9599*math.pow(10,10)
        ts  = math.sqrt(math.pow(rs,3)/(g*ms))
        day  = 60.*60.*24.
        au = 1.496*10**13

        # Use final snapshot to generate nbody calculation

        self.radius = [3.0, 3.0]
        self.density = [1000000., 1000000.]

        self.m = [planet.m for planet in snapshot.planets]
        self.a = [planet.aei['a'] for planet in snapshot.planets]
        self.e = [planet.aei['e'] for planet in snapshot.planets]
        self.I = [planet.aei['i']*180./math.pi for planet in snapshot.planets]
        self.g = [planet.aei['g']*180./math.pi for planet in snapshot.planets]
        self.n = [planet.aei['n']*180./math.pi for planet in snapshot.planets]
        self.M = [planet.aei['M']*180./math.pi for planet in snapshot.planets]

        self.x = [planet.r[0]*rs/au for planet in snapshot.planets]
        self.y = [planet.r[1]*rs/au for planet in snapshot.planets]
        self.z = [planet.r[2]*rs/au for planet in snapshot.planets]
        self.u = [planet.v[0]*rs/au/ts*day for planet in snapshot.planets]
        self.v = [planet.v[1]*rs/au/ts*day for planet in snapshot.planets]
        self.w = [planet.v[2]*rs/au/ts*day for planet in snapshot.planets]

        print 'nbody coordinates!'
        print self.x, self.y, self.z
        print self.u, self.v, self.w

        self.ang_x = [0. for planet in snapshot.planets]
        self.ang_y = [0. for planet in snapshot.planets]
        self.ang_z = [0. for planet in snapshot.planets]




    def run_calculation(self, src, dst):

        # Set up run directory

        print src
        print dst

        if not os.path.exists(dst): shutil.copytree(src, dst)

        #self.setup_planets_input(dst)
        self.setup_planets_input_xyz(dst)
#        self.write_pbs_script(dst)

        # Submit pbs file
        cwd = os.getcwd()
        os.chdir(dst)
        subprocess.call("./mercury6")
        os.chdir(cwd)
#        subprocess.call("msub " + dst + "/jobscript.pbs", shell = True)




    def setup_planets_input_xyz(self, dst):

        f = open(dst + '/big.in', 'w')
        f.write(')O+_06 Big-body initial data  (WARNING: Do not delete this line!!)\n')
        f.write(')---------------------------------------------------------------------\n')
        f.write(' style (Cartesian, Asteroidal, Cometary) = Cartesian\n')
        f.write(' epoch (in days) = 0.0\n')
        f.write(')---------------------------------------------------------------------\n')

        for i, m in enumerate(self.m):
            f.write('Kepler'+str(i)+'    m='+str(self.m[i])+' r='+str(self.radius[i])+' d='+str(self.density[i])+'\n')
            f.write(' '+str(self.x[i])+' '+str(self.y[i])+' '+str(self.z[i])+'\n'+' '+str(self.u[i])+' '+str(self.v[i])+' '+str(self.w[i])+'\n')
            f.write(str(self.ang_x[i])+' '+str(self.ang_y[i])+' '+str(self.ang_z[i])+'\n')

        f.close()




    def setup_planets_input(self, dst):

        print self.density, self.m, self.a, self.e

        f = open(dst + '/big.in', 'w')
        f.write(')O+_06 Big-body initial data  (WARNING: Do not delete this line!!)\n')
        f.write(')---------------------------------------------------------------------\n')
        f.write(' style (Cartesian, Asteroidal, Cometary) = Asteroidal\n')
        f.write(' epoch (in days) = 0.0\n')
        f.write(')---------------------------------------------------------------------\n')

        for i, m in enumerate(self.m):
            f.write('Kepler'+str(i)+'    m='+str(self.m[i])+' r='+str(self.radius[i])+' d='+str(self.density[i])+'\n')
            f.write(' '+str(self.a[i])+' '+str(self.e[i])+' '+str(self.I[i])+' '+str(self.g[i])+' '+str(self.n[i])+' '+str(self.M[i])+'\n')
            f.write(str(self.ang_x[i])+' '+str(self.ang_y[i])+' '+str(self.ang_z[i])+'\n')

        f.close()



    def write_pbs_script(self, dst):

        #________________________________________________RUN MERCURY

        # Write pbs file
        f = open(dst + '/jobscript.pbs', 'w')
        f.write('#interpret commands (except PBS directives) using the bash shell\n')
        f.write('#!/bin/bash\n')
        f.write('#Request one nodes using 1 cores per node\n')
        f.write('#PBS -l nodes=1:ppn=1\n')
        f.write('#PBS -l walltime=4:00:00:00\n')
        f.write('#PBS -N Mercury\n')
        f.write('#PBS -j oe\n')
        f.write('#Change the current working directory to the directory from which the script was submitted\n')
        f.write('cd $PBS_O_WORKDIR\n')
        f.write('#These Commands are what would be needed to run this job natively\n')
        f.write('module load intel/2011.3\n')
        f.write('./mercury6')
        f.close()


