import numpy as np
import os
import operator
import pickle
from snapshot import Snapshot


class Run:

    def __init__(self, args):

        # Create new snapshots
        # Create a dictionary of snapshots
        self.snapshot_keys = np.arange(args.nstart, args.numss, args.nskip)
        self.parent_numbers = {k: max([args.nstart, k-args.nskip]) for k in self.snapshot_keys}
        self.snapshots = {k: Snapshot(k, self.parent_numbers[k]) for k in self.snapshot_keys}
        self.empty_snapshots = {k: snapshot for k, snapshot in self.snapshots.iteritems()}

        # Repopulate snapshots that exist in restart file
        if args.restart == 1: self.restart(args.sdir + args.rname)

        # Create snapshots in ascending order, using the previous snapshot
        # parallelize this later
        for k, snapshot in sorted(self.empty_snapshots.items(), key=operator.itemgetter(0)):
            parent_snapshot = self.snapshots[self.parent_numbers[k]]
            snapshot.initialize(args, parent_snapshot)
            # Define innermost period of first snapshot
            self.find_innermost_period(snapshot)
            snapshot.set_unitless_time(self.innermost_period)
            self.write_restart(args.sdir + args.rname)

        self.create_summary_plots(args.fdir)
        self.write_summary_files()
        self.initialize_nbody_runs()
        self.initialize_thermal_evolution()



    def restart(self, f):

        # Check if there is a restart file and return pickled object
        # Load in old snapshots and update empty snapshot dictionary
        if os.path.isfile(f):
            with open(f, 'rb') as input:
                old_run = pickle.load(input)
                keys_a = set(self.snapshots.keys())
                keys_b = set(old_run.snapshots.keys())
                intersection = keys_a & keys_b
                for k in intersection:
                    self.snapshots[k] = old_run.snapshots[k]
                    self.empty_snapshots.pop(k, None)
        else:
            print 'restart file missing, proceeding from beginning'
            return None



    def write_restart(self, f):

        with open(f, 'wb') as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)



    def find_innermost_period(self, snapshot):

        if not hasattr(self, 'innermost_period'):
            for planet in snapshot.planets:
                if not hasattr(self, 'innermost_period'):
                    self.innermost_period = planet.aei['a']**1.5
                else:
                    if planet.aei['a']**1.5 < self.innermost_period:
                        self.innermost_period = planet.aei['a']**1.5



    def create_summary_plots(self, fdir):

        #planets_snapshots = [[snapshot.planets[i] for snapshot in snapshots] for i in range(len(snapshots[self.snapshot_keys[0]]))]

        plots.orbit_plots(fdir, self.snapshots)
        plots.mass_plots()






