
def orbit_plot(fdir, snapshots):

    c = ['b','r']

    t = np.array([snapshot.unitless_time for snapshot in snapshots])
    m_bound_gas = np.array([snapshot.bound_gas for snapshot in snapshots])
    m_unbound_gas = np.array([snapshot.unbound_gas for snapshot in snapshots])
    planet_energy = np.array([snapshot.planet_energy for snapshot in snapshots])
    a = [np.array([snapshot.get_planets_orbital_elements('a',i) for snapshot in snapshots]) for i in range(2)]
    e = [np.array([snapshot.get_planets_orbital_elements('e',i) for snapshot in snapshots]) for i in range(2)]
    m = [np.array([snapshot.get_planets_attribute('m',i) for snapshot in snapshots]) for i in range(2)]
    specific_e = [np.array([snapshot.get_planets_attribute('m',i) for snapshot in snapshots]) for i in range(2)]
    specific_h = [np.array([snapshot.get_planets_attribute('m',i) for snapshot in snapshots]) for i in range(2)]

    # Create orbit plots
    for i in range(2):
        am = np.ma.masked_where(a[i] != 0., a[i])
        em = np.ma.masked_where(a[i] != 0., e[i])
        plt.plot(t[i], am, color = c[i], linestyle = '-')  #semi-major axes
        plt.plot(t[i], (1.-em)*am, color = c[i], linestyle = ':')  #periapsis
        plt.plot(t[i], (1.+em)*am, color = c[i], linestyle = ':')  #apoapsis

    plt.xlim([0,max(t)])
    plt.xlabel(r'$t/\mathrm{P}_\mathrm{0}$', fontsize=15)
    plt.ylabel(r'$a\ [\mathrm{AU}]$', fontsize=15)

    fname = fdir + 'orbits.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    plt.clf()


    # Create mass plots
    for i in range(2):
        mm = np.ma.masked_where(m[i] != 0., m[i])
        plt.plot(t, mm-mm[0], color = c[i], linestyle = '-')  #change in mass

    plt.plot(t, m_bound_gas, color = 'k', linestyle = ':')   #mass bound to system
    plt.plot(t, m_unbound_gas, color = 'k', linestyle = '-.')  #mass in ejecta

    plt.axis([0,max(t),min(min(m[0]-m[0][0]),min(m[1]-m[1][0])),max(max(m[0]-m[0][0]),max(m[1]-m[1][0]))])
    plt.xlabel(r'$t/\mathrm{P}_\mathrm{0}$', fontsize=15)
    plt.ylabel(r'$\dot M/\mathrm{M}_\oplus$', fontsize=15)

    fname = fdir + 'mass_transfer.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    plt.clf()

    plt.plot(t, planet_energy, color = 'k', linestyle = '-')
    plt.plot([0,max(t)],[0.,0.],linestyle='--')
    plt.xlim([0,max(t)])
    plt.xlabel(r'$t/\mathrm{P}_\mathrm{0}$', fontsize=15)
    plt.ylabel(r'$E$', fontsize=15)

    fname = fdir + 'planet_energy.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    plt.clf()

    plt.plot(time,E[0,:], color = 'k', linestyle = '--')  #mass in planet 1
    plt.plot(time,E[1,:], color = 'k', linestyle = '-.')  #mass in planet 1
    plt.plot(time,E[0,:]+E[1,:], color = 'k', linestyle = '-')  #mass in planet 1
    plt.xlim([0,max(time)])
    plt.xlabel(r'$t/\mathrm{P}_\mathrm{0}$', fontsize=15)
    plt.ylabel(r'$E$', fontsize=15)

    fname = fdir + 'specific_orbital_energy.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    plt.clf()

    plt.plot(time,Mgas[0,:], color = 'k', linestyle = '--')  #mass in planet 1
    plt.plot(time,Mgas[1,:], color = 'k', linestyle = '-')   #mass in planet 2
    plt.xlim([0,max(time)])
    plt.xlabel(r'$t/\mathrm{P}_\mathrm{0}$', fontsize=15)
    plt.ylabel(r'$M_\mathrm{gas}$', fontsize=15)

    fname = fdir + 'Gas_Mass.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    plt.clf()

    plt.plot(time,h[0,:], color = 'k', linestyle = '--')  #mass in planet 1
    plt.plot(time,h[1,:], color = 'k', linestyle = '-.')  #mass in planet 1
    plt.plot(time,h[0,:]+h[1,:], color = 'k', linestyle = '-')  #mass in planet 1
    plt.xlim([0,max(time)])
    plt.xlabel(r'$t/\mathrm{P}_\mathrm{0}$', fontsize=15)
    plt.ylabel(r'$h$', fontsize=15)

    fname = fdir + 'specific_orbital_angularmomentum.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    plt.clf()



