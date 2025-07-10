# %% Initialization Cell

import time

print("Loading libraries...") #end and flush removed
t1 = time.time()

import sys, os, gc, copy
from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
#import pymp
import matplotlib as mpl
mpl.use('Agg') #for using mpl over ssh
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath("."))
from modules.user_tools import task, save_table_to_txt, create_gif
import modules.cubehelix as cubehelix

#%matplotlib inline
mpl.rcParams['figure.figsize'] = (6,4)
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12
mpl.rcParams['font.size'] = 12
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.unicode']=True
mpl.rcParams['figure.dpi'] = 100 #dpi when displaying figure
mpl.rcParams['savefig.dpi'] = 100 #dpi when calling np.save

t2 = time.time()
print("Libraries loaded in "+str( round(t2-t1, 3) )+"s.")

#clear any unused variables in RAM
gc.collect()

print("Working directory:", os.getcwd())

# %%

create_assembly_histories = False
make_fit_plots = False
plot_beta = False
median_bin_beta = True
tN_plots = False
make_gif = False

def time_to_snap_code(time):
    snap_code = str(np.int( np.round( time * 4096/13.8 ) ))
    return (4-len(snap_code))*"0" + snap_code

def snap_code_to_time(snap_code):
    return np.int(snap_code)*(13.8/4096.)

# simulations = task(
#     np.load,
#     start_text="Loading reduced halo data",
#     end_text="Loaded reduced halo data",
#     fail_text="Failed to load reduced halo data",
#     exit_on_fail=True
# )("reduced_time_series_data9.npy", encoding="latin1", allow_pickle=True).item()

if create_assembly_histories:
    stat_data = task(
        np.load,
        start_text="Loading M200 stat data",
        end_text="Loaded M200 stat data",
        fail_text="Failed to load M200 stat data",
        exit_on_fail=True
    )("consolidated_M200.npy", encoding="latin1", allow_pickle=True).item()

    DCJL_sims = ["h148","h229","h242","h329"]
    Marvel_sims = ["cptmarvel","rogue","elektra","storm"]
    sims_to_include = Marvel_sims

    path = '/home/christenc/Data/Sims/'
    trace_file_paths = {
        "cptmarvel": path + 'cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.M200.trace_back.hdf5'
        ,"rogue": path + 'rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096.M200.trace_back.hdf5'
        ,"elektra": path + 'elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096.M200.trace_back.hdf5'
        ,"storm": path + 'storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096.M200.trace_back.hdf5'
        ,"h148": path + 'h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/h148.cosmo50PLK.3072g3HbwK1BH.004096.M200.trace_back.hdf5'
        ,"h229": path + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/h229.cosmo50PLK.3072gst5HbwK1BH.004096.M200.trace_back.hdf5'
        ,"h242": path + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/h242.cosmo50PLK.3072gst5HbwK1BH.004096.M200.trace_back.hdf5'
        ,"h329": path + 'h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/h329.cosmo50PLK.3072gst5HbwK1BH.004096.M200.trace_back.hdf5'
    }

    times = []
    volumes = []
    halo_grps = []
    assembly_histories = []
    stellar_mass_histories = []
    Rvir_histories = []
    masses = []
    tracing_fractions = []
    times_per_sim = { sim: [] for sim in sims_to_include }

    for sim in sims_to_include:
        print(f"Getting Mstar and Mpeak for {sim}...")
        sim_trace = pd.read_hdf( trace_file_paths[sim] )

        # satellites_sim = satellites[sim]
        # v_band_mags_sim = v_band_mags[sim]

        halo_grps_sim = list(range(1000)) #list(v_band_mags[sim].keys())

        available_snapshots = list(np.sort(list(stat_data[sim].keys()))[::-1])
        traced_snapshots = list(sim_trace.columns)

        snapshots = [ snapshot for snapshot in traced_snapshots if snapshot[2:] in available_snapshots ]
        all_snapshots = [ "4096" ] + [ snapshot[2:] for snapshot in snapshots ]

        # times = [] #np.vectorize(snap_code_to_time)(snapshots)


        for halo_grp in halo_grps_sim:
            try:
                halo_grp_trace = [halo_grp] + list(sim_trace.loc[halo_grp, snapshots])
            except:
                continue

            # print(all_snapshots)
            # print(snapshots)
            # print(halo_grp_trace)
            # sys.exit()

            # Now get the grp from the earliest snapshot

            final_grp = halo_grp_trace[-1]
            if final_grp > 0 and sim != "elektra":
                if final_grp in inverted_cat.keys():
                    halo_grp_trace.append( int(np.min(inverted_cat[final_grp])) )
                    snaps.append( final_snapshots[sim] )

            # Now build Mstar and Mvir histories from stat_data

            time = []
            Mvir = []
            Mstar = []
            Rvir = []
            tracing_length = 0
            for i in range(len(all_snapshots)):
                grp = halo_grp_trace[i]
                if grp == -1: continue
                snapshot = all_snapshots[i]
                snap = stat_data[sim][ snapshot ]
                if grp not in snap.keys():
                    print(f"Halo {grp} is not in stat_data[{sim}][{snapshot}]")
                    continue
                # print(stat_data[sim][ snapshot ])
                Mvir.append( snap[ grp ][0] if grp >= 1 else 0 )
                Mstar.append( snap[ grp ][1] if grp >= 1 else 0 )
                Rvir.append( snap[ grp ][4] if grp >= 1 else 0 )
                time.append( snap_code_to_time(snapshot) )
                # if grp > 1:
                #     d_mw.append( float(stat_data[sim][ snapshot ][ grp ][3]) )
                if grp >= 1:
                    tracing_length += 1
            tracing_fraction = float(tracing_length)/len(all_snapshots)
            # time = np.vectorize(snapshot_to_time)(all_snapshots)

            # print(sim, halo_grp, d_mw)

            # Mmax_index = np.argmax(Mvir)
            # snapshot_Mpeak = all_snapshots[Mmax_index]
            # halo_grp_Mpeak = halo_grp_trace[Mmax_index]
            # Mmax = Mvir[Mmax_index]
            # z0_mass = Mvir[0]
            # stellar_mass_z0 = 0.4*Mstar[0]

            times.append( time[::-1] )
            volumes.append( sim )
            halo_grps.append( halo_grp )
            # times.append( time )
            masses.append( Mvir[-1] )
            assembly_histories.append( Mvir[::-1] )
            stellar_mass_histories.append( Mstar[::-1] )
            Rvir_histories.append( Rvir[::-1] )
            tracing_fractions.append( tracing_fraction )

    # column_titles = [ "Volume", "halo grp @ z=0", "M_halo @ z=0", "M_peak", "M_star @ z=0", "M_star @ M_peak", "V band mag" ]
    #
    # [ "Volume", "halo grp @ z=0", "halo grp @ M_peak", "snapshot @ M_peak", "M_peak", "M_halo @ z=0", "M_star @ z=0", "M_star,photometric @ z=0", "M_star @ M_peak", "M_star,photometric @ M_peak", "V band mag", "Satellite type" ]
    #
    table = {
        "Time": times,
        "Volume": volumes,
        "halo grp @ z=0": halo_grps,
        # "Age @ peak": peak_ages,
        "Assembly History": assembly_histories,
        "Stellar Mass History": stellar_mass_histories,
        "Rvir History": Rvir_histories,
        "M_halo @ z=0": masses,
        "Tracing Fraction": tracing_fractions
    }

    task(
        np.save,
        start_text=f"Saving data",
        end_text=f"Saved data",
        fail_text=f"Failed to save data",
        exit_on_fail=True,
        no_spinner=False
    )("/home/christenc/Code/Datafiles/assembly_histories.npy", table, allow_pickle=True)

    print("Done!")


if make_fit_plots:

    data = task(
        np.load,
        start_text="Loading assembly history data",
        end_text="Loaded assembly history data",
        fail_text="Failed to load assembly history data",
        exit_on_fail=True
    )("plots/calculate_beta/assembly_histories.npy", encoding="latin1", allow_pickle=True).item()



    assembly_histories = task(
        np.load,
        start_text="Loading M200 data",
        end_text="Loaded M200 data",
        fail_text="Failed to load M200 data",
        exit_on_fail=True
    )("plots/calculate_beta/M200_reduced_data.npy", encoding="latin1", allow_pickle=True).item()

    reduced_data = task(
        np.load,
        start_text="Loading reduced data",
        end_text="Loaded reduced data",
        fail_text="Failed to load reduced data",
        exit_on_fail=True
    )("reduced_time_series_data10.npy", encoding="latin1", allow_pickle=True).item()

    def beta_fit(t, log_m_0, beta):
        return (10.**log_m_0)*(t**beta)

    new_volumes, new_halo_grps, new_masses, peak_masses, times_at_stop, stellar_masses, stellar_mass_fractions, betas, m_0s, stds, sigma_betas, sigma_m_0s, fit_to_entire_history, halos_for_1e8_too_soon, rvirs = ( [] for _ in range(15) )

    times, volumes, halo_grps, stellar_mass_histories, Rvir_histories, masses, tracing_fractions = ( data[key] for key in ["Time","Volume","halo grp @ z=0","Stellar Mass History","Rvir History","M_halo @ z=0","Tracing Fraction"] )

    crossing_times = { f"t_{n}": [] for n in range(1,100,1) }

    number_of_halos = len(volumes)

    Marvel_volumes = ["elektra","rogue","storm","cptmarvel"]
    for i in range(number_of_halos):
        time, volume, halo_grp, stellar_mass_history, Rvir_history, tracing_fraction = ( item[i] for item in ( times, volumes, halo_grps, stellar_mass_histories, Rvir_histories, tracing_fractions ) )

        if volume not in Marvel_volumes: continue

        assembly_history = np.array( assembly_histories[volume][halo_grp]["M200 virial mass"]["data"] )[::-1]
        time_plus = np.array( assembly_histories[volume][halo_grp]["M200 virial mass"]["times"] )[::-1]

        if len(assembly_history) < 20:
            print(f"Skipping {volume} halo {halo_grp} due to truncated history.")
            continue

        mass = assembly_history[-1]
        peak_mass = np.max(assembly_history)
        stellar_mass = stellar_mass_history[-1]
        Rvir = Rvir_history[-1]
        stellar_mass_fraction = stellar_mass/mass

        # truncate data at M_peak or where M = 1e8 M_sun, whichever comes first
        peak_argmax = np.argmax(assembly_history)

        # print(assembly_history)

        arg_1e8 = len(assembly_history)-1
        for i in range(len(assembly_history)):
            if assembly_history[i] == np.nan:
                continue
            if assembly_history[i] >= 1.0e8:
                arg_1e8 = i
                break
        # print(f"{volume} {halo_grp} arg_1e8={arg_1e8}, {n} loops")

        arg_stop = peak_argmax if peak_argmax < arg_1e8 else arg_1e8
        fit_to_peak = peak_argmax < arg_1e8

        if peak_argmax < 3:
            print(f"{volume} halo {halo_grp} cannot be fit because M_peak happens within the first two time steps.")
            continue

        cut_off_at_peak = arg_1e8 <= 4
        if cut_off_at_peak:
            arg_stop = peak_argmax
            fit_to_peak = True

        # fteh = arg_1e8 <= 1
        fit_to_entire_history.append( fit_to_peak )
        halos_for_1e8_too_soon.append( cut_off_at_peak )
        # if fteh:
        #     arg_stop = len(assembly_history)-1
            # print(f"Skipping {volume} halo {halo_grp} because M_halo is never <= 1e8 M_sun.")
            # continue


        original_assembly_history = copy.copy(assembly_history)
        original_time = copy.copy(time_plus)

        assembly_history = assembly_history[:arg_stop+1] # +1 to include M_peak
        time_plus = time_plus[:arg_stop+1]

        time_at_stop = time_plus[-1]

        # normalize:
        # log_mass_history = np.log10(assembly_history)
        time_plus = np.array(time_plus)/13.8 # Gyr


        # Fit beta curve to assembly history
        # try:
        popt, pcov = curve_fit(
            beta_fit, time_plus, assembly_history,
            p0 = [ 9, 0.5 ],
            bounds = ([4,0],[12,20])
        )

        # print(volume, halo_grp, popt, len(assembly_history))
        # except:
        #     print(f"{volume} {halo_grp} failed to fit.")
        #     continue

        log_m_0, beta = tuple(popt)
        m_0 = 10.**log_m_0
        sigma_beta = pcov[1][1]**(0.5)
        sigma_m_0 = 10.**(pcov[0][0]**(0.5))

        # print(volume,halo_grp,beta)
        # print(mass_history_time)
        # print(mass_history)

        std = np.std( np.log10(assembly_history) - np.log10(beta_fit(time_plus, np.log10(m_0), beta)) )

        new_volumes.append(volume)
        new_halo_grps.append(halo_grp)
        times_at_stop.append(time_at_stop)
        new_masses.append(mass)
        peak_masses.append(peak_mass)
        stellar_masses.append(stellar_mass)
        stellar_mass_fractions.append(stellar_mass_fraction)
        betas.append(beta)
        m_0s.append(m_0)
        stds.append(std)
        sigma_betas.append(sigma_beta)
        sigma_m_0s.append(sigma_m_0)

        # print(volume, halo_grp, assembly_histories[volume][halo_grp].keys() )
        rvirs.append( reduced_data[volume][halo_grp]["M200 virial radius"]["data"][-1] )

        halo_type = "luminous" if stellar_mass > 0 else "dark"

        # calculate t_n for the assembly history
        # normalized_assembly_history = assembly_history/1.0e8
        for n in range(1,100,1):
            nth_mass_threshold = 1.0e8*(n/100.)
            try:
                arg_t_n = np.argwhere( assembly_history <= nth_mass_threshold )[-1][0]
            except:
                arg_t_n = 0
            t_n = time_plus[arg_t_n]*13.8
            crossing_times[f"t_{n}"].append( t_n )
            # crossing_times[f"t_{n}"]["halo_type"].append( halo_type )

        reason = (' because M=10^8M_sun happens too early' if cut_off_at_peak else (' because the halo never reaches M=10^8M_sun' if fit_to_peak else ''))

        # print(f"Fitting {volume} {halo_grp}: beta={beta}, Fit to: t<{('t_peak' if fit_to_peak else 't(M=10^8M_sun)')}{reason}.")

        continue # Comment out if you want to make the plots

        time_grid = np.linspace(0,1,100)
        mass_fit = beta_fit(time_grid,np.log10(m_0),beta)


        # plot fit
        f, ax = plt.subplots(1, 1, figsize=(6,4))

        ax.plot( original_time, original_assembly_history, 'g', linewidth=2, color="#2244CC", label="M200 Data", zorder=1 )
        ax.plot( time_grid*13.8, mass_fit, 'g', linewidth=2, color="#DD2244", label="Fit", zorder=0 )

        text = "\n".join([
                    r"$log_{10}(M_{halo,z=0}/M_\odot)="+str(np.round(np.log10(mass),3))+r"$",
                    r"$log_{10}(M_{peak}/M_\odot)="+str(np.round(np.log10(peak_mass),3))+r"$",
                    ( r"$log_{10}(M_{star}/M_\odot)="+str(np.round(np.log10(stellar_mass),3))+r"$" if stellar_mass > 0 else r"$M_{star}=0 \ M_\odot$" ),
                    r"$log_{10}(M_{0}/M_\odot)="+str(np.round(np.log10(m_0),3))+r"$",
                    rf"$\beta="+str(np.round(beta,3))+r"$",
                    r"$M=M_0 \left( t/(13.8 Gyr) \right)^{\beta}$",
                    r"Fit to "+( r"$t < t_{peak}$" if fit_to_peak else r"$t < t(M \leq 10^{8} M_{\odot})$" )
                ])
        ax.text( 0.95, 0.05, text, size=8, horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes)

        ax.set_title(f"{volume} {halo_grp}")
        ax.set_xlabel(r"Time [$Gyr$]")
        ax.set_ylabel(r"$M_{halo}$ [$M_\odot$]")

        ax.ticklabel_format(axis="both", style="sci", scilimits=(0,0))

        ax.set_xlim([0,13.8])

        # ax.set_yscale("log")

        ax.set_ylim([0,None])

        # plt.tight_layout()

        task(
            plt.savefig,
            start_text=f"Saving PNG of {volume} {halo_grp}",
            end_text=f"Saved PNG of {volume} {halo_grp}",
            fail_text=f"Failed to save PNG of {volume} {halo_grp}",
            exit_on_fail=False
        )(f"plots/calculate_beta/fits/{volume}/beta_fit_{volume}_{halo_grp}.png")

        plt.close()

    table = {}

    # new_volumes, new_halo_grps, masses, stellar_masses, stellar_mass_fractions, m_0s, betas, sigma_m_0s, sigma_betas, stds

    column_titles = ["Volume","grp @ z=0","M_halo @ z=0 [M_sun]","min( t_peak, t(M=1e8M_sun)) [Gyr]","M_peak [M_sun]","M_star [M_sun]","M_star/M_halo","M_0","beta","M_0 sigma","beta sigma","residual sigma","Fit to t<=t_peak?","crossed 1e8 too early","R_200 @ z=0"]

    mass_argsort = np.argsort(peak_masses)
    for item, column_title in zip([new_volumes, new_halo_grps, new_masses, times_at_stop, peak_masses, stellar_masses, stellar_mass_fractions, m_0s, betas, sigma_m_0s, sigma_betas, stds, fit_to_entire_history, halos_for_1e8_too_soon, rvirs],column_titles):
        table[column_title] = np.array(item)[mass_argsort]

    for n in range(1,100,1):
        crossing_times[f"t_{n}"] = np.array(crossing_times[f"t_{n}"])[mass_argsort]

    task(
        np.save,
        start_text="Saving betas.npy",
        end_text="Saved betas.npy",
        fail_text="Failed to save betas.npy",
        exit_on_fail=True
    )("plots/calculate_beta/betas.npy", table)

    task(
        np.save,
        start_text="Saving crossing_times.npy",
        end_text="Saved crossing_times.npy",
        fail_text="Failed to save crossing_times.npy",
        exit_on_fail=True
    )("plots/calculate_beta/crossing_times.npy", crossing_times)

    task(
        save_table_to_txt,
        start_text="Saving txt file",
        end_text="Saved txt file",
        fail_text="Failed to save txt file",
        exit_on_fail=True
    )( "plots/calculate_beta/betas.txt", table, column_titles=column_titles, column_padding=3 )


if plot_beta:

    table = task(
        np.load,
        start_text="Loading betas.npy",
        end_text="Loaded betas.npy",
        fail_text="Failed to load betas.npy",
        exit_on_fail=True
    )("plots/calculate_beta/betas.npy", encoding="latin1", allow_pickle=True).item()

    extra_data = task(
        np.load,
        start_text="Loading concentrations_and_einasto.npy",
        end_text="Loaded concentrations_and_einasto.npy",
        fail_text="Failed to load concentrations_and_einasto.npy",
        exit_on_fail=True
    )("plots/calculate_beta/concentrations_and_einasto.npy", encoding="latin1", allow_pickle=True).item()

    #["Volume","grp @ z=0","M_halo @ z=0 [M_sun]","min( t_peak, t(M=1e8M_sun)) [Gyr]","M_peak [Gyr]","M_star [M_sun]","M_star/M_halo","M_0","beta","M_0 sigma","beta sigma","residual sigma","Fit to t<t_peak?","crossed 1e8 too early"]

    volumes = table["Volume"]
    halo_grps = table["grp @ z=0"]
    masses = table["M_halo @ z=0 [M_sun]"]
    peak_masses = table["M_peak [M_sun]"]
    stellar_masses = table["M_star [M_sun]"]
    peak_ages = table["min( t_peak, t(M=1e8M_sun)) [Gyr]"]
    betas = table["beta"]
    stds = table["residual sigma"]
    rvirs = table["R_200 @ z=0"]
    #     ,"beta": betas
    #     ,"sigma": stds

    cmap = cubehelix.cmap(reverse=True, start=(3/4.)*np.pi, rot=0.65, maxLight=0.8, minLight=0.2, minSat=1.3, maxSat=1.3)
    # cmap.set_bad(color="#CCCCCC")
    concentrations = [ r_200/extra_data[volume][halo_grp]["scale radius"]["r_s"] for volume, halo_grp, r_200 in zip(volumes, halo_grps, rvirs) ]
    norm = mpl.colors.LogNorm(vmin=np.min(concentrations), vmax=np.max(concentrations))

    def con_zord(con):
        return int( (np.abs(np.log10(con) - 1))**(4/5.) + 10 )
    def con_alpha(con):
        return np.min( [0.8*(np.abs(np.log10(con) - 1))**(4/5.) + 0.2, 1.0] )


    f, ax = plt.subplots(1, 1, figsize=(6,4))

    for volume, halo_grp, peak_mass, stellar_mass, age, beta, r_200 in zip(volumes, halo_grps, peak_masses, stellar_masses, peak_ages, betas, rvirs):

        # print(extra_data[volume][halo_grp]["core einasto scale radius"].keys())
        r_s = extra_data[volume][halo_grp]["scale radius"]["r_s"] # scale radius
        # virial_radius = data[]
        concentration = r_200/r_s #extra_data[volume][halo_grp]["concentration"]

        markeredgecolor = "#FFFFFF"


        marker, markersize, markeredgewidth, zorder, alpha = ('*',12,0.6,20,1) if stellar_mass > 0 else ('o',4,0,con_zord(concentration),con_alpha(concentration))

        # print(f"{volume} {halo_grp} c={np.round(concentration,2)}, r_s={np.round(r_s,2)} kpc, r_200={np.round(r_200,2)}")
        # print(f"{norm(concentration)}")

        # print(f"{volume} {halo_grp} {concentration} {r_s} {r_200}")

        ax.plot( peak_mass, beta, marker, markersize=markersize, markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor, color=cmap(norm(concentration)), zorder=zorder, alpha=alpha )

    ax.axhline(y=1,linestyle="--",linewidth=2,color="#999999",alpha=0.5,zorder=40)

    ax.set_xlabel(r"$M_{peak}$ [$M_{\odot}$]")
    ax.set_ylabel(r"$\beta$")

    ax.set_xscale('log')

    ax.plot([],[],'*',markersize=12, markeredgewidth=0, label="Luminous halos")
    ax.plot([],[],'o',markersize=4, markeredgewidth=0, label="Dark halos")
    ax.legend(frameon=False)
    # ax.set_yscale('log')

    # ax.set_xlim([0,-2.1])
    ax.set_ylim([0,2])

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, shrink=1.0, ax=ax, pad=0.04) #, ticks=np.linspace(0,1,5))
    # cbar = plt.colorbar(sm, location="top", shrink=1.0, ax=axs[1], pad=0.04)

    cbar.set_label(r'Concentration $R_{200}/R_{scale}$', rotation=270)
    cbar.ax.get_yaxis().labelpad = 15
    # cbar.ax.invert_yaxis()

    plt.tight_layout()

    task(
        plt.savefig,
        start_text="Saving PNG",
        end_text="Saved PNG",
        fail_text="Failed to save PNG",
        exit_on_fail=True
    )("plots/calculate_beta/beta_vs_peak_mass.png")

    plt.close()


    f, ax = plt.subplots(1, 1, figsize=(6,4))

    for volume, halo_grp, std, stellar_mass, age, beta, r_200 in zip(volumes, halo_grps, stds, stellar_masses, peak_ages, betas, rvirs):

        r_s = extra_data[volume][halo_grp]["scale radius"]["r_s"] # scale radius

        concentration = r_200/r_s #extra_data[volume][halo_grp]["concentration"]
        # print(f"{volume} {halo_grp} c={np.round(concentration,2)}, r_s={np.round(r_s,2)} kpc, r_200={np.round(r_200,2)}")

        markeredgecolor = "#FFFFFF"
        marker, markersize, markeredgewidth, zorder, alpha = ('*',12,0.6,20,1) if stellar_mass > 0 else ('o',4,0,con_zord(concentration),con_alpha(concentration))
        ax.plot( std, beta, marker, markersize=markersize, markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor, color=cmap(norm(concentration)), zorder=zorder, alpha=alpha )

    ax.axhline(y=1,linestyle="--",linewidth=2,color="#999999",alpha=0.5,zorder=40)

    ax.set_xlabel(r"Residual $\sigma$")
    ax.set_ylabel(r"$\beta$")

    ax.set_xscale('log')
    # ax.set_yscale('log')

    ax.set_xlim([1e-3,1e0])
    ax.set_ylim([0,2])

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, shrink=1.0, ax=ax, pad=0.04) #, ticks=np.linspace(0,1,5))
    # cbar = plt.colorbar(sm, location="top", shrink=1.0, ax=axs[1], pad=0.04)

    cbar.set_label(r'Concentration $c=log_{10}(R_{200}/R_{scale})$', rotation=270)
    cbar.ax.get_yaxis().labelpad = 15
    # cbar.ax.invert_yaxis()

    ax.plot([],[],'o',markersize=2, markeredgewidth=0, label="Dark halos")
    ax.plot([],[],'*',markersize=12, markeredgewidth=0, label="Luminous halos")
    ax.legend(frameon=False)

    plt.tight_layout()

    task(
        plt.savefig,
        start_text="Saving PNG",
        end_text="Saved PNG",
        fail_text="Failed to save PNG",
        exit_on_fail=True
    )("plots/calculate_beta/beta_vs_std.png")

    plt.close()

    print("Done!")


if median_bin_beta:

    table = task(
        np.load,
        start_text="Loading betas.npy",
        end_text="Loaded betas.npy",
        fail_text="Failed to load betas.npy",
        exit_on_fail=True
    )("plots/calculate_beta/betas.npy", encoding="latin1", allow_pickle=True).item()

    extra_data = task(
        np.load,
        start_text="Loading concentrations_and_einasto.npy",
        end_text="Loaded concentrations_and_einasto.npy",
        fail_text="Failed to load concentrations_and_einasto.npy",
        exit_on_fail=True
    )("plots/calculate_beta/concentrations_and_einasto.npy", encoding="latin1", allow_pickle=True).item()

    volumes =  np.array(table["Volume"])
    halo_grps = np.array(table["grp @ z=0"])
    masses = np.array(table["M_halo @ z=0 [M_sun]"])
    peak_masses = np.array(table["M_peak [M_sun]"])
    stellar_masses = np.array(table["M_star [M_sun]"])
    peak_ages = table["min( t_peak, t(M=1e8M_sun)) [Gyr]"]
    betas = np.array(table["beta"])
    stds = table["residual sigma"]
    rvirs = np.array(table["R_200 @ z=0"])

    concentrations = np.array([ r_200/extra_data[volume][halo_grp]["scale radius"]["r_s"] for volume, halo_grp, r_200 in zip(volumes, halo_grps, rvirs) ])

    luminous_mass_bin_centers = []
    median_luminous_betas = []
    median_luminous_betas_25 = []
    median_luminous_betas_75 = []
    median_luminous_concentrations = []

    dark_mass_bin_centers = []
    median_dark_betas = []
    median_dark_betas_25 = []
    median_dark_betas_75 = []
    median_dark_concentrations = []

    mass_bin_edges = np.logspace(7, 11, 14)

    N=0
    for left_bin_edge, right_bin_edge in zip(mass_bin_edges[:-1],mass_bin_edges[1:]):
        mass_bin_args = np.argwhere( np.logical_and(peak_masses > left_bin_edge, peak_masses < right_bin_edge) )
        # print(len(masses > left_bin_edge), len([ item for item in (masses > left_bin_edge) if item ]))
        mass_bin_args = np.array([ int(mass_bin_args[i][0]) for i in range(len(mass_bin_args)) ])

        if len(mass_bin_args) == 0:
            continue

        mass_bin = peak_masses[ mass_bin_args ]
        stellar_mass_bin = stellar_masses[ mass_bin_args ]
        beta_mass_bin = betas[ mass_bin_args ]
        con_mass_bin = concentrations[ mass_bin_args ]

        mass_bin_center = 10.**( (np.log10(left_bin_edge) + np.log10(right_bin_edge))/2.0 )

        luminous_halo_args = np.argwhere( stellar_mass_bin > 0 )
        if len(luminous_halo_args) > 1:

            luminous_halo_args = np.array([ int(luminous_halo_args[i][0]) for i in range(len(luminous_halo_args)) ])
            luminous_mass_bin = mass_bin[ luminous_halo_args ]
            luminous_betas_mass_bin = beta_mass_bin[ luminous_halo_args ]
            luminous_con_mass_bin = con_mass_bin[ luminous_halo_args ]

            median_luminous_betas.append( np.median(luminous_betas_mass_bin) )
            median_luminous_betas_25.append( np.percentile( luminous_betas_mass_bin, 25 ) )
            median_luminous_betas_75.append( np.percentile( luminous_betas_mass_bin, 75 ) )
            median_luminous_concentrations.append( np.median(luminous_con_mass_bin) )

            luminous_mass_bin_centers.append( mass_bin_center )


        dark_halo_args = np.argwhere( stellar_mass_bin == 0 )
        if len(dark_halo_args) > 1:
            dark_halo_args = np.array([ int(dark_halo_args[i][0]) for i in range(len(dark_halo_args)) ])
            dark_mass_bin = mass_bin[ dark_halo_args ]
            dark_betas_mass_bin = beta_mass_bin[ dark_halo_args ]
            dark_con_mass_bin = con_mass_bin[ dark_halo_args ]

            median_dark_betas.append( np.median(dark_betas_mass_bin) )
            median_dark_betas_25.append( np.percentile( dark_betas_mass_bin, 25 ) )
            median_dark_betas_75.append( np.percentile( dark_betas_mass_bin, 75 ) )
            median_dark_concentrations.append( np.median(dark_con_mass_bin) )

            dark_mass_bin_centers.append( mass_bin_center )

        N += 1

    print("N=",N)

    plotted_concentrations = median_dark_concentrations + median_luminous_concentrations

    cmap = cubehelix.cmap(reverse=True, start=(3/4.)*np.pi, rot=0.65, maxLight=0.8, minLight=0.2, minSat=1.3, maxSat=1.3)
    norm = mpl.colors.LogNorm(vmin=2e0, vmax=2e1) #mpl.colors.LogNorm(vmin=np.min(plotted_concentrations), vmax=np.max(plotted_concentrations))

    f, ax = plt.subplots(1, 1, figsize=(6,4))

    # Plot the individual halos in the background:
    for peak_mass, stellar_mass, age, beta in zip(peak_masses, stellar_masses, peak_ages, betas):
        marker, markersize, markeredgewidth, zorder, alpha, color = ('*',9,0.6,20,1,"#666666") if stellar_mass > 0 else ('o',4,0,15-age,0.6,"#BBBBBB")
        ax.plot( peak_mass, beta, marker, markersize=markersize, markeredgewidth=0, color=color, zorder=zorder, alpha=0.3 )

    # # Now plot the medians:
    # ax.fill_between(
    #     luminous_mass_bin_centers, # x values
    #     median_luminous_betas_25, # y lower bound
    #     median_luminous_betas_75, # y upper bound
    #     facecolor="#FF6666",
    #     linewidth=0,
    #     zorder=30,
    #     alpha=0.3
    # )
    #
    # # Now plot the medians:
    # ax.fill_between(
    #     dark_mass_bin_centers, # x values
    #     median_dark_betas_25, # y lower bound
    #     median_dark_betas_75, # y upper bound
    #     facecolor="#0066FF",
    #     linewidth=0,
    #     zorder=29,
    #     alpha=0.3
    # )

    ax.axhline(y=1,linestyle="--",linewidth=2,color="#666666",zorder=31)

    ax.plot( dark_mass_bin_centers, median_dark_betas, 'g', linestyle="--", linewidth=1.5, color="#0066FF", zorder=31, alpha=1, label=r"Median $\beta$, dark halos" )
    for i in range(len(dark_mass_bin_centers)):
        ax.plot( dark_mass_bin_centers[i], median_dark_betas[i], 'o', color=cmap(norm(median_dark_concentrations[i])), markersize=6, zorder=32 )

    ax.plot( luminous_mass_bin_centers, median_luminous_betas, 'g', linestyle="--", linewidth=1.5, color="#FF6666", zorder=33, alpha=1, label=r"Median $\beta$, luminous halos" )
    for i in range(len(luminous_mass_bin_centers)):
        ax.plot( luminous_mass_bin_centers[i], median_luminous_betas[i], 'o', color=cmap(norm(median_luminous_concentrations[i])), markersize=6, zorder=34 )


    ax.set_xlabel(r"$M_{peak}$ [$M_\odot$]")
    ax.set_ylabel(r"$\beta$")

    ax.set_xscale('log')
    # ax.set_yscale('log')

    ax.set_xlim([1e7,1e11])
    ax.set_ylim([0,2])

    # sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    # sm.set_array([])
    # cbar = plt.colorbar(sm, shrink=1.0, ax=ax, pad=0.04) #, ticks=np.linspace(0,1,5))
    # cbar = plt.colorbar(sm, location="top", shrink=1.0, ax=axs[1], pad=0.04)

    # cbar.set_label(r'Age @ $M_{peak}$ [$Gyr$]', rotation=270)
    # cbar.ax.get_yaxis().labelpad = 15
    # cbar.ax.invert_yaxis()

    ax.plot([],[],'o',markersize=4, markeredgewidth=0, color="#BBBBBB", label="Dark halos")
    ax.plot([],[],'*',markersize=9, markeredgewidth=0, color="#666666", label="Luminous halos")
    ax.legend(frameon=False,prop={"size":8})

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, shrink=1.0, ax=ax, pad=0.04) #, ticks=np.linspace(0,1,5))
    # cbar = plt.colorbar(sm, location="top", shrink=1.0, ax=axs[1], pad=0.04)

    cbar.set_label(r'Median concentration $c=R_{200}/R_{scale}$', rotation=270)
    cbar.ax.get_yaxis().labelpad = 15

    plt.tight_layout()

    task(
        plt.savefig,
        start_text="Saving PNG",
        end_text="Saved PNG",
        fail_text="Failed to save PNG",
        exit_on_fail=True
    )("plots/calculate_beta/median_beta_vs_mpeak.png")

    plt.close()


if tN_plots:
    # ["Volume","grp @ z=0","M_halo @ z=0 [M_sun]","min( t_peak, t(M=1e8M_sun)","M_peak [Gyr]","M_star [M_sun]","M_star/M_halo","M_0","beta","M_0 sigma","beta sigma","residual sigma","Fit to t<t_peak?","crossed 1e8 too early"]

    table = task(
        np.load,
        start_text="Loading betas.npy",
        end_text="Loaded betas.npy",
        fail_text="Failed to load betas.npy",
        exit_on_fail=True
    )("plots/calculate_beta/betas.npy", encoding="latin1", allow_pickle=True).item()

    extra_data = task(
        np.load,
        start_text="Loading concentrations_and_einasto.npy",
        end_text="Loaded concentrations_and_einasto.npy",
        fail_text="Failed to load concentrations_and_einasto.npy",
        exit_on_fail=True
    )("plots/calculate_beta/concentrations_and_einasto.npy", encoding="latin1", allow_pickle=True).item()

    crossing_times = task(
        np.load,
        start_text="Loading crossing_times.npy",
        end_text="Loaded crossing_times.npy",
        fail_text="Failed to load crossing_times.npy",
        exit_on_fail=True
    )("plots/calculate_beta/crossing_times.npy", encoding="latin1", allow_pickle=True).item()

    # masses = np.array(table["M_halo @ z=0 [M_sun]"])
    peak_masses = np.array(table["M_peak [M_sun]"])
    stellar_masses = np.array(table["M_star [M_sun]"])
    peak_ages = table["min( t_peak, t(M=1e8M_sun)) [Gyr]"]
    betas = np.array(table["beta"])
    # stds = table["residual sigma"]

    # print(len(peak_masses), [ len(crossing_times[f"t_{n}"]["luminous"]) +  len(crossing_times[f"t_{n}"]["dark"]) for n in range(10,100,10) ])

    number_of_processors = 30
    with pymp.Parallel(number_of_processors, if_=True) as p: #set if_=False to disable parallelism
        for n in p.xrange(1,100,1):
        # for n in range(1,100,1):
            # Make one plot for every t_n comparing t_n for luminous halos vs dark halos: we expect to find that luminous halos have a higher median t_n

            cmap = cubehelix.cmap(reverse=True, start=(3/4.)*np.pi, rot=0.65, maxLight=0.8, minLight=0.2, minSat=1.3, maxSat=1.3)
            # cmap.set_bad(color="#CCCCCC")
            norm = mpl.colors.Normalize(vmin=0, vmax=1.5)

            f, ax = plt.subplots(1, 1, figsize=(6,4))

            colors = {
                "luminous": "#FF6666"
                ,"dark": "#0066FF"
            }

            for halo_type in ("luminous","dark"):

                tNs = np.array(crossing_times[f"t_{n}"])
                # tNs_halo_type = np.array(crossing_times[f"t_{n}"]["halo_type"])

                mass_bin_centers = []
                median_tNs = []
                median_tNs_25 = []
                median_tNs_75 = []

                number_of_halos_per_bin = 10
                # print(len(peak_masses))
                mass_bin_edges = peak_masses[::20] #np.logspace(7, 11, 20)
                # print(mass_bin_edges)
                for left_bin_edge, right_bin_edge in zip(mass_bin_edges[:-number_of_halos_per_bin],mass_bin_edges[number_of_halos_per_bin:]):
                    mass_bin_args = np.argwhere( np.logical_and(peak_masses > left_bin_edge, peak_masses < right_bin_edge) )

                    # if left_bin_edge > 2e8:
                    #     print(f"{len(mass_bin_args)} {halo_type} halos between ({np.round(np.log10(left_bin_edge),2)},{np.round(np.log10(right_bin_edge),2)}).")
                    # print(len(masses > left_bin_edge), len([ item for item in (masses > left_bin_edge) if item ]))
                    mass_bin_args = np.array([ int(mass_bin_args[i][0]) for i in range(len(mass_bin_args)) ]) #

                    if len(mass_bin_args) == 0:
                        continue

                    tNs_mass_bin = tNs[ mass_bin_args ]

                    mass_bin = peak_masses[ mass_bin_args ]
                    stellar_mass_bin = stellar_masses[ mass_bin_args ]

                    mass_bin_center = 10.**( (np.log10(left_bin_edge) + np.log10(right_bin_edge))/2.0 )

                    halo_args = np.argwhere( ( stellar_mass_bin > 0 if halo_type=="luminous" else stellar_mass_bin == 0 ) )
                    if len(halo_args) > 0:

                        halo_args = np.array([ int(halo_args[i][0]) for i in range(len(halo_args)) ])
                        mass_bin = mass_bin[ halo_args ]
                        tN_mass_bin = tNs_mass_bin[ halo_args ]

                        median_tNs.append( np.median(tN_mass_bin) )
                        median_tNs_25.append( np.percentile( tN_mass_bin, 25 ) )
                        median_tNs_75.append( np.percentile( tN_mass_bin, 75 ) )

                        mass_bin_centers.append( mass_bin_center )

                # Now plot the medians:
                # ax.fill_between(
                #     mass_bin_centers, # x values
                #     median_tNs_25, # y lower bound
                #     median_tNs_75, # y upper bound
                #     facecolor=colors[halo_type],
                #     linewidth=0,
                #     zorder=29,
                #     alpha=0.4
                # )
                #
                # ax.plot( mass_bin_centers, median_tNs, 'g', linewidth=2, color=colors[halo_type], zorder=31, alpha=1, label=r"Median $t_{"+str(n)+"}$, "+halo_type+" halos" )
                # ax.plot( mass_bin_centers, median_tNs, 'o', color=colors[halo_type], markersize=6, zorder=32 )

            # Plot the individual halos in the background:
            for peak_mass, stellar_mass, age, tN, beta in zip(peak_masses, stellar_masses, peak_ages, tNs, betas):
                marker, markersize, markeredgewidth, zorder, alpha, color = ('*',9,0.6,20,1,"#333333") if stellar_mass > 0 else ('o',4,0,15-age,0.6,"#BBBBBB")
                ax.plot( peak_mass, tN, marker, markersize=markersize, markeredgewidth=0, color=cmap(norm(beta)), zorder=zorder, alpha=0.6 )


            # ax.axhline(y=1,linestyle="--",linewidth=2,color="#999999",alpha=0.5,zorder=40)



            ax.set_xlabel(r"$M_{peak}$ [$M_\odot$]")
            ax.set_ylabel(r"$t_{"+str(n)+"}$ [$Gyr$]")
            ax.set_title(r"$M(t_{100}) = 10^8 M_{\odot}$")

            ax.set_xscale('log')
            # ax.set_yscale('log')

            ax.set_xlim([1e7,1e11])
            ax.set_ylim([0,14])

            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = plt.colorbar(sm, shrink=1.0, ax=ax, pad=0.04) #, ticks=np.linspace(0,1,5))
            # cbar = plt.colorbar(sm, location="top", shrink=1.0, ax=axs[1], pad=0.04)

            cbar.set_label(r'Beta', rotation=270)
            cbar.ax.get_yaxis().labelpad = 15
            # cbar.ax.invert_yaxis()

            ax.plot([],[],'o',markersize=4, markeredgewidth=0, color="#BBBBBB", label="Dark halos")
            ax.plot([],[],'*',markersize=9, markeredgewidth=0, color="#333333", label="Luminous halos")
            ax.legend(frameon=False, prop={'size': 8}, loc='upper right')

            plt.tight_layout()

            task(
                plt.savefig,
                start_text=f"Saving median_t{n}_vs_mpeak.png",
                end_text=f"Saved median_t{n}_vs_mpeak.png",
                fail_text=f"Failed to save median_t{n}_vs_mpeak.png",
                exit_on_fail=True
            )(f"plots/calculate_beta/tN_gif_plots/{n}.png")

            plt.close()


if make_gif:
    files = np.array(os.listdir("plots/calculate_beta/tN_gif_plots"))
    argsort_files = np.argsort([ "plots/calculate_beta/tN_gif_plots/"+( "0"*(6-len(file)) + file ) for file in files if file[-4:]==".png" ])

    files =  np.array([ "plots/calculate_beta/tN_gif_plots/"+file for file in files ])[ argsort_files ]

    total_duration = 4. #seconds
    t = total_duration/len(files)
    # times = np.linspace(  )
    dt = np.append( [t]*(len(files)-1), 0.5 )

    task(
        create_gif,
        start_text="Creating gif",
        end_text="Created gif",
        fail_text="Failed to create gif"
    )(files, list(dt), "plots/calculate_beta/tN_vs_Mpeak.gif")
