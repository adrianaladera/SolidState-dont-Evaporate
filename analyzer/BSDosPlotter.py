from utils import ElectronicStructure   
# from packages import Tools, Matplotlib, Pymatgen
from utils.packages import * 

def color_bs_dos_plotter(
    path,
    bandz,
    dos,
    band_xlim,
    zero_to_efermi=True,
    ylim=None,
    dos_xlim=None,
    smooth=False,
    vbm_cbm_marker=True,
    smooth_tol=0,
    smooth_k=3,
    smooth_np=100,
    bs_labels=None,
    marker_size=None,
    savefig=True
):
    """
    By Tess Smidt and Adriana Ladera. 
    Dos Plotter by Adriana Ladera.
    Get a matplotlib object for the bandstructures plot.
    Multiple bandstructure objs are plotted together if they have the
    same high symm path.
    Args:
        zero_to_efermi: Automatically subtract off the Fermi energy from
            the eigenvalues and plot (E-Ef).
        ylim: Specify the y-axis (energy) limits; by default None let
            the code choose. It is vbm-4 and cbm+4 if insulator
            efermi-10 and efermi+10 if metal
        smooth (bool or list(bools)): interpolates the bands by a spline cubic.
            A single bool values means to interpolate all the bandstructure objs.
            A list of bools allows to select the bandstructure obs to interpolate.
        smooth_tol (float) : tolerance for fitting spline to band data.
            Default is None such that no tolerance will be used.
        smooth_k (int): degree of splines 1<k<5
        smooth_np (int): number of interpolated points per each branch.
        bs_labels: labels for each band for the plot legend.
    """
    fuck, (a0, a1) = plt.subplots(1, 2, figsize=(12,8), dpi=600, gridspec_kw={'width_ratios': [5, 3]})
    vbm_line, cbm_line = None, None
    metals = ['Ag', 'Au', 'Hg', 'Cu']
    chalcs = ['Te', 'Se', "S"]

    ############ Plotting the DOS ############
    colors = ["#CC0000", "#FF7F50", "#FFD700", "#008000", "#2ACAEA", "#0000FF", "#8A2BE2", '#FF1493', "#666666", "#000000"]
    cunt = 4
    for orb in dos.get_dos_dict():
        x = dos.get_dos_dict()[orb]['densities']['1']
        y =  dos.get_dos_dict()[orb]['energies']
        if str(orb) in metals:
            a1.plot(x, y, label=orb, c=colors[0], linewidth=4)
        elif str(orb) in chalcs:   
            a1.plot(x, y, label=orb, c=colors[1], linewidth=4)
        elif str(orb) == "C":
            a1.plot(x, y, label=orb, c=colors[2], linewidth=4)
        elif str(orb) == "H":
            a1.plot(x, y, label=orb, c=colors[3], linewidth=4)
        else: 
            a1.plot(x, y, label=orb, c=colors[cunt], linewidth=4)
            cunt += 1
    a1.set_xlabel(r"DOS (E/eV)", fontsize=24)
    a1.set_xlim(dos_xlim[0], dos_xlim[1])
    if ylim is not None:
        a1.set_ylim(ylim[0], ylim[1])
        a1.set_ylabel(r"Energies (eV)", fontsize=24)
    else:
        a1.set_ylim(ylim[0], ylim[1])
    a1.set_xticklabels(np.arange(dos_xlim[0], dos_xlim[1]), fontsize=20)
    a1.set_yticklabels(np.arange(ylim[0], ylim[1]), fontsize=20)

    if not zero_to_efermi:
        a1.axhline(y = dos.get_dos_dict()['Ag']['efermi'], color = 'r', linestyle = '--', label='$E_{Fermi}$',linewidth=1.5)
        
    a1.grid()
    a1.legend(loc='upper right', fontsize=20)

    ################### Plotting the band structure ###############
    elements = ElectronicStructure.sort_elements(f"{path}/band/POSCAR")
    inorganics = elements[0]
    organics = elements[1]
    group_dict = [{'elements':inorganics,'color':[255, 70, 70]},{'elements':organics,'color':[0,0,0]}]
    
    if isinstance(smooth, bool):
        smooth = [smooth] * len(bandz._bs)

    handles = []
    vbm_min, cbm_max = [], []

    for ibs, bs in enumerate(bandz._bs):
        
        projections = bs.get_projection_on_elements()
        
        # Get projections into matrix form and order of elements in projections
        proj_array = ElectronicStructure.projections_to_array(projections) # [n_spin, n_bands, n_distances, n_elements]
        elem_keys = list(ElectronicStructure.get_element_projection_keys(projections)) # [n_elements]

        # Get groups and colors from group_dict
        color_matrix = np.array([np.array(g["color"]) for g in group_dict]) # [n_group, rgb]
        groups = [[elem_keys.index(elem) for elem in group['elements']] for group in group_dict]

        # Average color based on group occupation
        proj_totals = np.zeros(list(proj_array.shape[:-1]) + [len(groups)])
        for i, group in enumerate(groups):
            proj_totals[..., i] = proj_array[..., group].sum(axis=-1) / proj_array.sum(axis=-1)
        proj_colors = (proj_totals @ color_matrix) / 255 # [n_spin, n_bands, n_distances, 3]

        # set first bs in the list as ref for rescaling the distances of the other bands
        bs_ref = bandz._bs[0] if len(bandz._bs) > 1 and ibs > 0 else None

        if smooth[ibs]:
            # interpolation works good on short segments like branches
            data = bandz.bs_plot_data(zero_to_efermi, bs, bs_ref, split_branches=True)
        else:
            data = bandz.bs_plot_data(zero_to_efermi, bs, bs_ref, split_branches=False)
        
        # Rearrange colors based on shape of distances
        proj_colors_new = proj_colors.reshape(
            list(proj_colors.shape[:2]) 
            + list(np.array(data['distances']).shape) 
            + [proj_colors.shape[-1]])
        proj_colors_new_transpose = np.transpose(proj_colors_new, axes=[0, 2, 1, 3, 4])

        # remember if one bs is a metal for setting the ylim later
        one_is_metal = False
        if not one_is_metal and data["is_metal"]:
            one_is_metal = data["is_metal"]

        # remember all the cbm and vbm for setting the ylim later
        if not data["is_metal"] and data["vbm"] is not None and data["cbm"] is not None:
            bandgap = data["cbm"][0][1] - data["vbm"][0][1]
            vbm_line =  data["vbm"][0][1]
            cbm_line = data["cbm"][0][1]
            cbm_max.append(data["cbm"][0][1])
            vbm_min.append(data["vbm"][0][1])
        else:
            cbm_max.append(bs.efermi)
            vbm_min.append(bs.efermi)

        xticks = bandz.get_ticks()
        labels = xticks["label"]
        print("BAND LABELS:")
        
        band_labels = {}
        for i, l in enumerate(labels):
            band_labels[i] = l
        print(band_labels)

        for i, sp in enumerate(bs.bands):
            ls = "-" if str(sp) == "1" else "--"

            if bs_labels is None:
                bs_label = f"Band {ibs} {sp.name}"
            else:
                # assume bs_labels is Sequence[str]
                bs_label = f"{bs_labels[ibs]} {sp.name}"

            # handles.append(mlines.Line2D([], [], lw=2, ls=ls, color=colors[ibs], label=bs_label))

            distances, energies = data["distances"], data["energy"][str(sp)]
            colors = proj_colors_new_transpose[i]

            if smooth[ibs]:
                _, r = bandz._interpolate_bands(
                    distances,
                    colors[..., 0],
                    smooth_tol=smooth_tol,
                    smooth_k=smooth_k,
                    smooth_np=smooth_np,
                )

                _, g = bandz._interpolate_bands(
                    distances,
                    colors[..., 1],
                    smooth_tol=smooth_tol,
                    smooth_k=smooth_k,
                    smooth_np=smooth_np,
                )

                _, b = bandz._interpolate_bands(
                    distances,
                    colors[..., 2],
                    smooth_tol=smooth_tol,
                    smooth_k=smooth_k,
                    smooth_np=smooth_np,
                )

                distances, energies = bandz._interpolate_bands(
                    distances,
                    energies,
                    smooth_tol=smooth_tol,
                    smooth_k=smooth_k,
                    smooth_np=smooth_np,
                )
                # join all branches together
                distances = np.hstack(distances)
                energies = np.hstack(energies)
                colors = np.transpose(np.array([r, g, b]), [1, 2, 3, 0])
                # Interpolation can cause values to be outside valid rgb. Fix here.
                colors[colors < 0.] = 0.
                colors[colors > 1.] = 1.
            else:
                distances = np.array(distances).squeeze(0)
                energies = np.array(energies).squeeze(0)
                
            colors = np.transpose(colors, [1, 0, 2, 3])
            colors = colors.reshape(colors.shape[0], -1, 3)
            distances = np.repeat(distances.reshape(1, -1), energies.shape[0], axis=0)
            a0.scatter(distances.reshape(-1), energies.reshape(-1), c=colors.reshape(-1, 3), s=marker_size)

        vb_point, cb_point = [],[]
        # plot markers for vbm and cbm
        if vbm_cbm_marker and data["vbm"] is not None and data["cbm"] is not None:
            for cbm in data["cbm"]:
                cb_point.append([cbm[0], cbm[1]])
            for vbm in data["vbm"]:
                vb_point.append([vbm[0], vbm[1]])

        # Draw Fermi energy, only if not the zero
        if not zero_to_efermi:
            ef = bs.efermi
            # a0.axhline(ef, lw=2, ls="--", color=colors[ibs])
            a0.axhline(bs.efermi, lw=1.5, ls="--", color='r')
            a0.text(0.05, bs.efermi+0.1, f'E_f = {round(bs.efermi,3)} eV', color='r', fontsize=20)

    # defaults for ylim
    e_min = -4
    e_max = 4
    if one_is_metal:
        e_min = -10
        e_max = 10

    if ylim is None:
        if zero_to_efermi:
            if one_is_metal:
                # Plot A Metal
                a0.set_ylim(e_min, e_max)
            else:
                a0.set_ylim(e_min, max(cbm_max) + e_max)
        else:
            all_efermi = [b.efermi for b in bandz._bs]
            ll = min([min(vbm_min), min(all_efermi)])
            hh = max([max(cbm_max), max(all_efermi)])
            a0.set_ylim(ll + e_min, hh + e_max)
    else:
        a0.set_ylim(ylim)

    a0.set_xticks(xticks['distance'])
    a0.set_xticklabels(xticks['label'], fontsize=20)

    # setting the range to display the selected band labels
    if band_xlim is None:
        x_min_range = xticks['distance'][0]
        x_max_range = xticks['distance'][-1]
    else:
        x_min_range = xticks['distance'][band_xlim[0]]
        x_max_range = xticks['distance'][band_xlim[1]]
    a0.set_xlim(x_min_range,x_max_range)

    # returning labels to use for brillouin zone
    labels_list = []
    for i,j in zip(xticks["distance"], xticks["label"]):
        if x_min_range <= i <= x_max_range:
            labels_list.append(j)

    for tit in xticks['distance']:
        a0.axvline(x = tit, color = '#000000', linestyle = '-')
    a0.set_yticklabels(a0.get_yticks(), fontsize=20)

    # Main X and Y Labels
    a0.set_xlabel(r"$\mathrm{Wave\ Vector}$", fontsize=24)
    ylabel = r"$\mathrm{E\ -\ E_f\ (eV)}$" if zero_to_efermi else r"$\mathrm{Energy\ (eV)}$"
    a0.set_ylabel(ylabel, fontsize=24)

    # band gap information
    if len(cb_point) > 0 and len(vb_point) > 0:
        a0.axhline(y = cb_point[0][1], color = '#666666', linestyle = '-.', lw=3)
        a0.axhline(y = vb_point[0][1], color = '#666666', linestyle = '-.', lw=3)
        a0.annotate('', xy=(x_min_range+0.1,cb_point[0][1]), xytext=(x_min_range+0.1,vb_point[0][1]), arrowprops=dict(arrowstyle='<->', color='#000000', lw=3))
        a0.text(x_min_range+0.2, (cb_point[0][1]+vb_point[0][1])/2, f'{bandgap:.3f} eV', color='#666666', fontsize=20)

    for c in cb_point:
        a0.scatter(c[0], c[1], color="b", marker="o", s=100)
    for v in vb_point:
        a0.scatter(v[0], v[1], color="g", marker="o", s=100)

    # changing spine thickness
    for ax in [a0, a1]:
        for spine in ax.spines.values():
            spine.set_linewidth(3)

    fuck.tight_layout()

    if not os.path.exists(f"{path}/figures-and-data"):
        os.mkdir(f"{path}/figures-and-data")
    if savefig:
        fuck.savefig(f"{path}/figures-and-data/colored_bs_dos.png", dpi=600)

    return labels_list, vbm_line, cbm_line


def sexier_brillouin_plot(path, bsp, band_lbls, lw, fs, ms, savefig=True):
    '''Since the plot_brillouin() from Pymatgen produces
        the ugliest plot ever this is my version of plotting
        the Brillouin zone. It stores the k-path branches and
        labels from the bsp, and pulls the Brillouin zone from
        the primitive structure. Finally, the branches, labels,
        and Brillouin zone are plotted.
        
        path: the root containing the band/ directory
        bsp: a Pymatgen BSPlotter object
        band_lbls: the band range whose corresponding kpoint labels
        will be plotted, taken from xticks['labels'] in color_bs_dos_plotter()
        lw: linewidth
        fs: fontsize
        ms: marker size'''
    
    fig = plt.figure(figsize=(8,8), dpi=600)

    labels_str = []
    for l in band_lbls:
        if '$' in l:
            l2 = l.replace("$", "")
            if '\\mid' in l2:
                for kp in l2.split('\\mid'):
                    if kp not in labels_str:
                        labels_str.append(kp)
        else:
            labels_str.append(l)

    viewing_angles = [(-30,-120)]
    # titles = ['3D view', 'XY plane', 'XZ plane', 'YZ plane']
    axes = [fig.add_subplot(111, projection='3d')]

    # make labels and lines
    labels = {}
    for k in bsp._bs[0].kpoints:
        if k.label:
            labels[k.label] = k.frac_coords

    lines = []
    for branch in bsp._bs[0].branches:
        kpts = bsp._bs[0].kpoints
        start_idx, end_idx = branch["start_index"], branch["end_index"]
        lines.append([kpts[start_idx].frac_coords, kpts[end_idx].frac_coords])

    # matching up labels with the points that make up the line
    deeznuts = [] # list of all line info
    for line in lines:
        gottem = {} # dict of start and end points for line
        for key, value in labels.items():
            if np.array_equal(value, line[0]):
                gottem["start"] = [key, value] # start: [kpoint label, [array of 3 line coords]]
            if np.array_equal(value, line[1]): # end: [kpoint label, [array of 3 line coords]]
                gottem["end"] = [key, value]
        if "start" in gottem.keys() and "end" in gottem.keys():
            deeznuts.append(gottem)

    for ax, (elev, azim) in zip(axes, viewing_angles):

        # Plotting labels of path points in blue
        kpoints_list = []
        for l in labels:
            for kpoint in labels_str:
                if l == kpoint and l not in kpoints_list:
                    ax.scatter(labels[l][0], labels[l][1], labels[l][2], color='b', s=ms)
                    ax.text(labels[l][0], labels[l][1], labels[l][2], l, color='b', fontsize=fs)
                    kpoints_list.append(l) # only the labels indices entered will be plotted

        print(kpoints_list)
        print(labels_str)
        print(band_lbls)

        # Plotting branches in red
        for d in deeznuts:
            flag = 1
            for p in d.keys():
                if d[p][0] not in kpoints_list: # only the lines whose labels are in kpoints will be plotted
                    flag = 0
            if flag:
                x = [d["start"][1][0], d["end"][1][0]]
                y = [d["start"][1][1], d["end"][1][1]]
                z = [d["start"][1][2], d["end"][1][2]]
                ax.plot(x, y, z, color="#FF0000", lw=lw)
        
        structure = Structure.from_file(f"{path}/band/POSCAR")
        spg_analy = SpacegroupAnalyzer(structure) # should be the same as the prim_struct
        prim_struct = spg_analy.get_primitive_standard_structure(international_monoclinic=False)

        # Plotting the Brillouin zone in black
        for facet in prim_struct.lattice.get_brillouin_zone():
            x = [point[0] for point in facet]
            x.append(facet[0][0])
            y = [point[1] for point in facet]
            y.append(facet[0][1])
            z = [point[2] for point in facet]
            z.append(facet[0][2])
            ax.plot(x, y, z, color="#000000", lw=lw)
            
        minz, maxz = [],[]    
        for l in labels:
            minz.append(min(labels[l]))
            maxz.append(max(labels[l]))

        ax.grid(False)
        ax.set_axis_off()
        # ax.set_title(title, fontsize=16)
        ax.set_xlim(min(minz)*1.5, max(maxz)*1.5)
        ax.set_ylim(min(minz)*1.5, max(maxz)*1.5)
        ax.set_zlim(min(minz)*1.5, max(maxz)*1.5)
        ax.view_init(elev=elev, azim=azim)

    plt.tight_layout()
    plt.show()

    if not os.path.exists(f"{path}/figures-and-data"):
        os.mkdir(f"{path}/figures-and-data")
    if savefig:
        fig.savefig(f"{path}/figures-and-data/brillouin_zone.png", dpi=600)


def sexy_dos_plot(path, dos, ylim, xlim=None, savefig=True):
    '''Get element DOS and plot according to whether
        the element is a metal, chalcogen, or organic.
        Axes are scaled automatically according to the range
        of energies and densities. '''
    zero_to_efermi = dos.zero_at_efermi 

    fuck, ax = plt.subplots(1, figsize=(9,6), dpi=600)

    colors = ["#CC0000", "#FF7F50", "#FFD700", "#008000", "#2ACAEA", "#0000FF", "#8A2BE2", "#e0a7fc",'#FF1493', "#fca7dd", "#666666", "#000000"]
    cunt = 0
    for orb in dos.get_dos_dict():
        y = dos.get_dos_dict()[orb]['densities']['1']
        x =  dos.get_dos_dict()[orb]['energies']
        ax.plot(x, y, label=orb, c=colors[cunt], linewidth=3)
        cunt += 1
    ax.set_ylabel(r"DOS (E/eV)")
    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim)
    ax.set_xlabel(r"Energies (eV)")

    if not zero_to_efermi:
            ax.axhline(y = dos.get_dos_dict()['Ag']['efermi'], color = 'r', linestyle = '--', label='$E_{Fermi}$',linewidth=1.5)
    else:
        ax.axhline(y = 0.0, color = '#000000', linestyle = '--', label='$E_{Fermi}$',linewidth=1.5)
        
    ax.grid()
    ax.legend(loc='upper left')
 

    for spine in ax.spines.values():
        spine.set_linewidth(2)

    # ax.legend(loc='upper right', fontsize=20)
    fuck.tight_layout()
    if not os.path.exists(f"{path}/figures-and-data"):
        os.mkdir(f"{path}/figures-and-data")
    if savefig:
        fuck.savefig(f"{path}/figures-and-data/orbital_pdos.png", dpi=600)


def color_bs_plotter(
    path,
    bandz,
    zero_to_efermi=True,
    ylim=None,
    smooth=False,
    vbm_cbm_marker=True,
    smooth_tol=0,
    smooth_k=3,
    smooth_np=100,
    bs_labels=None,
    marker_size=None,
    savefig=True
):
    """
    By Tess Smidt and Adriana Ladera. 
    Dos Plotter by Adriana Ladera.
    Get a matplotlib object for the bandstructures plot.
    Multiple bandstructure objs are plotted together if they have the
    same high symm path.
    Args:
        zero_to_efermi: Automatically subtract off the Fermi energy from
            the eigenvalues and plot (E-Ef).
        ylim: Specify the y-axis (energy) limits; by default None let
            the code choose. It is vbm-4 and cbm+4 if insulator
            efermi-10 and efermi+10 if metal
        smooth (bool or list(bools)): interpolates the bands by a spline cubic.
            A single bool values means to interpolate all the bandstructure objs.
            A list of bools allows to select the bandstructure obs to interpolate.
        smooth_tol (float) : tolerance for fitting spline to band data.
            Default is None such that no tolerance will be used.
        smooth_k (int): degree of splines 1<k<5
        smooth_np (int): number of interpolated points per each branch.
        bs_labels: labels for each band for the plot legend.
    """
    fuck, ax = plt.subplots(1, figsize=(12,8), dpi=600)
    vbm_line, cbm_line = None, None

    elements = ElectronicStructure.sort_elements(f"{path}/band/POSCAR")
    inorganics = elements[0]
    organics = elements[1]
    group_dict = [{'elements':inorganics,'color':[255, 70, 70]},{'elements':organics,'color':[0,0,0]}]
    
    if isinstance(smooth, bool):
        smooth = [smooth] * len(bandz._bs)

    handles = []
    vbm_min, cbm_max = [], []

    for ibs, bs in enumerate(bandz._bs):
        
        projections = bs.get_projection_on_elements()
        
        # Get projections into matrix form and order of elements in projections
        proj_array = ElectronicStructure.projections_to_array(projections) # [n_spin, n_bands, n_distances, n_elements]
        elem_keys = list(ElectronicStructure.get_element_projection_keys(projections)) # [n_elements]

        # Get groups and colors from group_dict
        color_matrix = np.array([np.array(g["color"]) for g in group_dict]) # [n_group, rgb]
        groups = [[elem_keys.index(elem) for elem in group['elements']] for group in group_dict]

        # Average color based on group occupation
        proj_totals = np.zeros(list(proj_array.shape[:-1]) + [len(groups)])
        for i, group in enumerate(groups):
            proj_totals[..., i] = proj_array[..., group].sum(axis=-1) / proj_array.sum(axis=-1)
        proj_colors = (proj_totals @ color_matrix) / 255 # [n_spin, n_bands, n_distances, 3]

        # set first bs in the list as ref for rescaling the distances of the other bands
        bs_ref = bandz._bs[0] if len(bandz._bs) > 1 and ibs > 0 else None

        if smooth[ibs]:
            # interpolation works good on short segments like branches
            data = bandz.bs_plot_data(zero_to_efermi, bs, bs_ref, split_branches=True)
        else:
            data = bandz.bs_plot_data(zero_to_efermi, bs, bs_ref, split_branches=False)
        
        # Rearrange colors based on shape of distances
        proj_colors_new = proj_colors.reshape(
            list(proj_colors.shape[:2]) 
            + list(np.array(data['distances']).shape) 
            + [proj_colors.shape[-1]])
        proj_colors_new_transpose = np.transpose(proj_colors_new, axes=[0, 2, 1, 3, 4])

        # remember if one bs is a metal for setting the ylim later
        one_is_metal = False
        if not one_is_metal and data["is_metal"]:
            one_is_metal = data["is_metal"]

        # remember all the cbm and vbm for setting the ylim later
        if not data["is_metal"] and data["vbm"] is not None and data["cbm"] is not None:
            bandgap = data["cbm"][0][1] - data["vbm"][0][1]
            vbm_line =  data["vbm"][0][1]
            cbm_line = data["cbm"][0][1]
            cbm_max.append(data["cbm"][0][1])
            vbm_min.append(data["vbm"][0][1])
        else:
            cbm_max.append(bs.efermi)
            vbm_min.append(bs.efermi)

        xticks = bandz.get_ticks()
        labels = xticks["label"]
        print("BAND LABELS:")
        for i, l in enumerate(labels):
            print(i, l)

        for i, sp in enumerate(bs.bands):
            ls = "-" if str(sp) == "1" else "--"

            if bs_labels is None:
                bs_label = f"Band {ibs} {sp.name}"
            else:
                # assume bs_labels is Sequence[str]
                bs_label = f"{bs_labels[ibs]} {sp.name}"

            # handles.append(mlines.Line2D([], [], lw=2, ls=ls, color=colors[ibs], label=bs_label))

            distances, energies = data["distances"], data["energy"][str(sp)]
            colors = proj_colors_new_transpose[i]

            if smooth[ibs]:
                _, r = bandz._interpolate_bands(
                    distances,
                    colors[..., 0],
                    smooth_tol=smooth_tol,
                    smooth_k=smooth_k,
                    smooth_np=smooth_np,
                )

                _, g = bandz._interpolate_bands(
                    distances,
                    colors[..., 1],
                    smooth_tol=smooth_tol,
                    smooth_k=smooth_k,
                    smooth_np=smooth_np,
                )

                _, b = bandz._interpolate_bands(
                    distances,
                    colors[..., 2],
                    smooth_tol=smooth_tol,
                    smooth_k=smooth_k,
                    smooth_np=smooth_np,
                )

                distances, energies = bandz._interpolate_bands(
                    distances,
                    energies,
                    smooth_tol=smooth_tol,
                    smooth_k=smooth_k,
                    smooth_np=smooth_np,
                )
                # join all branches together
                distances = np.hstack(distances)
                energies = np.hstack(energies)
                colors = np.transpose(np.array([r, g, b]), [1, 2, 3, 0])
                # Interpolation can cause values to be outside valid rgb. Fix here.
                colors[colors < 0.] = 0.
                colors[colors > 1.] = 1.
            else:
                distances = np.array(distances).squeeze(0)
                energies = np.array(energies).squeeze(0)
                
            colors = np.transpose(colors, [1, 0, 2, 3])
            colors = colors.reshape(colors.shape[0], -1, 3)
            distances = np.repeat(distances.reshape(1, -1), energies.shape[0], axis=0)
            ax.scatter(distances.reshape(-1), energies.reshape(-1), c=colors.reshape(-1, 3), s=marker_size)

        vb_point, cb_point = [],[]
        # plot markers for vbm and cbm
        if vbm_cbm_marker and data["vbm"] is not None and data["cbm"] is not None:
            for cbm in data["cbm"]:
                cb_point.append([cbm[0], cbm[1]])
            for vbm in data["vbm"]:
                vb_point.append([vbm[0], vbm[1]])

        # Draw Fermi energy, only if not the zero
        if not zero_to_efermi:
            ef = bs.efermi
            # ax.axhline(ef, lw=2, ls="--", color=colors[ibs])
            ax.axhline(bs.efermi, lw=1.5, ls="--", color='r')
            ax.text(0.05, bs.efermi+0.1, f'E_f = {round(bs.efermi,3)} eV', color='r', fontsize=20)

    # defaults for ylim
    e_min = -4
    e_max = 4
    if one_is_metal:
        e_min = -10
        e_max = 10

    if ylim is None:
        if zero_to_efermi:
            if one_is_metal:
                # Plot A Metal
                ax.set_ylim(e_min, e_max)
            else:
                ax.set_ylim(e_min, max(cbm_max) + e_max)
        else:
            all_efermi = [b.efermi for b in bandz._bs]
            ll = min([min(vbm_min), min(all_efermi)])
            hh = max([max(cbm_max), max(all_efermi)])
            ax.set_ylim(ll + e_min, hh + e_max)
    else:
        ax.set_ylim(ylim)

    ax.set_xticks(xticks['distance'])
    ax.set_xticklabels(xticks['label'], fontsize=20)

    # setting the range to display the selected band labels
    x_min_range = xticks['distance'][0]
    x_max_range = xticks['distance'][-1]
    ax.set_xlim(x_min_range,x_max_range)

    # returning labels to use for brillouin zone
    labels_list = []
    for i,j in zip(xticks["distance"], xticks["label"]):
        if x_min_range <= i <= x_max_range:
            labels_list.append(j)

    for tit in xticks['distance']:
        ax.axvline(x = tit, color = '#000000', linestyle = '-')
    ax.set_yticklabels(ax.get_yticks(), fontsize=20)

    # Main X and Y Labels
    ax.set_xlabel(r"$\mathrm{Wave\ Vector}$", fontsize=24)
    ylabel = r"$\mathrm{E\ -\ E_f\ (eV)}$" if zero_to_efermi else r"$\mathrm{Energy\ (eV)}$"
    ax.set_ylabel(ylabel, fontsize=24)

    # band gap information
    if len(cb_point) > 0 and len(vb_point) > 0:
        ax.axhline(y = cb_point[0][1], color = '#666666', linestyle = '-.', lw=3)
        ax.axhline(y = vb_point[0][1], color = '#666666', linestyle = '-.', lw=3)
        ax.annotate('', xy=(x_min_range+0.1,cb_point[0][1]), xytext=(x_min_range+0.1,vb_point[0][1]), arrowprops=dict(arrowstyle='<->', color='#000000', lw=3))
        ax.text(x_min_range+0.2, (cb_point[0][1]+vb_point[0][1])/2, f'{bandgap:.3f} eV', color='#666666', fontsize=20)

    for c in cb_point:
        ax.scatter(c[0], c[1], color="b", marker="o", s=100)
    for v in vb_point:
        ax.scatter(v[0], v[1], color="g", marker="o", s=100)

    # changing spine thickness
    for spine in ax.spines.values():
        spine.set_linewidth(3)

    fuck.tight_layout()

    if not os.path.exists(f"{path}/figures-and-data"):
        os.mkdir(f"{path}/figures-and-data")
    if savefig:
        fuck.savefig(f"{path}/figures-and-data/colored_bandstructure.png", dpi=600)

    return vbm_line, cbm_line
