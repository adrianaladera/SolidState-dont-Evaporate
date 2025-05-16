from utils import ElectronicStructure   
import matplotlib.pyplot as plt
from collections import defaultdict
# from packages import Tools, Matplotlib, Pymatgen
from utils.packages import * 

METALS = ['Ag', 'Au', 'Hg', 'Cu']
CHALCS = ['Te', 'Se', "S"]

class ElectronicStructurePlotter:
    def __init__(self, path, savefig=False):
        if not os.path.exists(f"{path}/figures-and-data"):
            os.mkdir(f"{path}/figures-and-data")

        self.path = path
        self.savefig = savefig

    def color_bs_dos_plotter(
    self,
    bandz,
    dos,
    band_xlim,
    zero_to_efermi=True,
    # selected_indices=None,
    ylim=None,
    dos_xlim=None,
    smooth=False,
    vbm_cbm_marker=True,
    smooth_tol=0,
    smooth_k=3,
    smooth_np=100,
    bs_labels=None,
    marker_size=None
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
        path = self.path
        savefig = self.savefig
        fuck, (a0, a1) = plt.subplots(1, 2, figsize=(12,8), dpi=600, gridspec_kw={'width_ratios': [5, 3]})
        vbm_line, cbm_line = None, None

        ############ Plotting the DOS ############
        # colors = ["#DC267F", "#FFB000", "#FE6100", "#648FFF", "#785EF0", "#CC0000", "#FF7F50", "#FFD700", "#008000"] # red-green colorblind friendly!
        colors = ["#FF1493","#fc7b2b", "#FFB000",  "#8A2BE2", "#0269fa", "#e0a7fc"] # red-green colorblind friendly!
        # colors = ["#CC0000", "#FF7F50", "#FFD700", "#008000", "#2ACAEA", "#0000FF", "#8A2BE2", '#FF1493', "#666666", "#000000"]
        cunt = 4
        for orb in dos.get_dos_dict():
            x = dos.get_dos_dict()[orb]['densities']['1']
            norm_x = [(i - min(x)) / (max(x) - min(x)) for i in x]
            y =  dos.get_dos_dict()[orb]['energies']
            if str(orb) in METALS:
                a1.plot(norm_x, y, label=orb, c=colors[0], linewidth=4)
            elif str(orb) in CHALCS:   
                a1.plot(norm_x, y, label=orb, c=colors[1], linewidth=4)
            elif str(orb) == "C":
                a1.plot(norm_x, y, label=orb, c=colors[2], linewidth=4)
            elif str(orb) == "H":
                a1.plot(norm_x, y, label=orb, c=colors[3], linewidth=4)
            else: 
                a1.plot(norm_x, y, label=orb, c=colors[cunt], linewidth=4)
                cunt += 1
        a1.set_xlabel(r"Normalized DOS (E/eV)", fontsize=24)
        a1.set_xlim(-0.1, 1.0)
        xticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
        a1.set_xticks(xticks)
        if ylim is not None:
            a1.set_ylim(ylim[0], ylim[1])
            a1.set_ylabel(r"Energies (eV)", fontsize=24)
        else:
            a1.set_ylim(ylim[0], ylim[1])
        a1.set_xticklabels(xticks, fontsize=20)
        a1.set_yticklabels(np.arange(ylim[0], ylim[1]), fontsize=20)

        if not zero_to_efermi:
            a1.axhline(y = dos.get_dos_dict()['Ag']['efermi'], color = 'r', linestyle = '--', label='$E_{Fermi}$',linewidth=1.5)
            
        a1.grid()
        a1.legend(loc='upper right', fontsize=20)

        ################### Plotting the band structure ###############
        elements = ElectronicStructure.sort_elements(f"{path}/band/POSCAR")
        inorganics = elements[0]
        organics = elements[1]
        # group_dict = [{'elements':inorganics,'color':[255, 119, 36]},{'elements':organics,'color':[0,0,0]}]
        group_dict = [{'elements':inorganics,'color':[255, 119, 36]},{'elements':organics,'color':[0,0,0]}]
        print(elements)
        
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
                    print(f"cbm: {cbm}")
                    cb_point.append([cbm[0], cbm[1]])
                for vbm in data["vbm"]:
                    print(f"vbm: {vbm}")
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

    def sexier_brillouin_plot(self, bsp, band_lbls, lw, fs, ms, coords_are_cartesian=False):
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
        path = self.path
        savefig = self.savefig
        fig = plt.figure(figsize=(8,8), dpi=600)
        ax = fig.add_subplot(111, projection='3d')

        # Get all labels and lines
        all_labels = {}
        for k in bsp._bs[0].kpoints:
            if k.label:
                all_labels[k.label] = k.frac_coords

        all_lines = []
        for branch in bsp._bs[0].branches:
            kpts = bsp._bs[0].kpoints
            start_idx, end_idx = branch["start_index"], branch["end_index"]
            all_lines.append([kpts[start_idx].frac_coords, kpts[end_idx].frac_coords])

        # Selecting labels and corresponding plot coordinates
        print(all_labels)

        labels = {}
        for l in band_lbls:
            for a in all_labels:
                if '$' in l:
                    l2 = l.replace("$", "")
                    if '\\mid' in l2:
                        for kp in l2.split('\\mid'):
                            if str(kp) == str(a) and kp not in labels.keys():
                                # print(kp)
                                labels[kp] = all_labels[a]
                    elif str(l2) == str(a) and l2 not in labels.keys():
                        # print(l2)
                        labels[a] = all_labels[a]
                elif str(l) == str(a) and l not in labels.keys():
                    # print(l)
                    labels[a] = all_labels[a]

        print(labels)

        # Plot Brillouin zone lattice vectors in grey
        bz_lattice = bsp._bs[0].lattice_rec

        vertex1 = bz_lattice.get_cartesian_coords([0.0, 0.0, 0.0])
        directions = [
        bz_lattice.get_cartesian_coords([1.0, 0.0, 0.0]) - vertex1,  # X-direction
        bz_lattice.get_cartesian_coords([0.0, 1.0, 0.0]) - vertex1,  # Y-direction
        bz_lattice.get_cartesian_coords([0.0, 0.0, 1.0]) - vertex1   # Z-direction
        ]

        for direction in directions:
            ax.quiver(
                *vertex1, *direction,
                color='#7d7c7c', alpha=0.7, linewidth=1.5, arrow_length_ratio=0
            )
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

        # List of all line info
        deeznuts = []
        for line in all_lines:
            temp = []
            for pp in line:
                cunt = 0
                for key in labels.keys():
                    if key in all_labels:  # Ensure the key exists in labels
                        point = tuple(all_labels[key])
                        if np.array_equal(point, pp) and point not in temp:
                            temp.append(point)
                cunt += 1
            if len(temp) == 2 and tuple(temp) not in deeznuts:
                deeznuts.append(temp)

        # Plot the k-path in red
        lines = deeznuts
        if lines is not None:
            for line in lines:
                for k in range(1, len(line)):
                    vertex1 = line[k - 1]
                    vertex2 = line[k]
                    if not coords_are_cartesian:
                        if bz_lattice is None:
                            raise ValueError("coords_are_cartesian False requires the lattice")
                        vertex1 = bz_lattice.get_cartesian_coords(vertex1)
                        vertex2 = bz_lattice.get_cartesian_coords(vertex2)
                    ax.plot(*zip(vertex1, vertex2, strict=True), color="#FF0000", lw=lw)

        if labels is not None:
            for k, coords in labels.items():
                label = k
                if k.startswith("\\") or k.find("_") != -1:
                    label = f"${k}$"
                off = 0.02
                if coords_are_cartesian:
                    coords = np.array(coords)
                else:
                    if bz_lattice is None:
                        raise ValueError("coords_are_cartesian False requires the lattice")
                    coords = bz_lattice.get_cartesian_coords(coords)
                ax.text(*(coords + off), s=label, color='b', fontsize=fs)

            es.plot_points(
                labels.values(),
                bz_lattice,
                coords_are_cartesian=coords_are_cartesian,
                fold=False,
                ax=ax,
            )
        
        minz, maxz = [],[]    
        for l in labels:
            minz.append(min(labels[l]))
            maxz.append(max(labels[l]))

        ax.grid(False)
        ax.set_axis_off()
        
        ax.set_xlim(-0.7,0.7)
        ax.set_ylim(-0.7,0.7)
        ax.set_zlim(-0.7,0.7)
        ax.view_init(elev=25, azim=35)
        ax.axis("off")

        if not os.path.exists(f"{path}/figures-and-data"):
            os.mkdir(f"{path}/figures-and-data")
        if savefig:
            fig.savefig(f"{path}/figures-and-data/brillouin_zone.png", dpi=600)

    def sexy_dos_plot(self, dos, ylim=None, xlim=None):
        '''Get element DOS and plot according to whether
            the element is a metal, chalcogen, or organic.
            Axes are scaled automatically according to the range
            of energies and densities. '''
        path = self.path
        savefig = self.savefig
        zero_to_efermi = dos.zero_at_efermi 
        fuck, ax = plt.subplots(1, figsize=(6,4), dpi=600)

        colors = ["#FF1493","#fc7b2b", "#FFB000",  "#8A2BE2", "#0269fa", "#e0a7fc"] # red-green colorblind friendly!
        # colors = ["#DC267F", "#FFB000", "#FE6100", "#648FFF", "#785EF0", "#CC0000", "#FF7F50", "#FFD700", "#008000"] # first 5 are red-green colorblind friendly!
        # colors = ["#CC0000", "#FF7F50", "#FFD700", "#008000", "#2ACAEA", "#0000FF", "#8A2BE2", "#e0a7fc",'#FF1493', "#fca7dd", "#666666", "#000000"]
        cunt = 4
        ylims = []
        for orb in dos.get_dos_dict():
            y = dos.get_dos_dict()[orb]['densities']['1']
            x =  dos.get_dos_dict()[orb]['energies']
            ylims.append(max(y))
            if orb in METALS:
                ax.plot(x, y, label=orb, c=colors[0],linewidth=3)
            elif orb in CHALCS:
                ax.plot(x, y, label=orb, c=colors[1],linewidth=3)
            elif orb=="C":
                ax.plot(x, y, label=orb, c=colors[2],linewidth=3)
            elif orb=="H":
                ax.plot(x, y, label=orb, c=colors[3],linewidth=3)
            else:
                ax.plot(x, y, label=orb, c=colors[cunt],linewidth=3)
            # cunt += 1
        ax.set_ylabel(r"DOS (E/eV)")
        ax.set_xlim(xlim[0], xlim[1])
        ylim_max = max(ylims)
        if ylim is None:
            ax.set_ylim(-1, ylim_max+1)
        else:
            ax.set_ylim(ylim)
        ax.set_xlabel(r"Energies (eV)")

        if not zero_to_efermi:
                ax.axhline(y = dos.get_dos_dict()['Ag']['efermi'], color = 'r', linestyle = '--', label='$E_{Fermi}$',linewidth=1.5)
        else:
            ax.axhline(y = 0.0, color = '#000000', linestyle = '--', label='$E_{Fermi}$',linewidth=1.5)
            
        # ax.grid()
        ax.legend(loc='upper right')
    

        for spine in ax.spines.values():
            spine.set_linewidth(2)

        if savefig:
            fuck.savefig(f"{path}/element_dos.png", dpi=600)


    def sexy_orbital_plot(self, dos, ylim=None, species=None, xlim=None):
        '''Get element DOS and plot according to whether
            the element is a metal, chalcogen, or organic.
            Axes are scaled automatically according to the range
            of energies and densities. '''
        path = self.path
        savefig = self.savefig
        zero_to_efermi = dos.zero_at_efermi 

        fuck, ax = plt.subplots(1, figsize=(6,4), dpi=600)

        colors = [ "#8A2BE2",'#0269fa', "#FF1493"]
        markers = ["o", "^", "."]
        lw = [3,2,1]
        cunt = 0
        ylims = []
        for orb in dos.get_dos_dict():
            y = dos.get_dos_dict()[orb]['densities']['1']
            x =  dos.get_dos_dict()[orb]['energies']
            ylims.append(max(y))
            ax.plot(x, y, label=orb, c=colors[cunt], marker=markers[cunt], markersize=lw[cunt]+0.5,linewidth=lw[cunt])
            cunt += 1
        ax.set_ylabel(r"DOS (E/eV)")
        ax.set_xlim(xlim[0], xlim[1])
        ylim_max = max(ylims)
        if ylim is None:
            ax.set_ylim(-1, ylim_max+1)
        else:
            ax.set_ylim(ylim)
        ax.set_xlabel(r"Energies (eV)")
            
        # ax.grid()
        ax.legend(loc='upper right')

        for spine in ax.spines.values():
            spine.set_linewidth(2)

        # ax.legend(loc='upper right', fontsize=20)
        fuck.tight_layout()

        if savefig:
            fuck.savefig(f"{path}/figures-and-data/orbital_pdos_{species}.png", dpi=600)


    def color_bs_plotter(
        self,
        bandz,
        zero_to_efermi=True,
        xlim=None,
        ylim=None,
        smooth=False,
        vbm_cbm_marker=True,
        smooth_tol=0,
        smooth_k=3,
        smooth_np=100,
        bs_labels=None,
        marker_size=None
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
        path = self.path
        savefig = self.savefig
        fuck, ax = plt.subplots(1, figsize=(12,8), dpi=600)
        vbm_line, cbm_line = None, None

        elements = ElectronicStructure.sort_elements(f"{path}/band/POSCAR")
        inorganics = elements[0]
        organics = elements[1]
        group_dict = [{'elements':inorganics,'color':[255, 119, 36]},{'elements':organics,'color':[0,0,0]}]
        
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
        if xlim is None:
            x_min_range = xticks['distance'][0]
            x_max_range = xticks['distance'][-1]
        else:
            x_min_range = xticks['distance'][xlim[0]]
            x_max_range = xticks['distance'][xlim[1]]
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
        for ax in [ax]:
            for spine in ax.spines.values():
                spine.set_linewidth(3)

        fuck.tight_layout()

        if savefig:
            fuck.savefig(f"{path}/figures-and-data/colored_bandstructure.png", dpi=600)

        return labels_list, vbm_line, cbm_line



class BandAlignment:
    def __init__(self, path, matches, bs, savefig):
        self.path = path
        self.savefig = savefig
        self.matches = matches
        self.bs = bs

    def get_max_projections(self, proj_data, band_gap, is_vbm=True):
        """
        proj_data: list[list[dict[str, float]]]
            Outer list: bands
            Inner list: k-points
            Dict: {element: projection value}

        Returns:
            dict[element] = (max_value, kpoint_index, band_index)
        """
        max_proj = {}

        # Decrease from band_gap if is_vbm, else increase from band_gap
        band_range = range(band_gap, -1, -1) if is_vbm else range(band_gap, len(proj_data))

        for band in band_range:
            for k_idx, kpoint in enumerate(proj_data[band]):
                for element, value in kpoint.items():
                    if value > 0.1 and element not in max_proj:
                        max_proj[element] = (value, k_idx, band)
                    elif value > 0.05 and element not in max_proj:
                        max_proj[element] = (value, k_idx, band)

        return max_proj

    def get_band_proj_per_atom(self, matches, proj):
        '''
        Gets the indices of the VBM and CBM, then for both the 
        VBM and CBM, gets the projections for each atom onto the 
        band.

        Args:
            matches: (3D array) dims: n_kpoints x n_bands x len([eigenvalue, binary])
            proj: (dict) BSVasprun.get_band_structure().get_projection_on_elements()

        Returns:
            vbm_data: (dict) {(str) atomic species: (float) first corresponding projection
            > 0.1 decreasing from VBM value}
            cbm_data: (dict) {(str) atomic species: (float) first corresponding projection
            > 0.1 increasing from CBM value}

        '''
        vbm, cbm = 0, 0

        # Efficient band gap detection
        for kp in matches:
            for band in range(len(kp) - 1):
                e1, e2 = float(kp[band][1]), float(kp[band + 1][1])
                if e1 != 0 and e2 == 0:
                    vbm, cbm = band, band + 1
                    break
            if cbm > 0:  # break outer loop once found
                break

        print(f"VBM: {vbm}, CBM: {cbm}")

        vbm_data = self.get_max_projections(proj[Spin.up], band_gap=vbm, is_vbm=True)
        cbm_data = self.get_max_projections(proj[Spin.up], band_gap=cbm, is_vbm=False)

        print("VBM projections:", vbm_data)
        print("CBM projections:", cbm_data)

        # Remove elements we don't care about
        for unwanted in [0]: #, 'H']:
            vbm_data.pop(unwanted, None)
            cbm_data.pop(unwanted, None)

        for atom in vbm_data:
            if atom in cbm_data:
                print(f"{atom} — VBM: {vbm_data[atom]}, CBM: {cbm_data[atom]}")

        return vbm_data, cbm_data

    def get_band_alignment(self, vbm_data, cbm_data, bs):
        ''' 
        Gets the max projection in the VBM and min projection 
        in the CBM of all the atomic species, for both the 
        inorganic and organic species.

        Args:
            vbm_data: (dict) atomic species and their corresponding 
            projections from the VBM
            cbm_data: (dict) atomic species and their corresponding 
            projections from the CBM
            bs: (pymatgen.io.vasp.outputs.BSVasprun.get_band_structure())

        Returns: 
            tuple(
                vbm_final: (dict) {(str) "Inorganic": (float) max value of 
                projections, (str) "Organic": (float) max value of projections}
                cbm_final: (dict) {(str) "Inorganic": (float) min value of 
                projections, (str) "Organic": (float) min value of projections}
            )

        '''
        element_gaps = {}
        vbm_final, cbm_final = {}, {}

        vbm_inorg, cbm_inorg, vbm_org, cbm_org = [], [], [], []

        for atom in vbm_data:
            vbm_val = bs.bands[Spin.up][vbm_data[atom][2]][vbm_data[atom][1]]
            cbm_val = bs.bands[Spin.up][cbm_data[atom][2]][cbm_data[atom][1]]
            element_gaps[atom] = cbm_val - vbm_val

            if atom in METALS or atom in CHALCS:
                vbm_inorg.append(vbm_val)
                cbm_inorg.append(cbm_val)
            else:
                vbm_org.append(vbm_val)
                cbm_org.append(cbm_val)

        # Get frontier levels
        if vbm_inorg:
            vbm_final["Inorganic"] = max(vbm_inorg)
        if vbm_org:
            vbm_final["Organic"] = max(vbm_org)
        if cbm_inorg:
            cbm_final["Inorganic"] = min(cbm_inorg)
        if cbm_org:
            cbm_final["Organic"] = min(cbm_org)

        return vbm_final, cbm_final

    # def plot_band_alignment(self, vbm_gaps, cbm_gaps):
    def plot_band_alignment(self):
        ''' 
        Plots the band alignment for the inorganic and organic bands.

        Args:
            vbm_gaps: (dict) the VBM values for the inorganic and organic bands
            cbm_gaps: (dict) the CBM values for the inorganic and organic bands
        '''
        path = self.path
        savefig = self.savefig
        matches = self.matches
        bs = self.bs

        proj = bs.get_projection_on_elements()
        vbm_shit, cbm_shit = self.get_band_proj_per_atom(matches, proj)
        vbm_gaps, cbm_gaps = self.get_band_alignment(vbm_shit, cbm_shit, bs)

        atoms = list(vbm_gaps.keys())
        x_vals = list(range(len(atoms)))

        fig, ax = plt.subplots(figsize=(5,4), dpi=600)

        for i, atom in enumerate(atoms):
            vbm = vbm_gaps[atom]
            cbm = cbm_gaps[atom]
            bandgap = cbm - vbm

            # Plot VBM and CBM lines
            ax.hlines(y=vbm, xmin=i - 0.4, xmax=i + 0.4, colors='blue', label='VBM' if i == 0 else "")
            ax.hlines(y=cbm, xmin=i - 0.4, xmax=i + 0.4, colors='red', label='CBM' if i == 0 else "")
            
            # Annotate with band gap arrow
            ax.annotate(
                '', 
                xy=(i + 0.1, cbm), 
                xytext=(i + 0.1, vbm), 
                arrowprops=dict(arrowstyle='<->', color='black', lw=2)
            )
            ax.text(
                i + 0.2, 
                (vbm + cbm) / 2, 
                f'{bandgap:.3f} eV', 
                color='gray', 
                fontsize=10,
                verticalalignment='center'
            )

        ax.set_xticks(x_vals)
        ax.set_xticklabels(atoms, rotation=45)
        ax.set_ylabel("Energy (eV)")
        ax.set_title("Band Edges and Gaps")
        ax.legend()
        ax.grid(True, linestyle='--', alpha=0.5)
        plt.tight_layout()
        if savefig:
            plt.savefig(f"{path}/figures-and-data/band_alignment.png", dpi=600)
        plt.show()
