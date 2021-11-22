#!/usr/bin/env python
'''
DESCRIPTION

    Convert python notebook to script

FOR HELP

    python make_cuts_in_regions.py --help

AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
CREATED ON:  2020-12-15
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''


import sys, argparse, datetime
# Create count table:
from multiprocessing import Pool
import pysam
from glob import glob
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt
from collections import Counter, defaultdict
from itertools import product

import numpy as np
import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['figure.dpi'] = 300

from scipy.ndimage import gaussian_filter
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator
from singlecellmultiomics.features import FeatureContainer
from singlecellmultiomics.bamProcessing.bamBinCounts import get_binned_counts_prefixed

import os
from singlecellmultiomics.utils.pandas import createRowColorDataFrame
import singlecellmultiomics.features
import pyBigWig
from matplotlib import rcParams

import pickle

# %matplotlib inline
rcParams['font.family'] = 'Helvetica'
rcParams['font.sans-serif'] = ['Helvetica']

def hex_to_rgb(value, normalize=False):
    value = value.lstrip('#')
    lv = len(value)
    tup = tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))
    if normalize:
        tup = [i / 255. for i in tup]
    return tup

def create_gene_models(contig: str,
                       start :int,
                       end : int,
                       ax, 
                       features,
                       overlap_dist: int = 20_000,
                       gene_height:float = 0.002,
                       spacer:float = 0.035,
                       minlen:int=0, 
                       plot_gene_bodies:bool=True, 
                       plot_gene_names:bool=True, 
                       plot_tss_markers:bool=False, 
                       exon_height:float = 0.010):

    
    gene_y = {}
    ymax = 0

    # Coordinates for tss symbols:
    tss = {strand:{'x':[], 'y':[]} for strand in '+-'}

    
    for fs,fe,name,strand, feature_meta in features.findFeaturesBetween(contig, start, end):#features.features[contig]:


        if abs(fe-fs)<minlen:
            continue

        if not (((fs>=start or fe>=start) and (fs<=end or fe<=end))):
            continue


        feature_meta = dict(feature_meta)
        
        # Only plot genes: (type==gene)
        if feature_meta.get('type') == 'gene':

            if not 'gene_name' in  feature_meta or feature_meta.get('gene_name').startswith('AC'):
                continue

            # Determine g-y coordinate:
            #print(f"Determining Y-coord for {feature_meta.get('gene_name')}")
            #print(gene_y)
            gy_not_avail = set()
            for gene,(s,e,loc) in gene_y.items():


                if (s+overlap_dist>=fs and s-overlap_dist<=fe) or (e+overlap_dist>=fs and e-overlap_dist<=fe):
                    # Overlap:
                    #print(f'Overlapping with {loc}')
                    gy_not_avail.add(loc)

            gy = 0
            while gy in gy_not_avail:
                gy+=1


            gene_y[name] = (fs,fe,gy)
            y_offset = gy * spacer

            #print(f'Picked {gy}, {y_offset}, {-gene_height*0.5 + y_offset} ')

            ymax = max(y_offset+gene_height,ymax)

            r = Rectangle((fs,-gene_height*0.5 + y_offset), fe-fs, gene_height, angle=0.0, color='grey')

            # Store the direction and coordinates:
            tss[strand]['x'].append(fs if strand=='+' else fe)
            tss[strand]['y'].append( y_offset)

            if plot_gene_bodies:
                ax.add_patch( r )

            #The gene name might be outside the plot..
            if plot_gene_names:
                ax.text((fe+fs)*0.5,-1.6*exon_height + y_offset,
                        feature_meta.get('gene_name'),
                        horizontalalignment='center',
          verticalalignment='center',fontsize=4,clip_on=True)

    # Plot the tss markers:
    if plot_tss_markers:
        for strand, marker in (('+','>'),('-','<')):
            if len(tss[strand]['x'])>0:
                ax.scatter( tss[strand]['x'], tss[strand]['y'], marker=marker, s=10, zorder=3,c='k',edgecolor='white',linewidths=0.1)
            ax.set_xlim(start,end)

    # Format the axis:
    ax.set_ylim(-0.01,ymax+0.01)
    ax.set_yticks([])
    ax.set_xlabel(f'chr{contig} location bp', fontsize=6)

    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.ticklabel_format(useOffset=False, style='plain')
    ax.tick_params(axis='x', labelsize=4 )
    ax.tick_params(length=1)
                
def create_sc_region_plot(
                        bam_dict:dict,
                        normalize_to_counts:pd.DataFrame,
                          region:tuple,
                        selected_cells:list,
                        features:singlecellmultiomics.features.FeatureContainer=None,
                          outprefix:str='./plots/myplot',
                          bwpath:str =None,
                    region_bin_size:int = 250,
                    trace_sigma:float = 7,
                    sigma:float = 10,
                    sigma_cells:float = 0.00001,
                    sigma_bw:float = 4,
                    bw_log_mode:bool = False,
                    bw_ticks = None, # A list of ticks, an integer with the amount of ticks
                    PCT_COLOR:float = 97, 
                    second_bulk_axis_label:str='Second axis',
                    trace_normalizer:str  = 'mean', # or zscore
                    cell_normalizer:str = 'raw', # or norm
                    
                    # Meta data for row colors and grouping column name
                    meta_frame:pd.DataFrame=None,
                    bulk_track_meta_column_group:str = None,
                    color_cname:str = 'colorcodergb',
                    create_axis_grids=True,
                    bulk_lw:float = 1.5,
    
                    overlap_dist:int = 20_000,
                    gene_height:float = 0.002,
                    spacer:float = 0.035,

                    outsuffix = "pdf",
    
                    colorbar_label=None,
                    jfigsize:tuple = (5, 8),

                    norm_by = "mean",
                    skip_clustermap = False,

                    outtxtprefix = None,
                    
                    min_gene_len_to_plot:int = 0,
                    gene_models_height:float = 0.2,
                    plot_tss_markers:bool = False,
                    plot_gene_names:bool=True,
                    plot_gene_bodies:bool = True):

    sns.set(font_scale=0.5) 
    sns.reset_orig()
    mpl.rcParams['figure.dpi'] = 300

    contig,start,stop = region
    end = stop

    if bwpath is not None:

        p = pyBigWig.open(bwpath)
        bw_x = []
        bw_y = []
        for s,e,value in p.intervals('chr'+contig,start-10000,end+10000):
            bw_x.append( (s+e)*0.5 )
            bw_y.append(value)


    region = contig, start-1, stop-1

    # print(region)
    region_counts = get_binned_counts_prefixed(bam_dict, region_bin_size, regions=[ region ] )


    ## NORMALIZE
    region_counts = region_counts.fillna(0)
    if cell_normalizer=='raw':
        counts =  region_counts.T #
        if colorbar_label is None:
            caxlabel= 'counts'
    else:
        counts =(region_counts/normalize_to_counts.loc[region_counts.columns]).T
        if colorbar_label is None:
            caxlabel=f'spike-in normalised counts'
    
    if colorbar_label is not None:
        caxlabel = colorbar_label
    
    ##############

    contig, start, stop = region
    end = stop

    # Fill non intialized bins with zeros:

    # Add missing cells with zeroed-rows:
    add_cells = list( set(normalize_to_counts.index) .difference(set(counts.index)  ) )
    #print(add_cells)
    counts = counts.append([pd.Series(name=cell,dtype='float32') for cell in add_cells]).fillna(0)
    #counts.loc[add_cells,:] = 0

    add = []
    for i in np.arange(counts.columns[0][1], counts.columns[-1][1], region_bin_size):
        if not (contig,i) in counts.columns:
            add.append((contig,i))



    for a in add:
        counts[a] = 0

    counts = counts.sort_index(1)

    #selected_cells = normalize_to_counts.sum()[normalize_to_counts.sum()>5000].index
    selected_cells = [ cell for cell in selected_cells if cell in counts.index]

    order = counts.loc[ selected_cells ].index

    font = {'family' : 'Helvetica',
                    'weight' : 'normal',
                    'size'   : 7}

    mpl.rc('font', **font)
    
    qf = counts.loc[:, [(c,p) for c,p in counts if c==contig and p>=start and p<=end] ].sort_index()
    qf = qf.sort_index(1).sort_index(0)
    qf = qf.loc[order]

    # write output file matrix
    if outtxtprefix is not None:
        outmat = '.'.join([outtxtprefix, "mat.txt"])
        qf.to_csv(outmat)
    


    # Additive
    qf += pd.DataFrame(gaussian_filter(qf, sigma=(sigma_cells,sigma),mode='reflect'), index=qf.index, columns=qf.columns)
    # qf += pd.DataFrame(maximum_filter(qf, size=(sigma_cells,sigma),mode='reflect'), index=qf.index, columns=qf.columns)
    
    qf = qf.loc[meta_frame.index]
    # qf = qf.sort_index(1)
    # qf = qf.loc[order]

    counts = counts.sort_index()
    # counts = counts[meta_frame.index]

    # Create row colors
    if meta_frame is not None:
        row_colors = meta_frame[color_cname].to_frame()
        # lut = meta_frame[groupername]
        lut = {}
        lut[bulk_track_meta_column_group] = {meta_frame[bulk_track_meta_column_group][i]: tuple(meta_frame[color_cname][i]) for i in range(meta_frame.shape[0])}
        # print("Lut")
        # print(lut)
        # row_colors, lut = createRowColorDataFrame(meta_frame)
        # return(row_colors, meta_frame['colorcodergb'])
        # return(row_colors, lut)
    # return qf

    if not skip_clustermap:
        cm = sns.clustermap(qf,
               #z_score=0,
                row_cluster=False,
                col_cluster=False,
                vmax=np.percentile(qf,PCT_COLOR), #if cell_normalizer!='raw' else 1, #0.0005,
                vmin=0,
                dendrogram_ratio=0.1,
                row_colors=row_colors,
                figsize=jfigsize, cmap='Greys', cbar_kws={"shrink": .1},
                cbar_pos=(0.0, 0.5, 0.01, 0.16),)
    else:
        print("Skipping clustermap, but still intialize size")
        # cm = sns.clustermap(np.zeros(qf.shape), 
        #         row_cluster=False,
        #         col_cluster=False,
        #         vmax=np.percentile(qf,PCT_COLOR), #if cell_normalizer!='raw' else 1, #0.0005,
        #         vmin=0,
        #         dendrogram_ratio=0.1,
        #         # row_colors=row_colors,
        #         figsize=jfigsize, cmap='Greys', cbar_kws={"shrink": .1},
        #         cbar_pos=(0.0, 0.5, 0.01, 0.16),)

    # Add gene models
    if not skip_clustermap:
        fig = plt.gcf()
        ax = cm.ax_col_dendrogram
        ax.clear()
        ax.set_xticks([])
    else:
        fig = plt.Figure(figsize = jfigsize)
        ax = plt.subplot(111, aspect = 'auto')
        # ax = fig.add_axes(  (0, 0, jfigsize[0], jfigsize[1])  )
        # ax.clear()
        # ax.set_xticks([])
        # ax = None
    # Plot density:

    if trace_normalizer=='zscore':
        trace_norm_function = zscore
    elif trace_normalizer=='mean':
        trace_norm_function=lambda x: gaussian_filter([x],(trace_sigma,trace_sigma))[0]
    else:
        raise NotImplementedError()

    if meta_frame is None:
        bulk_track = qf.mean().droplevel(0)
        color= 'k'
        pd.Series(trace_norm_function(bulk_track), index=bulk_track.index).plot(color=color,lw=bulk_lw,ax=ax)
        
    
    grouper_name = bulk_track_meta_column_group
    
    # print("Before getting pseudobulk")
    if bulk_track_meta_column_group is not None:
        # print("Getting pseudobulk")
        grouper_name = bulk_track_meta_column_group #meta_frame.columns[0] # for example celltype
        
        for group in meta_frame[grouper_name].unique():
            '''
            if group != "Granulocytes":
                print("Skipping group" + group)
                continue
            '''
            color = lut[grouper_name][group]
            subset = qf[meta_frame[grouper_name]==group]
            meta_frame_subset = meta_frame[meta_frame[grouper_name]==group]
            if norm_by == "mean":
                print("Normalizing by mean")
                bulk_track = subset.mean().droplevel(0)
            else:
                print("Normalizing by colname:" + norm_by + ", group:" + group)
                norm_factor = meta_frame_subset[norm_by].sum()
                bulk_track = subset.sum().droplevel(0)
                print("Normalizing... with norm_factor")
                print(norm_factor)
                print("Normalizing...done")
                bulk_track = bulk_track / norm_factor
            # bulk_track = subset.mean().droplevel(0)
            # print("group: " + group)
            # print("color: ")
            # print(color)
            # print("subset: ")
            # print(subset.shape)
            # print(subset)
            pd.Series(trace_norm_function(bulk_track), index=bulk_track.index).plot(color=color,lw=bulk_lw,ax=ax)
            if outtxtprefix is not None:
                # write outpute bulk track
                outbulk = '.'.join([outtxtprefix, group, "bulktrack.txt"])
                bulk_track.to_csv(outbulk)

            
    
    # if not skip_clustermap:
    # Format the axis of the top bulk profiles:
    ax.tick_params(axis='y', labelcolor='black',labelsize=4)
    if trace_normalizer=='zscore':
        ax.set_ylabel('scChIC bulk\nz-score', color='grey', fontsize=4)  
        ax.set_yticks([-1,0,1,2])
    else:
        ax.set_ylabel('scChIC bulk\nmean counts', color='grey', fontsize=4)  
        #ax.set_yscale('log')
        if create_axis_grids:
            ax.grid(True, color='grey',which='minor')
            ax.grid(True, color='grey',which='major')
    ax.set_xlim(start,end)
    # else:
    #     print("Skip making axis of top bulk profiles because clustermap")
    

    # Plot bw track on secondary axis
    if bwpath is not None:
        ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
        color = 'tab:blue'
        ax2.set_ylabel(second_bulk_axis_label, color=color, fontsize=4)  # we already handled the x-label with ax1
        ax2.tick_params(axis='y', labelcolor=color,labelsize=4)
        ax2.plot(bw_x,gaussian_filter(np.array([bw_y]), sigma=(sigma_cells,sigma_bw),mode='reflect')[0], lw=bulk_lw,c='tab:blue')
        
        if bw_ticks is None:
            pass
        elif type(bw_ticks) is list:
            ax2.set_yticks([-2,-1,0,1,2])
        elif type(bw_ticks) is int:
            ax2.yaxis.get_major_locator().numtick = bw_ticks
        
        ax2.set_xlim(start,end)
        if create_axis_grids:
            ax2.grid(lw=0.3)
            ax2.set_axisbelow(True)
        if bw_log_mode:
            ax2.set_yscale('symlog')
        sns.despine(ax=ax2,top=True, right=False, left=False)

    sns.despine(ax=ax,top=True, right=False, left=False)

    if not skip_clustermap:
        cm.ax_heatmap.set_xticks([]) #np.arange(start,end, 1_000_000))
        cm.ax_heatmap.set_yticks([])
        cm.ax_heatmap.set_ylabel(f'{qf.shape[0]} single cells', fontsize=6)
        cm.ax_heatmap.tick_params(length=0.5)
        cm.ax_heatmap.set_xlabel(None)


        cm.cax.set_ylabel(caxlabel,fontsize=4)
        cm.cax.tick_params(labelsize=4)

        heatmap_start_x,heatmap_start_y, heatmap_end_x, heatmap_end_y = cm.ax_heatmap.get_position().bounds

        width = heatmap_end_x #-heatmap_start_x
        height = gene_models_height if features is not None else 0.05
        ax = fig.add_axes(  (heatmap_start_x, heatmap_start_y-height-0.02, width, height)  )
        ax.ticklabel_format(axis='x',style='sci')
        ax.set_yticks([])
        # despine the gene map
        sns.despine(ax=ax,left=True,right=True,top=True)
    else:
        print("Skipped adjusting heatmap because skip_clustermap")


    if not skip_clustermap: 
        if features is not None:
            create_gene_models(contig,start,end,features=features, ax=ax,minlen=min_gene_len_to_plot,plot_tss_markers=plot_tss_markers,
                           plot_gene_names=plot_gene_names,plot_gene_bodies=plot_gene_bodies,                    
        
                        overlap_dist = overlap_dist,
                        gene_height =gene_height,
                        spacer = spacer)
    else:
        print("Skipped making gene models because skip_clustermap")

    plt.savefig(f'{outprefix}_{region[0]}_{region[1]}_{region[2]}.{outsuffix}',bbox_inches='tight',format=outsuffix)
    
    # if meta_frame is not None:
    #     sns.palplot(lut[grouper_name].values())
    #     plt.xticks(range(len(lut[grouper_name])),lut[grouper_name].keys())
    #     plt.savefig(f'{outprefix}_{region[0]}_{region[1]}_{region[2]}_legend.png',bbox_inches='tight')
    
    return counts, qf



def main():
    parser = argparse.ArgumentParser(description='Convert python notebook to script')
    parser.add_argument('-infbam', metavar='PATH',
                        help='bam_dict_pickle')
    parser.add_argument('-inftotalcounts', metavar='PATH',
                        help='total_dict')
    parser.add_argument('-infnorm', metavar='PATH',
                        help='normalize_dict')
    parser.add_argument('-infgenelocs', metavar='PATH',
                        help='gene_dict')
    parser.add_argument('-inffeatures', metavar='PATH',
                        help='features_dict')
    parser.add_argument('-outdir', metavar='PATH',
                        help='Outdir. Will add gene, region and pdf automatically.')
    parser.add_argument('-radiusleft', metavar='basepairs', type=int, default=100000,
                        help='Add radius to start and end')
    parser.add_argument('-radiusright', metavar='basepairs', type=int, default=100000,
                        help='Add radius to start and end')
    parser.add_argument('-gene', metavar='gene name',
                        help='Gene name')
    parser.add_argument('-infmeta', metavar='infmeta',
                        help='infmeta')
    parser.add_argument('-mark', metavar='mark',
                        help='H3K4me1, H3K4me3 or H3K27me3')
    parser.add_argument('-sigma_cells', metavar = 'sigma', type=float, default=0.00001, help="Sigma smoothing across cells, maybe mark dependent")
    parser.add_argument('-sigma', metavar = 'sigma', type=float, default=6, help="Sigma smoothing across regions")
    parser.add_argument('-jsep', metavar = 'separator', default=",", help="Separator for meta dat, sometimes , sometimes \t")
    parser.add_argument('-jname', metavar = 'Name', default="_", help="Add name between mark and sigma, default _")
    parser.add_argument('-width', metavar = 'Width', default="5", type=int, help="Width probably inches")
    parser.add_argument('-height', metavar = 'Height', default="8", type=int, help="Height probably inches")
    parser.add_argument('-norm_by', metavar = 'mean or colname', default="mean", help="Colname from metadata for normalizing. Default is mean, uses bam mean")
    parser.add_argument('-colorcname', metavar = 'colname', default="colorcode", help="Color code usually colorcode or clustercol")
    parser.add_argument('-outfiletype', metavar = 'ImageType', default="pdf", help="Either pdf or png")
    parser.add_argument('-outtxtprefix', metavar = 'Output prefix', default=None, help="Path to write text files of input")
    parser.add_argument('-percentile', metavar='Percentile', type=float, default=99.5,
                        help='Perceentile. If too low background turns grey')
    parser.add_argument('-trace_sigma', metavar='trace_sigma', type=float, default=2,
                        help='Smoothing for bulk track')
    parser.add_argument('--rserver2hpc_prefix', action='store_true',
                        help='Swap prefix')
    parser.add_argument('--skip_clustermap', action='store_true',
                        help='Skip clustermap just show the bulk')
    parser.add_argument('--is_region', action='store_true',
                        help='Switch gene to region, skips using gene dict and assigns regions directly')
    parser.add_argument('-subdickey', default = "scchix",
                        help='Subdic used to access dictionary eg bam_dict_dict["K9m3"][subdickey]')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Suppress some print statements')
    parser.add_argument('--logfile', '-l', metavar='LOGFILE', default = None,
                        help='Write arguments to logfile')
    args = parser.parse_args()

    # store command line arguments for reproducibility
    CMD_INPUTS = ' '.join(['python'] + sys.argv)    # easy printing later
    # store argparse inputs for reproducibility / debugging purposes
    args_dic = vars(args)
    # ARG_INPUTS = ['%s=%s' % (key, val) for key, val in args_dic.iteritems()]  # for python2
    ARG_INPUTS = ['%s=%s' % (key, val) for key, val in args_dic.items()]  # for python3
    ARG_INPUTS = ' '.join(ARG_INPUTS)

    # Print arguments supplied by user
    if not args.quiet:
        if args.logfile is not None:
            sys.stdout = open(args.logfile, "w+")
        print(datetime.datetime.now().strftime('Code output on %c'))
        print('Command line inputs:')
        print(CMD_INPUTS)
        print ('Argparse variables:')
        print(ARG_INPUTS)

    if args.outfiletype == "pdf":
        print("Setting matploblib to pdf")
        matplotlib.use('pdf')
    elif args.outfiletype == "png":
        print("Setting matploblib to Agg")
        matplotlib.use('Agg')
    else:
        print(args.outfiletype, " outfiletype must be pdf or png" )

    percentile = args.percentile
    jmark = args.mark
    # assert jmark in ['H3K4me1', 'H3K4me3', 'H3K27me3', 'H3K9me3']
    print(jmark)

    print("Loading objects") 
    bam_dict_dict = pickle.load(open(args.infbam, "rb"))
    # swap prefixes
    if args.rserver2hpc_prefix:
        print("Replacing prefixes")
        infvec = []
        for b in bam_dict_dict[jmark][args.subdickey]:
            print(b)
            infvec.append(b.replace("/home/jyeung/hub_oudenaarden", "/hpc/hub_oudenaarden"))
        bam_dict_dict[jmark][args.subdickey] = infvec
        print("After:")
        print(infvec)
        
    total_count_per_cell_dict = pickle.load(open(args.inftotalcounts, "rb"))
    normalize_to_counts_dict = pickle.load(open(args.infnorm, "rb"))
    gene_locations = pickle.load(open(args.infgenelocs, "rb"))
    features = pickle.load(open(args.inffeatures, "rb"))
    print("Loading objects...done") 


    regions = []
    # extra_radius = 100000
    # extra_radius = args.radius
    radiusleft = args.radiusleft
    radiusright = args.radiusright
    genesvec = [args.gene]
    print(genesvec)

    # output_folder = '/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_jupyter/scchic_cuts_visualizations/' + jmark + "_with_pseudobulk_pdfs"
    output_folder = args.outdir

    for gene in genesvec :
        if not args.is_region:
            contig,start,end,strand =  gene_locations[gene]

            if strand == "+":
                jstart = start - radiusleft
                jend = start + radiusright
            elif strand == "-":
                jstart = end - radiusleft
                jend = end + radiusright
            else:
                print("strand either + or -, setting start and end:"  + strand)
                jstart = start - radiusleft
                jend = end + radiusright
                # regions.append( [contig,start-extra_radius, end+extra_radius] )
            assert jend - jstart > 0

        else:
            contig = gene.split(":")[0]
            jstartend = gene.split(":")[1]
            jstart = int(jstartend.split("-")[0])
            jend = int(jstartend.split("-")[1])
            strand = "+"

        # jname = gene
        regions.append( [contig, jstart, jend] )


    # Different thresholds for percentiles:
    # percentile = 99.5 # Higher values yield lower contrast


    # for jmark in jmarksvec:

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # infmeta = '/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/metadata_umap_celltype_cuts.' + jmark + '.txt'

    infmeta = args.infmeta
    jakeframe = pd.read_csv(infmeta, sep = args.jsep)
    jakeframe.index = pd.MultiIndex.from_tuples( (args.subdickey,cell) for cell in jakeframe['cell'] )

    jakeframe['colorcodergb'] = [hex_to_rgb(x, normalize=True) for x in jakeframe[args.colorcname]]

    selected_cells_meta = jakeframe['cell'].index

    jname = args.jname

    for region in regions:

        print(region)
        jstart = region[2]
        jend = region[1]
        # jname = args.gene

        print(jend)
        print(jstart)
        print(jend - jstart)
        # if (jend - jstart) > 400000:
        #     print("Region too big, skipping")
        #     next

        cell_subset = None #  [((m, plate), cell) for (m, plate), cell in total_count_per_cell.index if m==mark]

        bwpath = None
        jsigmacells = args.sigma_cells
        jsigma = args.sigma

        target_dir = f'{output_folder}'
        if not os.path.exists(target_dir):      
            os.makedirs(target_dir)
        counts, qf = create_sc_region_plot(
                            bam_dict_dict[jmark],
                            region=region,
                            selected_cells=selected_cells_meta,
                            # outprefix=f'{target_dir}/{jname}_{jmark}_{jsigmacells}_{jsigma}',
                            outprefix=f'{target_dir}/{jmark}_{jname}_{jsigmacells}_{jsigma}',
                            bwpath = bwpath,
                            region_bin_size = 500, # The size of the cells in the heatmap
                            normalize_to_counts=total_count_per_cell_dict[jmark], #.loc[cell_subset],
                            trace_sigma = args.trace_sigma, # Sigma for the bulk trace (gauss)
                            sigma = args.sigma, # Sigma for the single cell dots (gauss)
                            # sigma_cells = 0.00001, # Sigma for smoothing across cells. Very low values basically disable it.
                            sigma_cells = args.sigma_cells, # Sigma for smoothing across cells. Very low values basically disable it.
                            sigma_bw = 10, # Smoothing on the bigwig
                            bw_ticks=12, # Amount of yticks on the bigwig axis
                            PCT_COLOR = percentile, 
                            trace_normalizer  = 'mean', # or zscore

                            second_bulk_axis_label = None if bwpath is None else ('ChIP\nlog2 fold change' if 'log2' in bwpath  else 'ChIP coverage') ,
                            bw_log_mode = False,

                            meta_frame = jakeframe,
                            color_cname = 'colorcodergb',
                            bulk_track_meta_column_group = 'cluster',

                            outsuffix = args.outfiletype,
                            norm_by = args.norm_by,
                            skip_clustermap = args.skip_clustermap,
                            jfigsize = (args.width, args.height),

                            outtxtprefix = args.outtxtprefix,

                            features=features,
                            gene_height=0.0001,
                            gene_models_height=0.05,
                            create_axis_grids=False,
                            cell_normalizer = 'raw', # or norm
                            plot_tss_markers = True,
                            plot_gene_names=True,
                            min_gene_len_to_plot=1000,
                            plot_gene_bodies = True)

        plt.close('all')


if __name__ == '__main__':
    main()
