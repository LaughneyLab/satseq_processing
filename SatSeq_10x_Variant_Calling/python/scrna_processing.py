import warnings

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
from sklearn.linear_model import LinearRegression


## from plot.xtick_vertical
def plot_xtick_vertical(ax=None):
    """set xticklabels on ax to vertical instead of the horizontal default orientation"""
    if ax is None:
        ax = plt.gca()
    ax.tick_params(axis='x', labelrotation=90)
    if False: # old code
        xt = ax.get_xticks()
        if np.all(xt.astype(int) == xt):  # ax.get_xticks() returns floats
            xt = xt.astype(int)
        ax.set_xticklabels(xt, rotation='vertical')

## copied from filter.low_count
def filter_low_count(molecules, is_invalid, plot=False, ax=None, adjust_inflection_pt=0.9, num_rolling=10):
    """
    updates is_invalid to reflect cells whose molecule counts are below the inflection
    point of an ecdf constructed from cell molecule counts. Typically this reflects cells
    whose molecule counts are approximately <= 100.

    :param molecules: scipy.stats.coo_matrix, molecule count matrix
    :param is_invalid:  np.ndarray(dtype=bool), declares valid and invalid cells
    :param bool plot: if True, plot a summary of the filter
    :param ax: Must be passed if plot is True. Indicates the axis on which to plot the
      summary.
    :return: is_invalid, np.ndarray(dtype=bool), updated valid and invalid cellsg
    """

    # copy, sort, and normalize molecule sums
    ms = np.ravel(molecules.tocsr()[~is_invalid, :].sum(axis=1))
    idx = np.argsort(ms)[::-1]  # largest cells first
    norm_ms = ms[idx] / ms[idx].sum()  # sorted, normalized array

    # identify inflection point from second derivative
    cms = np.cumsum(norm_ms)
    d1 = np.diff(pd.Series(cms).rolling(num_rolling).mean()[num_rolling:])
    d2 = np.diff(pd.Series(d1).rolling(num_rolling).mean()[num_rolling:])
    try:
        # throw out an extra 5% of cells from where the inflection point is found.
        # these cells are empirically determined to have "transition" library sizes
        # that confound downstream analysis
        inflection_pt = np.min(np.where(np.abs(d2) == 0)[0])
        inflection_pt = int(inflection_pt * adjust_inflection_pt)
    except ValueError as e:
        if e.args[0] == ('zero-size array to reduction operation minimum which has no '
                         'identity'):
            print('Low count filter passed-through; too few cells to estimate '
                  'inflection point.')
            return is_invalid  # can't estimate validity
        else:
            raise

    vcrit = ms[idx][inflection_pt]

    is_invalid = is_invalid.copy()
    is_invalid[ms < vcrit] = True

    if plot and ax:
        cms /= np.max(cms)  # normalize to one
        ax.plot(np.arange(len(cms))[:inflection_pt], cms[:inflection_pt])
        ax.plot(np.arange(len(cms))[inflection_pt:], cms[inflection_pt:], c='indianred')
        ax.hlines(cms[inflection_pt], *ax.get_xlim(), linestyle='--')
        ax.vlines(inflection_pt, *ax.get_ylim(), linestyle='--')
        #ax.set_xticklabels([])
        ax.set_xlabel('putative cell')
        ax.set_ylabel('ECDF (Cell Size)')
        ax.set_title('Cell Size')
        ax.set_ylim((0, 1))
        ax.set_xlim((0, inflection_pt*3))

    return is_invalid


## from seqc.filter.high_mitochondrial_rna; with some modifications:
## - is it takes an argument "mt_prefix" which was hard-coded as 'MT-' in seqc but should be 'mt-' for the CIN samples
## - remove argument mini_summary_d which returns some summary statistics
def filter_high_mitochondrial_rna(molecules, gene_ids, is_invalid, mt_prefix='mt-', max_mt_content=0.2, 
                                  plot=False, ax=None, filter_on=True):
    """
    Sets any cell with a fraction of mitochondrial mRNA greater than max_mt_content to
    invalid.

    :param molecules: scipy.stats.coo_matrix, molecule count matrix
    :param gene_ids: np.ndarray(dtype=str) containing string gene identifiers
    :param is_invalid:  np.ndarray(dtype=bool), declares valid and invalid cells
    :param mt_prefix: string that identifies mitochondrial genes
    :param max_mt_content: float, maximum percentage of reads that can come from
      mitochondria in a valid cell
    :param bool plot: if True, plot a summary of the filter
    :param ax: Must be passed if plot is True. Indicates thenumcores axis on which to plot the
      summary.
    :return: is_invalid, np.ndarray(dtype=bool), updated valid and invalid cells
    """
    # identify % genes that are mitochondrial
    mt_genes = np.fromiter(map(lambda x: x.startswith(mt_prefix), gene_ids), dtype=bool)
    mt_molecules = np.ravel(molecules.tocsr()[~is_invalid, :].tocsc()[:, mt_genes].sum(
        axis=1))
    ms = np.ravel(molecules.tocsr()[~is_invalid, :].sum(axis=1))
    ratios = mt_molecules / ms

    if filter_on:
        failing = ratios > max_mt_content
        is_invalid = is_invalid.copy()
        is_invalid[np.where(~is_invalid)[0][failing]] = True
    else:
        is_invalid = is_invalid.copy()

    if plot and ax:
        if ms.shape[0] and ratios.shape[0]:
            plot_scatter_continuous(ms, ratios, colorbar=False, ax=ax, s=3)
        else:
            return is_invalid  # nothing else to do here
        if filter_on and (np.sum(failing) != 0):
            ax.scatter(ms[failing], ratios[failing], c='indianred', s=3)  # failing cells
        xmax = np.max(ms)
        ymax = np.max(ratios)
        ax.set_xlim((0, xmax))
        ax.set_ylim((0, ymax))
        ax.hlines(max_mt_content, *ax.get_xlim(), linestyle='--', colors='indianred')
        ax.set_xlabel('total molecules')
        ax.set_ylabel('mtRNA fraction')
        if filter_on:
            ax.set_title(
                'mtRNA Fraction: {:.2}%'.format(np.sum(failing) / len(failing) * 100))
        else:
            ax.set_title('mtRNA Fraction')
        plot_xtick_vertical(ax=ax)

    return is_invalid


## from plot.scatter.continuous
def plot_scatter_continuous(x, y, c=None, ax=None, colorbar=True, randomize=True,
                   remove_ticks=False, **kwargs):
        """
        wrapper for scatter wherein the coordinates x and y are colored according to a
        continuous vector c
        :param x, y: np.ndarray, coordinate data
        :param c: np.ndarray, continuous vector by which to color data points
        :param remove_ticks: remove axis ticks and labels
        :param args: additional args for scatter
        :param kwargs: additional kwargs for scatter
        :return: ax
        """

        if ax is None:
            ax = plt.gca()

        if c is None:  # plot density if no color vector is provided
            x, y, c = plot_density_2d(x, y)

        if randomize:
            ind = np.random.permutation(len(x))
        else:
            ind = np.argsort(c)

        sm = ax.scatter(x[ind], y[ind], c=c[ind], **kwargs)
        if remove_ticks:
            ax.xaxis.set_major_locator(plt.NullLocator())
            ax.yaxis.set_major_locator(plt.NullLocator())
        if colorbar:
            cb = plt.colorbar(sm)
            #cb.ax.xaxis.set_major_locator(plt.NullLocator())
            #cb.ax.yaxis.set_major_locator(plt.NullLocator())
        return ax
    
## from plot.scatter.density_2d
def plot_density_2d(x, y):
        """return x and y and their density z, sorted by their density (smallest to largest)

        :param x, y: np.ndarray: coordinate data
        :return: sorted x, y, and density
        """
        xy = np.vstack([np.ravel(x), np.ravel(y)])
        z = gaussian_kde(xy)(xy)
        return np.ravel(x), np.ravel(y), np.arcsinh(z)


#Here is a simpler alternative version of filter_low_coverage
#The seqc version seemed to choose a fairly arbitrary cut-off. Instead, just get mean/sd of coverage, 
#and filter cells with coverage < mean - cutoff*sd
def filter_low_coverage(molecules, reads, is_invalid, cutoff=3, plot=False, ax=None, filter_on=True):
    """
    Fits a two-component gaussian mixture model to the data. If a component is found
    to fit a low-coverage fraction of the data, this fraction is set as invalid. Not
    all datasets contain this fraction.

    For best results, should be run after filter.low_count()

    :param molecules: scipy.stats.coo_matrix, molecule count matrix
    :param reads: scipy.stats.coo_matrix, read count matrix
    :param is_invalid:  np.ndarray(dtype=bool), declares valid and invalid cells
    :param bool plot: if True, plot a summary of the filter
    :param ax: Must be passed if plot is True. Indicates the axis on which to plot the
      summary.
    :param filter_on: indicate whether low coverage filter is on
    :return: is_invalid, np.ndarray(dtype=bool), updated valid and invalid cells
    """
    ms = np.ravel(molecules.tocsr()[~is_invalid, :].sum(axis=1))
    rs = np.ravel(reads.tocsr()[~is_invalid, :].sum(axis=1))

    if ms.shape[0] < 10 or rs.shape[0] < 10:
        log.notify(
            'Low coverage filter passed-through; too few cells to calculate '
            'mixture model.')
        return is_invalid

    # get read / cell ratio, filter out low coverage cells
    ratio = rs / ms
    
    mean_ratio = np.mean(ratio)
    sd_ratio = np.std(ratio)
    
    if filter_on:
        failing = np.where(ratio < mean_ratio - sd_ratio * cutoff)[0]

        # set smaller mean as invalid
        is_invalid = is_invalid.copy()
        is_invalid[np.where(~is_invalid)[0][failing]] = True
        
    if plot and ax:
        logms = np.log10(ms)
        plot_scatter_continuous(logms, ratio, colorbar=False, ax=ax, s=3)
        ax.set_xlabel('log10(molecules)')
        ax.set_ylabel('reads / molecule')
        if filter_on:
            ax.set_title('Coverage: {:.2}%'.format(np.sum(failing) / len(is_invalid) * 100))
        else:
            ax.set_title('Coverage')
        xmin, xmax = np.min(logms), np.max(logms)
        ymin, ymax = np.min(ratio), np.max(ratio)
        print(f'ymax={ymax}')
        ax.set_xlim((xmin, xmax))
        ax.set_ylim((ymin, ymax))
        plot_xtick_vertical(ax=ax)

        # plot 1d conditional densities of two-component model
        # todo figure out how to do this!!

        # plot the discarded cells in red, like other filters
        if filter_on:
            ax.scatter(
                logms[failing], ratio[failing],
                s=4, c='indianred')
    return is_invalid


## from filter.low_gene_abundance
## residual_cutoff was hard-coded as 0.15 but seems too stringent
def filter_low_gene_abundance(molecules, is_invalid, residual_cutoff = 0.15, plot=False, ax=None, filter_on=True):
    """
    Fits a linear model to the relationship between number of genes detected and number
    of molecules detected. Cells with a lower than expected number of detected genes
    are set as invalid.

    :param molecules: scipy.stats.coo_matrix, molecule count matrix
    :param is_invalid:  np.ndarray(dtype=bool), declares valid and invalid cells
    :param bool plot: if True, plot a summary of the filter
    :param ax: Must be passed if plot is True. Indicates the axis on which to plot the
      summary.
    :return: is_invalid, np.ndarray(dtype=bool), updated valid and invalid cells
    """

    ms = np.ravel(molecules.tocsr()[~is_invalid, :].sum(axis=1))
    genes = np.ravel(molecules.tocsr()[~is_invalid, :].getnnz(axis=1))
    x = np.log10(ms)[:, np.newaxis]
    y = np.log10(genes)

    if not (x.shape[0] or y.shape[0]):
        return is_invalid

    # get line of best fit
    with warnings.catch_warnings():  # ignore scipy LinAlg warning about LAPACK bug.
        warnings.simplefilter('ignore')
        regr = LinearRegression()
        regr.fit(x, y)

    # mark large residuals as failing
    yhat = regr.predict(x)
    residuals = yhat - y
    failing = residuals > residual_cutoff

    is_invalid = is_invalid.copy()
    if filter_on:
        is_invalid[np.where(~is_invalid)[0][failing]] = True

    if plot and ax:
        m, b = regr.coef_, regr.intercept_
        plot_scatter_continuous(x, y, ax=ax, colorbar=False, s=3)
        xmin, xmax = np.min(x), np.max(x)
        ymin, ymax = np.min(y), np.max(y)
        lx = np.linspace(xmin, xmax, 200)
        ly = m * lx + b
        ax.plot(lx, np.ravel(ly), linestyle='--', c='indianred')
        if filter_on:
            ax.scatter(x[failing], y[failing], c='indianred', s=3)
        ax.set_ylim((ymin, ymax))
        ax.set_xlim((xmin, xmax))
        ax.set_xlabel('molecules (cell)')
        ax.set_ylabel('genes (cell)')
        if filter_on:
            ax.set_title('Low Complexity: {:.2}%'.format(np.sum(failing) / len(failing) * 100))
        else:
            ax.set_title('Low Complexity')
        plot_xtick_vertical(ax=ax)

    return is_invalid


# Function to optimize the number of PCs to retain based on the linear deviation method
# A linear deviation threshold of 1e-5 is the default, as this has empirically worked well for most cases
# Returns the optimal number of PCs to keep (sensitive to starting number of PCs)
# The optimize_dev_thr function can be used to test the threshold used here.

def optimize_pcs(adata, pca_name='pca', dev_thr=1e-5, save=True, fn='PCA_linear_deviation_optimization.png', output_dir=None):

    plt.style.use('default')
    
    if pca_name not in adata.uns.keys():
        exit(f'ERROR: {pca_name} not in adata.')
    
    x = range(len(adata.uns[pca_name]['variance_ratio']))
    y = np.cumsum(adata.uns[pca_name]['variance_ratio'])
    
    deviation = []
    slopes = []
    
    #dev_thr = 1e-5
    thr_hit=False
    
    opt_n_pcs_linear_dev=0
    
    fig = plt.subplots(figsize=(6,6))
    
    for i in range(len(x)-2):
        
        tmp_slope = (y[i+1] - y[i]) / (x[i+1] - x[i])
        tmp_int = y[i] - tmp_slope*x[i]
        
        plt.scatter(x, y, s=10, color='blue')
        
        y_line = [tmp_slope*x + tmp_int for x in x]
        
        tmp_dev = abs(y_line[i+2] - y[i+2])
        
        deviation.append(tmp_dev)
        
        slopes.append(tmp_slope)
        
        plt.plot(x, y_line, color='grey', linewidth=0.5)
        
        if (tmp_dev < dev_thr) and (thr_hit==False):
            #print(f'Optimal nPCs = {i+1}')
            #print(f'Optimal PCs index = {i}')
            #print(f'Deviation = {tmp_dev}')
            
            opt_n_pcs_linear_dev = i
            
            plt.scatter(x[i], y[i], s=15, color='red')
            plt.axvline(x[i], color='red', linestyle='--', linewidth=1)
            
            thr_hit=True
        
        plt.xlabel('PCs Ranked')
        plt.ylabel('Cumulative Variance')
        plt.title(f'Cumulative Explained Variance, i = {opt_n_pcs_linear_dev} | nPCs = {opt_n_pcs_linear_dev+1}')
        
        plt.ylim([y[0]*0.9, y[-1]*1.1])

    if output_dir is not None:
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        fn = f'{output_dir}/{fn}'
        plt.tight_layout()
        print(fn)
        plt.savefig(fn, dpi=400)

    return opt_n_pcs_linear_dev+1
        
# Function to test different linear deviation thresholds for optimizing nPCs

# We want to choose the deviation threshold that is the inflection point in the nPC -by- threshold plot. 
# The number of PCs is stable for a wide range of thresholds (representing another valid choice) and nPCs shoots up after becoming slightly more stringent. 
# We take the point right before the number of PCs begins to shoot up
def optimize_dev_thr(adata, pca_name, save=False):
    
    plt.style.use('default')
    if pca_name not in adata.uns.keys():
        exit(f'ERROR: {pca_name} not in adata.')
    
    x = range(len(adata.uns[pca_name]['variance_ratio']))
    y = np.cumsum(adata.uns[pca_name]['variance_ratio'])
    
    end_line_slope = (y[-1] - y[-2]) / (x[-1] - x[-2])
    end_line_int = y[-1] - end_line_slope*x[-1]
    
    deviation = []
    slopes = []
    dev_thr_list = [0.001, 0.0009, 0.0008, 0.0007, 0.0006, 0.0005, 0.0004, 0.0003, 0.0002,
                    0.0001, 0.00009, 0.00008, 0.00007, 0.00006, 0.00005, 0.00004, 0.00003, 0.00002, 
                    0.00001, 0.000009, 0.000008, 0.000007, 0.000006, 0.000005, 0.000004, 0.000003, 0.000002,
                    0.00000001, 0.000000001, 0.0000000001]
    
    thr_hit=False
    
    opt_n_pcs_linear_dev_list=[]
    
    for dev_thr in dev_thr_list:
        
        print(dev_thr)
        thr_hit=False
        
        for i in range(len(x)-2):
            
            tmp_slope = (y[i+1] - y[i]) / (x[i+1] - x[i])
            tmp_int = y[i] - tmp_slope*x[i]
            
            plt.scatter(x, y, s=10, color='blue')
            
            y_line = [tmp_slope*x + tmp_int for x in x]
            
            tmp_dev = abs(y_line[i+2] - y[i+2])
            
            deviation.append(tmp_dev)
            
            slopes.append(tmp_slope)
            
            plt.plot(x, y_line, color='grey', linewidth=0.5)
            
            if (tmp_dev < dev_thr) and (thr_hit==False):
                print(f'Optimal nPCs = {i+1}')
                print(f'Optimal PCs index = {i}')
                print(f'Deviation = {tmp_dev}')
                
                opt_n_pcs_linear_dev = i
                opt_n_pcs_linear_dev_list.append(i)
                
                plt.scatter(x[i], y[i], s=15, color='red')
                plt.axvline(x[i], color='red', linestyle='--', linewidth=1)
                
                thr_hit=True
        
        if thr_hit==False:
            print(f'Optimal nPCs = {len(x)}')
            print(f'Optimal PCs index = {len(x)-1}')
            print(f'Deviation = {tmp_dev}')
            
            opt_n_pcs_linear_dev = len(x)-1
            opt_n_pcs_linear_dev_list.append(len(x)-1)
            
            plt.scatter(x[i], y[i], s=15, color='red')
            plt.axvline(x[i], color='red', linestyle='--', linewidth=1)
            
            thr_hit=True
            
        plt.xlabel('PCs Ranked')
        plt.ylabel('Cumulative Variance')
        plt.title(f'Dev Thr = {dev_thr}, i = {opt_n_pcs_linear_dev} | nPCs = {opt_n_pcs_linear_dev+1}')
        
        plt.ylim([y[0]*0.9, y[-1]*1.1])

        if save:
            fn =  + f'{output_dir}/PCA_linear_deviation_optimization_dev_thr_{dev_thr}.png'
            plt.tight_layout()
            print(fn)
            plt.savefig(fn, dpi=400)

        plt.show()
    
    fig, ax = plt.subplots(figsize=(6,4))
    
    plt.plot(dev_thr_list, opt_n_pcs_linear_dev_list)
    plt.scatter(dev_thr_list, opt_n_pcs_linear_dev_list, color='red', s=10)
    
    ax.invert_xaxis()
    plt.ylabel('Optimal PCs')
    plt.xlabel('Deviation Threshold')

    if save:
        fn = f'{output_dir}/PCA_linear_deviation_optimization_dev_thr_stability.png'
        plt.tight_layout()
        print(fn)
        plt.savefig(fn, dpi=400)
    
    fig, ax = plt.subplots(figsize=(6,4))
    
    plt.plot(dev_thr_list[9:], opt_n_pcs_linear_dev_list[9:])
    plt.scatter(dev_thr_list[9:], opt_n_pcs_linear_dev_list[9:], color='red', s=10)
    
    ax.invert_xaxis()
    plt.ylabel('Optimal PCs')
    plt.xlabel('Deviation Threshold')

    if save:
        fn = f'{output_dir}/PCA_linear_deviation_optimization_dev_thr_stability_2.png'
        plt.tight_layout()
        print(fn)
        plt.savefig(fn, dpi=400)
    return pd.DataFrame({'dev_thresh':dev_thr_list,
                 'n_pcs':opt_n_pcs_linear_dev_list})
