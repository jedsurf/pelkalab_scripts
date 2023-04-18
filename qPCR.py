import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot_mc_96(csv, title=None):
    mc = pd.read_csv(csv, encoding="latin_1", skiprows=20, header=0).set_index('Derivative').transpose()
    fig, ax = plt.subplots()
    plt.xlabel('T\u2098')
    ax.set_title(title)
    T = [int(i) for i in mc.index.tolist()[:-1]]
    for col in mc.columns:
        ax.plot(T, list(mc[col])[:-1])
    plt.xlim(min(T), max(T))
    ax.legend(mc.columns)
    plt.show()
    return ax
 
def read_ddct(filename, encoding='UTF-8'):
    ddct = pd.read_csv(filename, encoding=encoding, skiprows=18, header=0)
    ddct.dropna(axis=1, inplace=True, how='all')
    ddct.dropna(axis=0, inplace=True, subset=['GOI', 'Reference gene'], how='all')
    return ddct

def plot_ct_96(filename, title = None, hline = None):
    df = read_ddct(filename)
    df.replace("No Ct", 41.0, inplace=True)
    df['Ct GOI'] = df['Ct GOI'].astype(float)
    df['Ct Ref. gene'] = df['Ct Ref. gene'].astype(float)
    colors = dict(boxes='dimgrey', whiskers='k', medians='r', caps='k')
    boxproperties = dict(linestyle='-', linewidth=1.5,  facecolor='dimgrey', color='dimgrey')
    box = df.boxplot(column="Ct GOI", by="Sample name", rot=75, color=colors, patch_artist=True, boxprops=boxproperties, showfliers=True)
    if hline is not None:
        if isinstance(hline, list):
            for Ct in hline:
                plt.axhline(y=Ct, alpha=0.5, label="C\u209C = " + str(Ct))
        else:
            plt.axhline(y=hline, alpha=0.5, label="C\u209C = " + str(hline))
    if title is not None:
        plt.title(title)
    else:
        plt.title(filename)
    plt.suptitle('')
    plt.ylabel("C\u209C", size=20)
    plt.legend()
#     plt.ylim(0, 40) # would rather have it figure out what the max is, but Ct = 45 should be fine since we only run 45 cycles
    return box

def plot_dct_96(filename):
    df = read_ddct(filename)
    df.replace("No Ct", 41.0, inplace=True)
    df['dCt (Ref.Gen – GOI)'] = df['dCt (Ref.Gen – GOI)']*-1
    colors = dict(boxes='dimgrey', whiskers='k', medians='r', caps='k')
    boxproperties = dict(linestyle='-', linewidth=1.5,  facecolor='dimgrey', color='dimgrey')
    box = df.boxplot(column="dCt (Ref.Gen – GOI)", by="Sample name", rot=75, color=colors, patch_artist=True, boxprops=boxproperties, showfliers=True)
    plt.title(filename)
    plt.suptitle('')
    plt.ylabel("\u0394C\u209C (GOI - Ref. gene)", size=20)
    return box

# 384-well stuff
def read_ct_384(csv, replace_no_ct=46):
    ct = pd.read_csv(csv, encoding='latin_1')
    ct.dropna(axis=1, inplace=True, how='all')
    ct["Cq Mean"].replace(0, replace_no_ct, inplace=True)
    return ct


if __name__ == "__main__":
    GOI_list = ['GAPDH', np.nan, 'B2m', np.nan, 'LIF', np.nan, 	'CCL20', np.nan, 'CX3CL1', 	np.nan,
                'RNAse I', 	np.nan, 'STAT1', np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
                np.nan, np.nan, np.nan, np.nan, np.nan]
    cytokines = {
        'A': 'IL1β',
        'B': np.nan,
        'C': 'IL1β',
        'D': np.nan,
        'E': 'IL1β',
        'F': np.nan,
        'G': 'Unstimulated',
        'H': np.nan,
        'I': 'Unstimulated',
        'J': np.nan,
        'K': 'Unstimulated',
        'L': np.nan,
        'M': "RT (-) ctrl",
        'N': np.nan,
        'O': 'Water control',
        'P': np.nan
    }

def define_genes_cytokines(ct_df, genes_list, cytokines_dict):
    GOI_wells = []
    cytokines_wells = []
    for well in ct_df['Well'].values.tolist():
        GOI_wells.append(genes_list[int(well[1:3])-1])
        cytokines_wells.append(cytokines_dict.get(well[0]))
    ct_df['GOI'] = GOI_wells
    ct_df['Cytokine'] = cytokines_wells
    return ct_df

def get_ct_384(processed_df, replace_no_ct=46.000000):
    df = processed_df.loc[:, ["Well", "Cq Mean", "GOI", "Cytokine"]]
    df.replace(replace_no_ct, "No Ct", inplace=True)
    cytokine_df = df.loc[:, "Cytokine"]
    cytokine_df.drop_duplicates(inplace=True)
    df_list = []
    for cytokine in cytokine_df.values.tolist():
        df_list.append(df.loc[df["Cytokine"] == cytokine])
    ct384 = pd.concat(df_list)
    ct384.sort_values(by = "GOI", inplace=True)
    return ct384

def plot_ct_384(df, cytokine):
    Cytokines_wells = df[df["Cytokine"] == cytokine]
    colors = dict(boxes='dimgrey', whiskers='k', medians='r', caps='k')
    boxproperties = dict(linestyle='-', linewidth=1.5,  facecolor='dimgrey', color='dimgrey')
    box = Cytokines_wells.boxplot(column="Cq Mean", by="GOI", rot=75, color=colors, patch_artist=True, boxprops=boxproperties, showfliers=True)
    plt.title(cytokine)
    plt.suptitle('')
    plt.ylabel("C\u209C", size=20)
    plt.ylim(0, 45) # would rather have it figure out what the max is, but Ct = 45 should be fine since we only run 45 cycles
    plt.show()
    return box

def plot_ct_384_bygene(df, gene):
    genes_wells = df[df["GOI"] == gene]
    colors = dict(boxes='dimgrey', whiskers='k', medians='r', caps='k')
    boxproperties = dict(linestyle='-', linewidth=1.5,  facecolor='dimgrey', color='dimgrey')
    box = genes_wells.boxplot(column="Cq Mean", by="Cytokine", rot=75, color=colors, patch_artist=True, boxprops=boxproperties, showfliers=True)
    plt.title(gene)
    plt.suptitle('')
    plt.ylabel("C\u209C", size=20)
    plt.ylim(0, 45) # would rather have it figure out what the max is, but Ct = 45 should be fine since we only run 45 cycles
    return box

def plot_dct_384_bygene(df, gene):
    gene_ct = np.array(df[df["GOI"]==gene]["Cq Mean"].values.tolist())
    refgene_ct = np.array(df[df["GOI"]=="GAPDH"]["Cq Mean"].values.tolist())
    # Check that the refgene and GOI are in the same row for a given sample before calculating ∆Ct
    if df[df["GOI"]==gene]["Cytokine"].values.tolist() == df[df["GOI"]=="GAPDH"]["Cytokine"].values.tolist():
        dCt = np.subtract(gene_ct, refgene_ct)
        dCt_df = df.copy()[df["GOI"]==gene]
        dCt_df["∆Ct"] = dCt
    else:
        print("Incorrect plate layout!")
#     Coloring, etc.
    colors = dict(boxes='dimgrey', whiskers='k', medians='r', caps='k')
    boxproperties = dict(linestyle='-', linewidth=1.5,  facecolor='dimgrey', color='dimgrey')
#     Making the plot
    box = dCt_df.boxplot(column="∆Ct", by="Cytokine", rot=75, color=colors, patch_artist=True, boxprops=boxproperties, showfliers=True)
    plt.title(gene)
    plt.suptitle('')
    plt.axhline(y=0, c='k')
    plt.ylabel("∆C\u209C", size=20)
    return box

def process_mc_384(mc_csv, genes_list, cytokines_dict):
    mc384 = pd.read_csv(mc_csv)
    mc384.dropna(axis=1, inplace=True, how='all')
    mc384 = mc384.T
    mc384.reset_index(inplace=True)
    GOI_wells = []
    cytokines_wells = []
    for well in mc384['index'].values.tolist()[1:]:
        GOI_wells.append(genes_list[int(well[1:3]) - 1])
        cytokines_wells.append(cytokines_dict.get(well[0]))
    mc384['GOI'] = ['n/a'] + GOI_wells
    mc384['Cytokine'] = ['n/a'] + cytokines_wells
    mc384.dropna(axis=0, inplace=True, subset=['GOI', "Cytokine"], how='any')
    return mc384


def plot_mc_384(mc384, gene):
    T = mc384.iloc[0].values.tolist()[1:-2]
    genes_wells_mc = mc384[mc384["GOI"] == gene]
    fig, ax = plt.subplots()
    for well in genes_wells_mc.iloc[0:].values.tolist():
        ax.plot(T, well[1:-2])
    ax.set_title(gene)
    ax.legend(genes_wells_mc["index"].values.tolist())
    plt.show()
    return ax