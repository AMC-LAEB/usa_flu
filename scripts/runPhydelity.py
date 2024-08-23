import os
import pandas as pd
from ete3 import Tree
from multiprocessing.pool import ThreadPool


num  = 40
tp = ThreadPool(num)

pct = 0.1
time = 0.083
files = os.listdir('../time_trees_subsampled')


def runClusters(path, time, pct):
    if 'H3N2' in path:
        subtype = 'H3N2'
    if 'H1N1' in path:
        subtype = 'H1N1'
    if 'Yam' in path:
        subtype = 'Yam'
    if 'Vic' in path:
        subtype = 'Vic'


    dates = pd.read_csv("../../" + subtype + "_dates.csv")
    base_path = path[:-4]
    tree_path = "../../time_trees_subsampled/" + base_path + ".nwk"
    tree = Tree(tree_path,format = 5)
    leaves = [node.name for node in tree.get_leaves()]
    max_date = max(dates.loc[dates['Strain'].isin(leaves)]['Date'])

    call = ['phydelity_usa.py', "--tree", tree_path, "--wcl", str(1.5), "--date_last_sample", str(max_date), "--min_date_cluster", str(-10000),
    "--coal_time", str(time), "--coal_percentage", str(pct)]
	    
    os.system(" ".join(call))


dir = "time_"+str(time)+"_pct_"+str(pct)
if not os.path.isdir(dir):
    os.mkdir(dir)
os.chdir(dir)


for file in files:#

    tp.apply_async(runClusters, (file,time,pct))

tp.close()
tp.join()

