import pandas as pd

SA_df = pd.read_excel('SADelta_5k_Delta_ref_5k_metadata_no_duplicates_n10582.xlsx', index_col=0)
SA_df=SA_df[SA_df['country']=='South Africa'].reset_index()

print(SA_df)

print(SA_df.columns)

phylo_df=pd.read_excel('Phylotype_results_country.xlsx',index_col=0)
phylo_df=phylo_df.dropna()
phylo_df=phylo_df.reset_index()
print(phylo_df)

print(phylo_df.columns)

def strains_list(x):
    seq_list=x.split()
    seq_list_final=[]
    for i in range(0,len(seq_list)):
        seq=seq_list[i].split(' ')
        #seq_name=seq[0]+'/'+seq[1]+'/'+seq[5]
        seq_list_final.append(seq)
    return seq_list_final

def define_cluster(x):
    cluster=str(int(x))
    return cluster

phylo_df['seq_list']=phylo_df['strains'].apply(lambda x: strains_list(x))
phylo_df['ID']=phylo_df['ID'].apply(lambda x: define_cluster(x))
print(phylo_df)

cluster_id_list=phylo_df['ID'].tolist()
print(cluster_id_list)


def flatten(x):
    flat_list = [item for sublist in x for item in sublist]
    return flat_list

cluster_df=pd.DataFrame()

cluster_df['cluster_ID']=phylo_df['ID']
cluster_df['sequence_list']=phylo_df['seq_list']
cluster_df['sequence_list']=cluster_df['sequence_list'].apply(lambda x: flatten(x))
cluster_df=cluster_df.set_index('cluster_ID')
print(cluster_df)


cluster_dict=cluster_df.to_dict()
#print(cluster_dict['sequence_list']['8922'])

def classify_cluster(x):
    print(x)
    y=[]
    for i in range(0, len(cluster_id_list)):
        if x in cluster_dict['sequence_list'][cluster_id_list[i]]:
            y=cluster_id_list[i]
            return y
        else:
            continue
        y=[]

def num_clusters(x):
    try:
        return len(x)
    except:
        print('oops')

SA_df['cluster']=SA_df['strain'].apply(lambda x: classify_cluster(x))
#SA_df['cluster_num']=SA_df['cluster'].apply(lambda x: num_clusters(x))
print(SA_df['cluster'].tolist())
#print(SA_df['cluster_num'].tolist())

SA_df.to_csv('SADelta_5k_Delta_ref_5k_metadata_no_duplicates_n10582_clusters.csv',header=True)
