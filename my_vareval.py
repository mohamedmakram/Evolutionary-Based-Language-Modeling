import pandas as pd
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import pairwise2
from Bio import SeqIO
from Bio import SearchIO
import os 
import shutil
import random
import math
from typing import Dict, List, AnyStr
import argparse



parser = argparse.ArgumentParser(description='Generate evaluation report for generared vs natural sequences.')
parser.add_argument('-n','--nat', help='path to natural fasta sequences file', required=True)
parser.add_argument('-g','--gen', help='path to generated fasta sequences file', required=True)
parser.add_argument('-d','--db', help='path to cath database file .. something like /path/to/cath-S35-hmm3-v4_1_0.lib', required=True)
parser.add_argument('-m','--hmmer', help='path to hmmer binary, ', required=True)

args = vars(parser.parse_args())

def fasta2dict (infile:AnyStr) -> Dict[AnyStr, List[AnyStr]]:
    """Parse fasta file and convert it into dict of {'ID' : [all ids], 'Seq' : [all sequences]}

    Args:
        infile (AnyStr): path to fasta file with >= 1 sequences and ids.

    Returns:
        Dict[AnyStr, AnyStr]: dict of the format 'ID' : [all ids], 'Seq' : [all sequences].
    """
    tmp = SeqIO.parse(infile,'fasta')
    fasta_dict = {'ID':[],'Seq':[]}
    for record in tmp:
        fasta_dict['ID'].append(record.id)
        fasta_dict['Seq'].append(str(record.seq))
    return fasta_dict  
  
## Function to filter duplicated sequences in fasta file ##
## Outputs dictionary for the filtered seqs and path of the filtered fasta file
def filter_duplicates(infile:AnyStr, natural_dict: Dict[AnyStr, List[AnyStr]] , verbose:bool=True):
    if verbose: 
        print('Filtering duplicates...\n')
   
    in_dict = fasta2dict(infile)
    filter_dict={'ID':[], 'Seq':[]}
    for i in range(len(in_dict['ID'])):
        if (in_dict['Seq'][i] not in filter_dict['Seq']) and (in_dict['Seq'][i] not in natural_dict['Seq']):
            filter_dict['ID'].append(in_dict['ID'][i])
            filter_dict['Seq'].append(in_dict['Seq'][i])

    filename = infile.split('/')[-1]
    filename = '.'.join(filename.split('.')[:-1]) 
    outname = filename+'_filtered.fasta'
    filtered_fasta = open(filename+'/'+outname,'w')

    c=1
    for i in range(len(filter_dict['ID'])):
        filtered_fasta.write(">"+filter_dict['ID'][i])
        filtered_fasta.write('\n')
        filtered_fasta.write(filter_dict['Seq'][i])
        filtered_fasta.write('\n')
        c=c+1
    filtered_fasta.close()
  
    return [filter_dict, outname]
  
## Function to get identity scores between natural and generated seqs ##
def nat_gen_ident(nat_dict,gen_dict,filename):
    print('Calculating pairwise identities between natural and generated seqs..')
    
    Ident_dict = {}
    c = 0

    for gen_seq in gen_dict['Seq']:
        ind = 0
        IDs = ''
        scores=[]
        for nat_seq in nat_dict['Seq']:
            alns = pairwise2.align.globalxx(nat_seq, gen_seq,score_only=True)
            score = (alns/len(nat_seq))*100
            scores.append(score)
  
        ind = np.argmax(scores)
        IDs = str(gen_dict['ID'][c])+ ' : '+ str(nat_dict['ID'][ind])
        Ident_dict[IDs] = max(scores)
        c = c+1
        ind=0

    x=[]
    for i in Ident_dict:
        x.append(Ident_dict[i])

    sns.set(style="darkgrid")

    sns.histplot(data=x, color="blue",kde=True)
    print('Plotting..\n')
    plt.xlabel('Identity, %')
    plt.ylabel('sequences')
    plt.title('Generated Against Natural')
    #plt.show()
    plt.savefig(filename+'/gen_vs_nat.png')
    plt.clf()
    return Ident_dict
  
## A function to plot diversity Histograms ##
def diversity_histo(nat_dict, gen_dict,filename):
    print('Calculating internal pairwise identities for generated seqs..')
    gen_iden = []
    c = 0

    for i in range(len(gen_dict['Seq'])):
        scores=[]
        for j in range(len(gen_dict['Seq'])):
            if i != j:
                alns = pairwise2.align.globalxx(gen_dict['Seq'][i], gen_dict['Seq'][j],score_only=True)
                score = (alns/len(gen_dict['Seq'][i]))*100
                scores.append(score)
            else:
                score = 0
        gen_iden.append(max(scores))
        c = c+1
    
    nat_iden = []
    c = 0
    print('Calculating internal pairwise identities for natural seqs..')
    for i in range(len(natural_dict['Seq'])):
        scores=[]
        for j in range(len(natural_dict['Seq'])):
            if i != j:
                alns = pairwise2.align.globalxx(natural_dict['Seq'][i], natural_dict['Seq'][j],score_only=True)
                score = (alns/len(natural_dict['Seq'][i]))*100
                scores.append(score)
            else:
                score = 0
        nat_iden.append(max(scores))
        c = c+1

    sns.set(style="darkgrid")
    sns.histplot(data=gen_iden, color="red", label="Generated seq", kde=True)
    sns.histplot(data=nat_iden, color="skyblue", label="Natural seq", kde=True)
    print('Plotting..')
    plt.xlabel('Identity, %')
    plt.ylabel('sequences')
    plt.title('Identities Distributions')
    plt.legend() 
    #plt.show()
    plt.savefig(filename+'/internal_identities.png')
    plt.clf()

## A function to generate MSA files from FASTA files ##  
def msa(fasta_path,filename):
    print('Generating Multiple Sequence Alignment files..')
    cwd = os.getcwd()
    outpath= fasta_path.split('/')[-1]
    outpath = outpath.split('.')[0]
    outpath = outpath+'_msa.fasta'
    home =os.path.expanduser('~')
    os.chdir(home)
    os.system('chmod u+x clustalo')
    cmd = 'clustalo -i "'+fasta_path+ '" --iter=2 --full-iter'+' -o '+outpath+' -v'
    os.system(cmd)
    os.replace(outpath, cwd+'/'+filename+'/'+outpath)
    os.chdir(cwd)
    return(cwd+'/'+filename+'/'+outpath)

## Shannon Entropy Functions ##
def parseMSA(msa_file, alnformat="fasta"):
    """Parse in the MSA file using Biopython's AlignIO"""
    print('Parsing MSA...')
    from Bio import AlignIO
    alignment = AlignIO.read(msa_file, alnformat)

    # Do a little sanity checking:
    seq_lengths_list = []
    for record in alignment:
        seq_lengths_list.append(len(record))

    seq_lengths = set(seq_lengths_list)

    print("Alignment length is:" + str(list(seq_lengths)))

    if len(seq_lengths) != 1:
        sys.stderr.write("Your alignment lengths aren't equal. Check your alignment file.")
        sys.exit(1)

    index = range(1, list(seq_lengths)[0]+1)

    return alignment, list(seq_lengths), index
def shannon_entropy(list_input, normalized=False):
    """Calculate Shannon's Entropy per column of the alignment (H=-\sum_{i=1}^{M} P_i\,log_2\,P_i)"""
    
    import math
    unique_base = set(list_input)
    M = len(list_input)
    entropy_list = []
    # Number of residues in column
    for base in unique_base:
        n_i = list_input.count(base) # Number of residues of type i
        P_i = n_i/float(M) # n_i(Number of residues of type i) / M(Number of residues in column)
        entropy_i = P_i*(math.log(P_i,2)) if not normalized else P_i*(math.log(P_i,M))
        entropy_list.append(entropy_i)

        sh_entropy = -(sum(entropy_list))

    return sh_entropy
def shannon_entropy_list_msa(alignment, drop_col_if_gaps_more_than:float=None, normalized=False):
    """Calculate Shannon Entropy across the whole MSA"""
    # alginment is a 2D array where each row is a sequence, each column is position in the msa.
    print('Calculating shannon entropies..')
    shannon_entropy_list = []
    for col_no in range(len(list(alignment[0]))):

        list_input = list(alignment[:, col_no])  # this column from all sequences

        if drop_col_if_gaps_more_than is not None:
            if len([x for x in list_input if x == '-'])/len(list_input) > drop_col_if_gaps_more_than:
                continue
        
        shannon_entropy_list.append(shannon_entropy(list_input, normalized=normalized))

    return shannon_entropy_list
def running_mean(l, N):
    sum = 0
    result = list(0 for x in l)

    for i in range( 0, N ):
        sum = sum + l[i]
        result[i] = sum / (i+1)

    for i in range( N, len(l) ):
        sum = sum - l[i-N] + l[i]
        result[i] = sum / N

    return result
def plot(shannon_df,filename):
    """"Create a quick plot via matplotlib to visualise"""
    print("Plotting data...")
    plt.figure(figsize=(15,8))
    sns.set(style="darkgrid")
    sns.lineplot(data=shannon_df)

    print('Plotting...\n')
    plt.xlabel('MSA Position Index', fontsize=16)
    plt.ylabel('Shannon Entropy', fontsize=16)
    plt.title('Positional Shannon Entropies')
    plt.legend() 
    plt.savefig(filename+'/'+'Shannon_entropies.png')
  
def random_combinations (natural_fasta, filtered_fasta, filename):
    filtered_dict = fasta2dict(filtered_fasta)
    natural_dict = fasta2dict(natural_fasta)
    print('Generating random combinations fasta files...')
    cwd = os.getcwd()
    outpath = cwd+'/'+filename+'/cath_search'
    os.mkdir(outpath)

    num_files = 0
    natural_seqs = len(natural_dict['Seq'])
    generated_seqs = len(filtered_dict['Seq'])
    large_dataset = dict()
    small_dataset = dict()
    if (natural_seqs > generated_seqs):
        large_dataset = dict(zip(natural_dict['ID'], natural_dict['Seq']))
        small_dataset = dict(zip(filtered_dict['ID'], filtered_dict['Seq']))
        if ((natural_seqs%generated_seqs) == 0):
            num_files = natural_seqs/generated_seqs
        else:
            num_files = int(math.floor(natural_seqs/generated_seqs)+1)
    else:
        large_dataset = dict(zip(filtered_dict['ID'], filtered_dict['Seq']))
        small_dataset = dict(zip(natural_dict['ID'], natural_dict['Seq']))
        if ((generated_seqs%natural_seqs)==0):
            num_files = generated_seqs/natural_seqs
        else:
            num_files = int(math.floor(generated_seqs/natural_seqs)+1)
            #print(num_files, natural_seqs,generated_seqs)
    ## Shuffling ##
    large_key_list = [*large_dataset]
    #print(large_key_list)
    random.seed(5)
    random.shuffle(large_key_list)
    large_shuffled = dict()
    for key in large_key_list:
        large_shuffled[key] = large_dataset[key]

    small_key_list = [*small_dataset]
    random.seed(5)
    random.shuffle(small_key_list)
    small_shuffled = dict()
    for key in small_key_list:
        small_shuffled[key] = small_dataset[key]

    counter = 0
    paths_list=[]
    for f in range (int(num_files)):
        c=0
        paths_list.append(cwd+'/'+filename+'/cath_search'+'/combination_'+str(f)+'.fasta')
        with open (filename+'/cath_search'+'/combination_'+str(f)+'.fasta','w') as out_fasta:
            for id_small in small_shuffled:
                out_fasta.write('>'+id_small+'\n') 
                out_fasta.write(small_shuffled[id_small]+'\n')
            while c <= len([*small_dataset]):
                if counter < len([*large_dataset]):
                    out_fasta.write('>'+[*large_shuffled][counter]+'\n')
                    out_fasta.write(large_shuffled[[*large_shuffled][counter]]+'\n')
                
                else:
                    counter = 0
                    out_fasta.write('>'+[*large_shuffled][counter]+'\n')
                    out_fasta.write(large_shuffled[[*large_shuffled][counter]]+'\n')
                counter = counter+1
                c=c+1
    return paths_list
  
def cath_search (paths_list, cath_database, hmmer_path, filename):
    print('Looking for functional domains on CATH database...')
    cwd = os.getcwd()
    outpath = cwd+'/'+filename+'/cath_search'
    i=0
    parsed_paths = []
    for path in paths_list:  
        os.chdir(hmmer_path+'/src')  
        search_cmd = './hmmscan --domtblout ${outpath}/combination_${i}_domain.out --tblout ${outpath}/combination_${i}_seq.out -E 1e-6 --domE 1e-6 ${cath_database} ${path} > ${outpath}/combination_${i}.log'
        search_cmd = search_cmd.replace('${outpath}',outpath)
        search_cmd = search_cmd.replace('${i}',str(i))
        search_cmd = search_cmd.replace('${cath_database}', cath_database)
        search_cmd = search_cmd.replace('${path}',path)
        os.system(search_cmd)
        os.chdir(cwd)

        parse_cmd = 'bash hmmscan-parser.sh ${outpath}/combination_${i}_domain.out > ${outpath}/combination_${i}_final.txt'
        parse_cmd = parse_cmd.replace('${outpath}', outpath)
        parse_cmd = parse_cmd.replace('${i}',str(i))
        os.system(parse_cmd)
        
        parsed_paths.append(outpath+'/combination_'+str(i)+'_final.txt')
        i=i+1
    return parsed_paths

def parse_cath (parsed_paths,filename):
    tmp_df = pd.DataFrame()
    df = pd.read_csv(parsed_paths[0],sep='\t',header=None)
    print('Parsing CATH output..')
    for path in parsed_paths[1:]:
        tmp_df=pd.read_csv(path,sep='\t',header=None)
        df = df.append(tmp_df)
    
    cath_df = pd.DataFrame()
    cath_df['target_ID']= df[0]
    cath_df['t_len']= df[1]
    cath_df['query_ID']= df[2]
    cath_df['query_len']=df[3]
    cath_df['i-evalue']=df[4]
    cath_df['hmm_from']=df[5]
    cath_df['hmm_to']=df[6]
    cath_df['ali_from']=df[7]
    cath_df['ali_to']=df[8]
    cath_df['acc']= df[9]
    cath_df = cath_df.drop_duplicates(subset=['query_ID'])

    nat_list = []
    gen_list = []
    query_list = list(cath_df['query_ID'].values)
    for query in query_list:
        if query.startswith('g'):
            gen_list.append(query)
        elif query.startswith('n'):
            nat_list.append(query)
        else:
            raise ValueError(f'Sequence id shoudl start with either "g" for generated sequences of "n" for natural sequences ..got : {query}')

    df_nat = cath_df[cath_df['query_ID'].isin(nat_list)]
    df_gen = cath_df[cath_df['query_ID'].isin(gen_list)]
    df_nat.to_csv(filename+'/'+'natural_cath_summary.csv',index=False)
    df_gen.to_csv(filename+'/'+'generated_cath_summary.csv',index=False)
    print('Writing summaries to csv..\n')
    domains_nat = dict()
    for t in list(df_nat['target_ID'].values):
        domains_nat[t] = 0
    for t in list(df_nat['target_ID'].values):
        domains_nat[t] = domains_nat[t]+1
    
    domains_gen = dict()
    for t in list(df_gen['target_ID'].values):
        domains_gen[t] = 0
    for t in list(df_gen['target_ID'].values):
        domains_gen[t] = domains_gen[t]+1

    return [domains_nat,domains_gen]


## Input fasta processing ##
generated_fasta = args['gen']
#natural_fasta = input('Enter the path of the natural seqs FASTA: ')
natural_fasta = args['nat']
## Creating dictionnaries ##
generated_dict = fasta2dict(generated_fasta)
natural_dict = fasta2dict(natural_fasta)

filename = generated_fasta.split('/')[-1]
filename = '.'.join(filename.split('.')[:-1]) 
os.mkdir(filename)
## Filtering duplicates ##
filtered_dict = filter_duplicates(generated_fasta,natural_dict)[0] ##new dictionary for filtered sequences
cwd=os.getcwd()
filtered_fasta =cwd+'/'+filename+'/'+ filter_duplicates(generated_fasta,natural_dict)[1] ##path of the filtered sequences FASTA file

## Natural Vs Generated identity scores ##
nat_gen_identities = nat_gen_ident(natural_dict,filtered_dict,filename)
## Diversity histograms ##
diversity_histo(natural_dict,filtered_dict,filename)

## Generate MSA files ##
with open(filename+'/gen_nat.fasta','w') as wfd:
    for f in [natural_fasta,filtered_fasta]:
        with open(f,'r') as fd:
            shutil.copyfileobj(fd, wfd)
        wfd.write('\n')
combined_fasta = cwd+'/'+filename+'/gen_nat.fasta'
combined_msa = msa(combined_fasta,filename)

gen_msa = cwd+'/'+filename+'/gen_msa.fasta'
nat_msa = cwd+'/'+filename+'/nat_msa.fasta'
nat_copy_bool = True
with open (combined_msa,'r') as comb:
    with open (gen_msa,'w') as gen:
        with open (nat_msa,'w') as nat:
            for line in comb:
                if line.startswith('>'):
                    if 'g_' in line:
                        nat_copy_bool = False
                    elif 'n_' in line:
                        nat_copy_bool = True
                    else:
                        raise ValueError(f"Sequence id must have either 'g_' for generated of 'n_' for nature sequences.. got: {line}")
                if nat_copy_bool == True: 
                    nat.write(line)
                else:
                    gen.write(line)

##Shannon Entropies ##
alnformat = "fasta" # args.alnformat
makeplot = True # args.makeplot

    # Start calling functions to do the heavy lifting

generated_alignment, generated_seq_lengths, generated_index = parseMSA(gen_msa, alnformat)
generated_sel = shannon_entropy_list_msa(generated_alignment)
generated_sel = running_mean(generated_sel, 1)

natural_alignment, natural_seq_lengths, natural_index = parseMSA(nat_msa, alnformat)
natural_sel = shannon_entropy_list_msa(natural_alignment)
natural_sel = running_mean(natural_sel, 1)

index=0
if generated_seq_lengths[0] >= natural_seq_lengths[0]:
    index=generated_index
    for i in range(generated_seq_lengths[0]-natural_seq_lengths[0]):
        natural_sel.append(np.NaN)
else:
    index=natural_index
    for i in range(natural_seq_lengths[0]-generated_seq_lengths[0]):
        generated_sel.append(np.NaN)

shannon_df = pd.DataFrame()
shannon_df['generated']=generated_sel
shannon_df['natural']=natural_sel
plot(shannon_df,filename)

out_df = pd.DataFrame()
out_df['Position'] = index
out_df['gen_entropy'] = shannon_df['generated']
out_df['nat_entropy']=shannon_df['natural']
out_df.to_csv(filename+'/'+'shannon_entropies.csv',index=False)

shannon_df=shannon_df.dropna()
#shannon_df.to_csv('shannon_nona.csv', index=False)
pearson=stats.pearsonr(shannon_df['generated'],shannon_df['natural'])
mse = np.square(np.subtract(shannon_df['generated'],shannon_df['natural']).mean())

## CATH SEARCH ##
#cath_database = input('Enter the path to the CATH domains database: ')
#hmmer_path = input('Enter the path to the HMMER3 directory: ')
cath_database = args['db']
hmmer_path = args['hmmer']
natural_fasta=os.path.abspath(natural_fasta)
combination_paths = random_combinations(natural_fasta,filtered_fasta,filename)
parsed_paths = cath_search(combination_paths,cath_database,hmmer_path,filename)
all_domains = parse_cath(parsed_paths,filename)
nat_domains = all_domains[0]
gen_domains = all_domains[1]


## Output Report ##
print('Writing stats report...')
with open (filename+'/'+'stats_report.txt','w') as report:
  report.write('GENERAL STATISTICS\n')
  report.write('Number of natural sequences = '+str(len(natural_dict['Seq']))+'\n')
  report.write('Number of generated sequences = '+str(len(generated_dict['Seq']))+'\n')
  report.write('Number of unique generated sequences = '+str(len(filtered_dict['Seq']))+'\n')
  report.write('-------\n\nSHANNON ENTROPY STATISTICS\n')
  report.write('Pearson correlation r-value = '+str(pearson[0])+'\n')
  report.write('P-value = '+str(pearson[1])+'\n')
  report.write('Mean Square Error (MSE) = '+str(mse)+'\n')
  report.write('-------\n\nCATH FUNCTIONAL DOMAINS SEARCH SUMMARY\n')
  report.write('\n[1]  NATURAL SEQUENCES: \n')
  for domain in nat_domains:
    report.write('*'+str(nat_domains[domain])+' sequences matched with cath domain: '+domain+'\n')

  report.write('\n[2]  GENERATED SEQUENCES: \n')
  for domain in gen_domains:
    report.write('*'+str(gen_domains[domain])+' sequences matched with cath domain: '+domain+'\n')
print('DONE\n------------------')


