"""
Integration of miRDP, miRDIP, and VarElect for miRNA target gene analysis.
Authors: 
Emily Smith (miRDB integration, Mapping combined results to VarElect)
Pax Bosner (miRDIP integration, VarElect API access)
"""


import requests
import csv
import urllib.request, urllib.parse
import pandas as pd
import mygene

# ------ miRDB Integration ------  (ES)
# miRNA Target Gene Analysis using miRDB locally stored data

# load miRDB locally 
miRDB = r"PATH TO FILE"   # REPLACE WITH PATH TO MIRDB FILE
df = pd.read_csv(miRDB, sep="\s+", usecols=['miRNA', 'GeneID', 'TargetScore']) # TargetScore changed from how file is downloaded (Target Score)
    
df['TargetScore'] = pd.to_numeric(df['TargetScore'], errors='coerce')
df = df.dropna(subset=['TargetScore'])


 # MiRNAs  
#      - Follow the notation: "hsa-miR-xxx-3/5p"
miRNAs = ("hsa-miR-122-5p", "hsa-miR-155-5p", "hsa-miR-140-3p") # miRNAs to search miRDB for
g = 10 # number of target genes to return for each miRNA

chosenmiRNA = df[df['miRNA'].isin(miRNAs)]
targets = []
for miRNA, group in chosenmiRNA.groupby('miRNA'):
    if miRNA in miRNAs:
        targets.append(group.nlargest(g, 'TargetScore'))
targets = pd.concat(targets, ignore_index=True)

NM = [(gene_id) for gene_id in targets['GeneID'].unique() if (gene_id).startswith('NM_')]
print(NM)

# Using mygene to convert NMXXXX Ids to gene symbols
mg = mygene.MyGeneInfo()
results = mg.querymany(NM, scopes='refseq', fields='symbol', species='human')
mapping = {item['query']: item.get('symbol', item['query']) for item in results}

targets['Gene Name'] = targets['GeneID'].map(mapping).fillna(targets['GeneID'])
targets = targets.drop('GeneID', axis=1) 
print(targets)

df_miRDB_final = targets.rename(columns={
    'miRNA': 'miRNA',
    'Gene Name': 'Gene Name',
    'TargetScore': 'miRNA-Gene Score',
})[['miRNA', 'Gene Name', 'miRNA-Gene Score']]
df_miRDB_final['Database'] = 'miRDB' # add database column
columns_final = df_miRDB_final[['miRNA', 'Gene Name', 'miRNA-Gene Score', 'Database']]
columns_final.to_csv('miRDB_Targets_Emily.csv', index=False) # SAVE TARGETS TO CSV - RENAME BASED ON NEED


gene_symbols = columns_final['Gene Name'].unique().tolist() # gene symbols (miRDB/mygene) for varelect
m = columns_final.values.tolist() # miRDIP Gene name column
c = columns_final['Database'].unique().tolist()  # confirm database


# ------ miRDIP Integration ------ (PB)

class mirDIP_Http:
    mapScore = {
        'Very High': '0',
        'High': '1',
        'Medium': '2',
        'Low': '3'
    }
    url = "http://ophid.utoronto.ca/mirDIP"
        
    map = {}  # results will be here

    def __init__(self):
        self.url = "http://ophid.utoronto.ca/mirDIP"
        return

    def unidirectionalSearchOnMicroRNAs(self, microRNAs, minimumScore):
        self.sendPost(self.url + "/Http_U", '', microRNAs, 
                      self.mapScore[minimumScore])
        return

    def sendPost(self, url_, geneSymbols, microrna, minimumScore, sources='', occurrances='1'):
        params = {
            'genesymbol': geneSymbols,
            'microrna': microrna,
            'scoreClass': minimumScore,
            'dbOccurrences': occurrances,
            'sources': sources
        }
        params = bytes(urllib.parse.urlencode(params).encode())
        
        try:
            handler = urllib.request.urlopen(url_, params)
        except Exception:
            self.print_exc()
        else:
            self.response = handler.read().decode('utf-8')
            self.makeMap()
        return

    def makeMap(self):     
        ENTRY_DEL = 0x01
        KEY_DEL = 0x02
        arr = self.response.split(chr(ENTRY_DEL))
            
        for str in arr:
            arrKeyValue = str.split(chr(KEY_DEL))
            if len(arrKeyValue) > 1: 
                self.map[arrKeyValue[0]] = arrKeyValue[1]
        return

    def getResuls(self): 
        if "results" in self.map: 
            return self.map["results"]
        else: 
            return ''


 # MiRNAs 
#      - Comma delimited. 
#      - Follow the notation: "hsa-miR-xxx-3/5p"
#      - Edit desired number of gene targets line 162

microRNAs = ("hsa-miR-122-5p,hsa-miR-155-5p,hsa-miR-140-3p")

 # Minimum Score 
#      - Choices: 'Very High', 'High', 'Medium', 'Low'. 
minimumScore = "Very High"

o = mirDIP_Http()  
o.unidirectionalSearchOnMicroRNAs(microRNAs, minimumScore)
dfs=[]
qqq=[]
scsv=o.getResuls()     
reader = csv.reader(scsv.split('\n'), delimiter=',')
for row in reader:
    qqq.append('\t'.join(row))
for q in range (1,(len(qqq)-1)):
    for w in range (8):
        line = qqq[q]
        values = line.split("\t")
        for ww in range (len(values)):
            values[ww].strip("'")
    dfs.append(values)
    
df= pd.DataFrame(dfs, columns=['Gene', 'Uniprot', 'Pseudogene','MicroRNA',
                               'IntegratedScore','Number of Sources',
                               'Score Class','Sources'])
df['IntegratedScore'] = pd.to_numeric(df.IntegratedScore, errors='coerce')
df= df.drop(labels=['Sources', 'Pseudogene', 'Number of Sources',
                    'Score Class'], axis=1)
miRna = microRNAs.split(",")
miRna = [nu.strip(' ') for nu in miRna]
l=[]
m=[]
b=[]
e=10     ##INPUT NUMBER OF TARGET GENES PER miRNA ***
 
###sort out number (e) of most relevant targets based on the score###
for k in range(len(miRna)):
    c = 0  
    for z in range(len(df.index)):
        if miRna[k] == df.iat[z, 2] and c < e:  
            b.append(miRna[k])  
            b.append(df.iat[z, 0]) 
            b.append(df.iat[z, 3])  
            l.append(df.iat[z, 0])  
            c += 1  
            m.append(b)  
            b = []  
        elif c >= e: 
            break
            
    
dv= pd.DataFrame (m, columns=['miRNA', 'Gene', 'miRNA-Gene Score'])
dv.index += 1
print(dv)
print("Saving DV")
dv.to_csv('researchResultsMirDIP.csv')  # rename based on need

# ------ Combine miRDB and miRDIP ------ (ES)

df_miRDB_final['Database'] = 'miRDB'  
df_miRDB_final = df_miRDB_final.rename(columns={
    'miRNA': 'miRNA',
    'Gene Name': 'Gene Name',
    'miRNA-Gene Score': 'Score'
})

dv['Database'] = 'miRDIP'  
dv = dv.rename(columns={'Gene': 'Gene Name', 'miRNA-Gene Score': 'Score'})
dv = dv[['miRNA', 'Gene Name', 'Score', 'Database']]  

# Combine miRDB and miRDIP results
combined = pd.concat([df_miRDB_final, dv], ignore_index=True)
combined = combined.sort_values(by=['miRNA', 'Score'], ascending=[True, False]) # Sort by miRNA & Score
final_combined = combined.groupby(['miRNA', 'Database']).head(e).reset_index(drop=True)  # Limit to set amount of gene targets per DB
final_combined.to_csv('Combined_miRDB_miRDIP_Final.csv', index=False)


# ------ VarElect Integration ------ (ES & PB)

# Use the combined gene names for VarElect
gene_symbols = final_combined['Gene Name'].unique().tolist()  # Extract unique gene names from combined results

# Prepare the VarElect query - API KEY REQUIRED FROM VARELECT GENECARDS
query = {
    'UserName': 'USERNAME', 
    'Key': 'KEY',
    'Phenotype': 'RV',  # Phenotype to map gene names to
    'Symbols': gene_symbols,
    'Customer': 'P',
    'RequestID': '1'
}


response = requests.post('https://ve.genecards.org/Api/ve/Analyze', json=query)
r = response.json()

dd = r.get('Data', {})
dire = dd.get('Direct', [])
for i in range(len(dire)):
    a = dire[i]
    if 'GiftScore' in a:
        a['GiftScore'] = float(a['GiftScore'])

varele = pd.DataFrame(dire, columns=['Rank', 'PValue', 'Symbol', 'Name', 'Category', 'GiftScore', 'Score', 'MatchedPhenotypes'])  
# varele = varele.drop(labels="Name", axis=1, errors='ignore') 
print(varele)
print("Saving VarElect to CSV")
varele.to_csv('VarElect_Results_Combined_Final.csv', index=False) # Change file name based on need
