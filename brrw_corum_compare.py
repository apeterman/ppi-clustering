#! /usr/bin/python
# lookup uniprot ids and convert to gene name
# compare gene names in corum cluster (n) with brrw results (m)
# n x m comparisons, find overlaps in each comparison and return top match for each cluster

import csv
import re
import itertools

new_cluster_list =[]
fullresults = []
results = []

mapping = list(csv.reader(open('/Users/andrewpetermen/Downloads/_apstuff/uniprot_mapping_new.txt', 'rU'), delimiter='\t'))

brrw_clusters = list(csv.reader(open('/Users/andrewpetermen/Downloads/_apstuff/brrw_clusters.txt', 'rU'), delimiter=','))

corum_data = list(csv.reader(open('/Users/andrewpetermen/Downloads/_apstuff/corum_human_clusters.txt', 'rU'), delimiter='\t'))

#update uniprot with gene name in cluster column
for l in corum_data:
		new_list = []
		cluster_string = l[2].replace('(','').replace(')','')
		cluster_list = re.split(r',+', cluster_string)
		#print(cluster_list)

		for i in cluster_list:

			for t in mapping:
				if i == t[0]:
					new_list.append(t[1])
					break

		#print(new_list)
		new_cluster_list.append(new_list)
		#meta.append(l,new_cluster_list)
		#print(meta)

#search brrw clusters that have corum overlap for corum > 2 proteins
i = 0
print('Corum ID' + '|' + 'Corum Description' + '|' + 'Corum Cluster' + '|' + 'TBRRW Cluster' + '|' + 'Corum/TBRRW Overlap' + '|' + 'Precision' + '|' + 'Recall' + '|' + 'Fmeasure')

for cluster in new_cluster_list:
	
	if (len(cluster) > 2):
		
		#compare corum cluster to all brrw clusters
		result = [(filter(lambda x: x in cluster, sublist),sublist) for sublist in brrw_clusters]

		#remove empty resulting lists of overlaps
		result = [x for x in result if x[0] != []]

		#separate the overlapping protein list from the tbrrw cluster
		overlap = [x[0] for x in result]
		tbrrw = [x[1] for x in result]

		if result != []:
			#get best cluster by max overlap cluster length
			overlap = max(overlap,key=len)

			#find the matching tbrrw cluster
			tbrrw = [x[1] for x in result if set(x[0]) == set(overlap)]

			#if more than one, pick the smallest (most accurate)
			tbrrw = min(tbrrw,key=len)

			#calculate stats
			precision = len(overlap)*1.0/len(tbrrw)*1.0
			recall = len(overlap)*1.0/len(cluster)*1.0
			Fmeasure = 2*(precision*recall)/(precision+recall)

			print(corum_data[i][0] + '|' + corum_data[i][1] + '|' + ', '.join(cluster) + '|' + ', '.join(tbrrw) + '|' + ', '.join(overlap) + '|' + str(precision) + '|' + str(recall) + '|' + str(Fmeasure))

	i+=1
