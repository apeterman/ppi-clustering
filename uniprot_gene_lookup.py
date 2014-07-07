#! /usr/bin/python
import urllib
import csv
import re

#read in manual mapping
mapping = list(csv.reader(open('/Users/andrewpetermen/Downloads/_apstuff/uniprot_mapping.txt', 'rU'), delimiter='\t'))

#write out lookup values
with open('/Users/andrewpetermen/Downloads/_apstuff/uniprot_mapping_new2.txt', 'w') as fp:
	a = csv.writer(fp, delimiter='\t')

	for u in mapping:
		code = u[0]
		#code = 'Q6FI13'
		data = urllib.urlopen("http://www.uniprot.org/uniprot/" + code + ".txt").read()
		lines = data.split('\n')

		for line in lines:
			if line[:2] == 'GN':
				#get name
				if (line.find('Name=') > 0):
					start = line.index('Name=') + 5
					gene = line[start:].split(';')[0]
					print (code + '\t' + gene)
					#a.writerows(list(code + '\t' + gene))
				# #get synonym(s)	
				if (line.find('Synonyms=') > 0):
					start = line.index('Synonyms=') + 9
					synonyms = line[start:].split(';')[0]
					synonym = synonyms.split(', ')
					for s in synonym:
						print(code + '\t' + s)			
						#a.writerows(list(code + '\t' + s))

		