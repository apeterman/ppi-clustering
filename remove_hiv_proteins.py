#remove the HIV proteins from the dataset (also removes human proteins if orphaned)
#! /usr/bin/python
import csv

origProteins = ["CA","GAG","GP120","GP41","GP160","IN","MA","NC","NEF","POL","P6","P1","PR","REV","RT","TAT","VIF","VPR","VPU"]
newProteins = []


with open ('/Users/andrewpetermen/Downloads/_apstuff/HIV_H_MSVM.csv','rU') as csvin:
	csvin = csv.reader(csvin, delimiter=',') 
	with open ('/Users/andrewpetermen/Downloads/_apstuff/HIV_H_MSVM_HUMAN.csv', 'wb') as csvout:
		csvout = csv.writer(csvout)

		for row in csvin:
			protein0 = row[0]
			protein1 = row[1]
			#print(protein2)

			if (protein0.upper() not in origProteins and protein1.upper() not in origProteins):
				csvout.writerow(row)

