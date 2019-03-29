#coding utf-8

'''
usage python geneSquence.py allgenesfile lookgenesfile outputfile
	-----------genes file----------------
	>Gbar_A01G000010.1
	ATGAGCTACCAATCAAACCCAGAGATTGAACTTCCTGGCTTTAGATTCCACCCAACGGAGGAAGAATTG
	TAGAATTTTACCTCAGAAACATCGTTTATGGCAAGATGTTGAGTTCTGATGTAATTGGATTTCTCAACAT
	-------------gene id-----------------"
	Gbar_A01G000010.1
	Gbar_A01G000110.1q
'''

import sys,getopt,re
from tqdm import tqdm

def usage():
	print("!! tqdm module is require s !!")
	print("usage:\n")
	print("\t"+"-G|--genes="+"\t"+"all genes file")
	print("\t"+"-q|--query="+"\t"+"query genes ID")
	print("\t"+"-o|--out="+"\t"+"output fasta file")
	print("\t"+"-h|--help"+"\t"+"help information")

try:
	opts,args=getopt.getopt(sys.argv[1:],"hG:q:o:",["help","genes=","query=","out="])
except getopt.GetoptError:
	print ("getopt moudle is error")
	sys.exit()
if len(opts)==3:
	for name,value in opts:
		if name in ("-h","--help"):
			usage()
			sys.exit()
		elif name in("-G","--genes"):
			genesfile=value
		elif name in("-q","--query"):
			querygenesfile=value
		elif name in("-o","--out"):
			outfile=value
else:
	print(">>>>>>>>>>>>>>>>>>>>>>>"+"Attention! lack parameters "+"<<<<<<<<<<<<<<<<<<<<<")
	usage()
	sys.exit()


with open(genesfile,'r') as allgenes:
	genes=allgenes.readlines()
with open(querygenesfile,'r') as lookgens:
	yourgenes=lookgens.readlines()



genelist={}

#pattern=re.compile("\.[0-9]([0-9])*")
for i in tqdm(range(0,len(genes)),desc="read genes file"):
	if re.match("^>",genes[i]):
		genes[i]=genes[i].strip("\n").strip(">")
		index1=i+1
		sequence=""
		while re.match("^[^>]",genes[index1]):
			sequence+=genes[index1].strip("\n")
			index1+=1
			#最后一行是序列，在less 过一行就不存在了
			if index1==len(genes):
				break
		genelist[genes[i]]=sequence
	else:
		continue

with open(outfile,'w') as out:
	for i in tqdm(range(0,len(yourgenes)),desc="write file"):
		yourgenes[i]=yourgenes[i].strip("\n")
		out.write(">"+yourgenes[i]+"\n"+genelist[yourgenes[i]]+"\n")



	