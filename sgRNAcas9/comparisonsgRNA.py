#coding:utf-8
'''
读取基因的gff数据
[geneID:[flag,start,end]]
1.需要基因组gff注释文件
2.ot2gtf脚本处理并且过滤之后的文件
3.每个基因名字的长度信息
eg： Ghir_A09G006360 填15
4.sgRNA的结果文件 sgRNAcas9_report.xls
5.输出文件
'''
def usage():
	print("usage:\n")
	print("\t"+"-h|--help"+"\t"+"print help information")
	print("\t"+"-g|--gff="+"\t"+"gff file path way")
	print("\t"+"-s|--sgRNA="+"\t"+"sgRNA file path way")
	print("\t"+"-l|--genelength="+"\t"+"length of gene")
	print("\t"+"-r|--sequence="+"\t"+"sgRNA sequence path way")
	print("\t"+"-o|--outfile="+"\t"+"output file path way")

import sys,re,getopt

try:
	opts,args=getopt.getopt(sys.argv[1:],"hg:s:l:o:r:",["help","gff=","sgRNA=","genelength=","outfile=","sequence="])
except getopt.GetoptError:
	print("commd parameters is wrong!\n")
	sys.exit()

for name,value in opts:
	if name in ("-h","--help"):
		usage()
		sys.exit()
	if name in("-g","--gff"):
		gfffile=value
	if name in("-s","--sgRNA"):
		sgRNAfile=value
	if name in("-l","--genelength"):
		genelength=int(value)
	if name in("-o","--outfile"):
		outfile=value
	if name in("-r","--sequence"):
		sgRNAsquencefile=value


with open(gfffile,'r') as genegff:
	list1=genegff.readlines()

genelist={}
for i in range(0,len(list1)):
	list1[i]=list1[i].strip('\n')
	tmp=list1[i].split()
	if tmp[2]!="gene":
		continue
	else:
		#获取对应的基因编号
		genelist[tmp[8][3:18]]=[tmp[6],tmp[3],tmp[4]]
		
with open(sgRNAfile,'r') as sgRNA:
	list2=sgRNA.readlines()

sgRNAlist={}
###给每个基因赋一个初始值
#with open("11111111111",'w') as out:
for i in range(0,len(list2)):
	list2[i]=list2[i].strip('\n')
	tmp=list2[i].split()
	sgRNAlist[tmp[0][0:genelength]]=[0 for i in range(4)]
	# 获取靶向基因的信息
	if genelist[tmp[-1]][0]=="+":
		sgRNAlist[tmp[0][0:genelength]][0]=tmp[0]
		sgRNAlist[tmp[0][0:genelength]][1]=int(tmp[-3])-int(genelist[tmp[-1]][1])
		#out.write(tmp[0]+"\t"+tmp[-3]+"\t"+tmp[28]+"\t"+tmp[-1]+"\t"+genelist[tmp[0][0:genelength]][0]+"\t"+genelist[tmp[0][0:genelength]][1]+"\t"+genelist[tmp[0][0:genelength]][2]+"\n")
	else:
		sgRNAlist[tmp[0][0:genelength]][0]=tmp[0]
		sgRNAlist[tmp[0][0:genelength]][1]=int(genelist[tmp[-1]][2])-int(tmp[-3])
##找最靠近5'的那个sgRNA
for i in range(0,len(list2)):
	tmp=list2[i].split()
	if (genelist[tmp[-1]][0]=="+") and (int(tmp[-3])-int(genelist[tmp[-1]][1]))<sgRNAlist[tmp[0][0:genelength]][1]:
		sgRNAlist[tmp[0][0:genelength]][0]=tmp[0]
		sgRNAlist[tmp[0][0:genelength]][1]=int(tmp[-3])-int(genelist[tmp[-1]][1])
		###初始化第三第四位，用于放第二个sgRNA
	elif (genelist[tmp[-1]][0]=="-") and (int(genelist[tmp[-1]][2])-int(tmp[-3]))<sgRNAlist[tmp[0][0:genelength]][1]:
		sgRNAlist[tmp[0][0:genelength]][0]=tmp[0]
		sgRNAlist[tmp[0][0:genelength]][1]=int(genelist[tmp[-1]][2])-int(tmp[-3])
	else:
		continue

###找下一个间隔100bp~200bp的sgRNA
for i in range(0,len(list2)):
	tmp=list2[i].split()
	if  sgRNAlist[tmp[0][0:genelength]][0]==tmp[0]:
		continue
	elif (genelist[tmp[-1]][0]=="+") and (int(tmp[-3])-int(genelist[tmp[-1]][1])-sgRNAlist[tmp[0][0:genelength]][1])>=50 and (int(tmp[-3])-int(genelist[tmp[-1]][1])-sgRNAlist[tmp[0][0:genelength]][1])<=200:
		sgRNAlist[tmp[0][0:genelength]][2]=tmp[0]
		sgRNAlist[tmp[0][0:genelength]][3]=int(tmp[-3])-int(genelist[tmp[-1]][1])
	elif (genelist[tmp[-1]][0]=="-") and (int(genelist[tmp[-1]][2])-int(tmp[-3]))-sgRNAlist[tmp[0][0:genelength]][1]>=50 and (int(genelist[tmp[-1]][2])-int(tmp[-3])-sgRNAlist[tmp[0][0:genelength]][1])<=200:
		sgRNAlist[tmp[0][0:genelength]][2]=tmp[0]
		sgRNAlist[tmp[0][0:genelength]][3]=int(genelist[tmp[-1]][2])-int(tmp[-3])

###输出数据
# with open(sys.argv[3],'w') as out:
# 	for k in sgRNAlist:
# 		out.write(str(k)+"\t"+str(sgRNAlist[k][0])+"\t"+str(sgRNAlist[k][1])+"\t"+str(sgRNAlist[k][2])+"\t"+str(sgRNAlist[k][3])+"\n")

#提取sgRNA序列
#'Ghir_D06G001960_S_1\t3\t25\tGGCTTCTACGAGGAAAGATATGG\t23\t45.0 %\t#\t2\t0\t0\t5\t118\t195\t319\t#\t2\t0\t0\t0\t0\t1\t2\tRepeat_sites_or_bad?\n'
with open(sgRNAsquencefile,'r') as sgRNAseq:
	list3=sgRNAseq.readlines()

seq={}
for i in range(1,len(list3)):
	seq[list3[i].strip("\n").split("\t")[0]]=list3[i].strip("\n").split("\t")[3]

#反向互补函数
def complement(s):
	basecomplement={"A":"T","T":"A","G":"C","C":"G","a":"t","t":"a","g":"c","c":"g"}
	letters=[basecomplement[base] for base in s]
	return ''.join(letters)
##对应到每一个基因
# genename={}
# with open("./../Tf_CDS_and_Tf",'r') as name:
# 	list4=name.readlines()

# for i in range(0,len(list4)):
# 	if re.match("^>",list4[i]):
# 		genename[list4[i].strip("\n").split()[1][0:genelength]]=list4[i].strip("\n").split()[0]
# 	else:
# 		continue


#输出文件
with open(outfile,'w') as out:
	for k in sgRNAlist:
		# 仅仅只找到一个sgRNA
		if sgRNAlist[k][2]==0:
			seq1="Ttctagctctaaaac"+complement(seq[sgRNAlist[k][0]][0:20][::-1])+"tgcaccagccgggaat"
			out.write(">"+str(k)+"\t"+seq[sgRNAlist[k][0]]+"\n"+seq1+"\n")
		else:
			seq1="Ttctagctctaaaac"+complement(seq[sgRNAlist[k][0]][0:20][::-1])+"tgcaccagccgggaat"
			seq2="Ttctagctctaaaac"+complement(seq[sgRNAlist[k][2]][0:20][::-1])+"tgcaccagccgggaat"
			out.write(">"+str(k)+"\t"+seq[sgRNAlist[k][0]]+"\t"+seq1+"\t"+seq[sgRNAlist[k][2]]+"\t"+seq2+"\n")
