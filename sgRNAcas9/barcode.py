import pandas as pd
import os
import sys
import time
import re
'''
+ R1压缩文件路径
+ R2压缩文件路径
+ fastx_toolki软件斌目录
+ allBarcode文件路径
+ 结果文件路径
'''

def alter(file):
    '''
    替换文件中的head并且改变输出格式
    '''
    pattern = re.compile(r'^>')
    fasta_file = open(file, "r")
    # 获取头部内容
    head = fasta_file.readline()
    # 文件指针回到初始位置
    fasta_file.seek(0)
    # 将头部分成两部
    head_1 = ":".join(head.split("\n")[0].split(" ")[0].split(":")[0:4])+":"
    head_2 = head.split("\n")[0].split(" ")[1]
    file_data = "id\tseq\n"
    for line in fasta_file.readlines():
        if pattern.match(line):
            # 正则表达将字符串修改
            line = line.split("\n")[0].replace(
                head_1, '').replace(head_2, '\t')
        file_data += line
    with open(file, 'w') as f:
        f.write(file_data)


def rev(str):
    str = str.upper()[::-1]
    basecomplement = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G"
    }
    letters = list(str)
    letters = [basecomplement[base] for base in letters]
    return ''.join(letters)


'''
#arguments
+ R1 sequence file sml-3_FDDP190723711-1a_1.clean.fq.gz 
+ R2 sequence file sml-3_FDDP190723711-1a_2.clean.fq.gz
+ fastx software path
+ all barcode path 
+ out file path
'''
tmp_dir = "tmp"+time.strftime("%Y%m%d%a%H%M%S", time.localtime(time.time()))
R1_file = sys.argv[1]
R2_file = sys.argv[2]
fastx_toolki_path = sys.argv[3]  # ~/software/
all_barcode = sys.argv[4]
out_file = sys.argv[5]

os.system("mkdir "+tmp_dir)  # 存放临时文件的文件夹
print(">>Begin uncompress the file:"+R1_file)
os.system("gzip -d "+R1_file)
print(">>Begin uncompress the file:"+R2_file)
os.system("gzip -d "+R2_file)
print(">> transform fastq to fasta the file:"+os.path.splitext(R1_file)[0])
# os.path.splitext(R1_file)[0] 去除文件后缀的字符串
os.system(os.path.join(fastx_toolki_path, "fastq_to_fasta") + " -Q 33 -i "+os.path.splitext(R1_file)[0]+" -o "+tmp_dir+"/R1.fasta")
print(">> transform fastq to fasta the file:"+os.path.splitext(R2_file)[0])
os.system(os.path.join(fastx_toolki_path, "fastq_to_fasta") +" -Q 33 -i "+os.path.splitext(R2_file)[0]+" -o "+tmp_dir+"/R2.fasta")
print(">> reverse R2 file:"+os.path.splitext(R2_file)[0])
os.system(os.path.join(fastx_toolki_path, "fastx_reverse_complement") + "  -i  "+tmp_dir+"/R2.fasta"+" -o "+tmp_dir+"/RC_R2.fasta")
# 切换word dir
# print(">> Change work dir")
# os.system("cd "+tmp_dir)
# 读取文件第一行
# fasta_file = open("R1.fasta")
# R2_fasta = open("RC_R2.fasta")
## 去除头部信息
print(">>delete head")
alter(tmp_dir+"/R1.fasta")
alter(tmp_dir+"/RC_R2.fasta")

print(">>Merge R1 and R2")
# 合并R1与R2
a = pd.read_csv(tmp_dir+"/R1.fasta", sep="\t")
b = pd.read_csv(tmp_dir+"/RC_R2.fasta", sep="\t")
merge = pd.merge(left=a, right=b, left_on='id', right_on='id')
merge.to_csv(tmp_dir+"/R1_R2_oneline.fasta", index=False)


# pattern1=re.compile(r'[ATGC]{9}TATAAGCGAAAGAAGCATCAGATGGGCAAACAAAGCACCAGTGGTCTAGTGGTAGAATAGTACCCTGCCACGGTACAGACCCGGGTTCGATTCCCGGCTGGTGCA.*TAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTTTTAGCGCGTGCATGCCTGCAGGTCCACAAATTCGGGTC[ATGC]{6}')
# 获取带barcode序列的文件
print(">>filter sequence!")
pattern2 = re.compile(r'[ATCG]*([ATCG]{9}TATAAGCGAAAGAAGCATCAGATGGGCAAACAAAGCACCAGTGGTCTAGTGGTAGAATAGTACCCTGCCACGGTACAGACCCGGGTTCGATTCCCGGCTGGTGCA.*TAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTTTTAGCGCGTGCATGCCTGCAGGTCCACAAATTCGGGTC[ATGC]{6})[ATCG]*')
file_data = ''
with open(tmp_dir+"/R1_R2_oneline.fasta", 'r') as f:
    for line in f.readlines():
        # 获取seq序列
        line = (line.split(",")[1]+line.split(",")[2]).split("\n")[0]
        if pattern2.match(line):
            line = pattern2.match(line)[1]+"\n"
            file_data += line
        else:
            continue
# 将过滤后的数据写入文件
with open(tmp_dir+"/Filter.fasta", 'w') as f:
    f.write(file_data)

print(">>extract barcode!")
barcode = ''
# 提取barcode以及输出最终文件
with open(tmp_dir+"/Filter.fasta", 'r') as f:
    for line in f.readlines():
      line=line.split("\n")[0]
      barcode+=line[0:9]+rev(line)[0:6]+"\t"+line[114:134]+"\n"

with open(tmp_dir+"/barcode_sgRNA",'w') as f:
  f.write(barcode)

print(">>run bash script!") 
os.system("bash barcode.sh "+all_barcode+" "+out_file+" "+tmp_dir)
print("clean work dir")
os.system("rm -rf "+tmp_dir)
print("successful!")

