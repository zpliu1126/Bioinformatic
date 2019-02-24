########编写一个判断sgRNA是否更靠近5‘端的awk脚本


#
#该脚本的输入文件格式为
#Ghir_A08G003390_A_80POT24841    - - - - - - - - -   0M      Ghir_A08        3394208 +       Ghir_A08G003390 Ghir_A08G003390 3392637 3396503 -
#也就是ot2gtf.pl脚本输出文件加上了gff文件中对应的 基因的编号 起始位点 终止位点 正负链
#
#
#输出文件的在输入文件最后面加上了距离基因5‘端的距离
#------------------------------------------------------------
########对该脚本得到的结果取最靠近基因5‘端的sgRNA

##awk -F "\t"  '{print $1"\t"substr($1,1,15)"\t"$13}' awk_outputfile|sort -k2,3|awk '{print $1"\t"$3"\t"$2}'|uniq -f2 >End_out_file
##########得到sgRNA的编号之后就利用grep去抓取输出文件就行了

{
	gene_start=$10;
	gene_end=$11;
	gene_tag=$12;
	sgRNA_position=$6;
	sgRNA_tag=substr($1,17,1);
}
{if(gene_tag=="+"&&sgRNA_tag=="S")
		print $0"\t"sgRNA_position-gene_start;
}
{if(gene_tag=="+"&&sgRNA_tag=="A")
		print $0"\t"sgRNA_position-gene_start-23;
}
{if(gene_tag=="-"&&sgRNA_tag=="S")
		print $0"\t"gene_end-sgRNA_position;
	}
{if(gene_tag=="-"&&sgRNA_tag=="A")
		print $0"\t"gene_end-sgRNA_position-23;
	}

	




