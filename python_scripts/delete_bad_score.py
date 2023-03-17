#!/usr/bin/env python

import math, os, sys, random
try:
	from optparse import OptionParser
except:
	from optik import OptionParser

def main():
	(inlst1, thres, output) =  parse_command_line()
	g = open(inlst1, "r")
	num=0
	oo=open(output,"w")
	inlst_line1=g.readlines()
        k=3
        lst_line=len(inlst_line1[4].split())
        for j in range(2,(lst_line)):
            tmp=str(inlst_line1[4].split('\t')[j].split('=')[0])
            if(tmp=="score"):
	        index=int(j)

	for i in range(3,len(inlst_line1)):
                #cx = float(inlst_line1[i].split('\t')[6].split('=')[1].split(',')[0])
                #cy = float(inlst_line1[i].split('\t')[6].split('=')[1].split(',')[1])
                #score=abs(float(inlst_line1[i].split('\t')[7].split('=')[1]))
                score=abs(float(inlst_line1[i].split('\t')[index].split('=')[1]))
                #delta_x = abs(75-cx)
                #delta_y = abs(75-cy)
                if(score>thres):
                      oo.write(inlst_line1[i])
                                   
        g.close()
	oo.close()

	
def parse_command_line():
	usage="%prog <input1> <thres> <output>"
	parser = OptionParser(usage=usage, version="%1")
	
	if len(sys.argv)<3: 
		print "<input1> <thres> <output>"
		sys.exit(-1)
	
	(options, args)=parser.parse_args()
	
	inlst1 = args[0]
        thres = float(args[1])
        output = args[2]
	return (inlst1, thres, output)

if __name__== "__main__":
	main()
