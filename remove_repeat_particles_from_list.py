#!/usr/bin/env python

import math, os, sys, random
try:
	from optparse import OptionParser
except:
	from optik import OptionParser

def main():
	(inlst1, num, center, euler_thres, output) =  parse_command_line()
	g = open(inlst1, "r")
        inlst_line1=g.readlines()
        PI=3.14159265359
        m_o=[0]*9
        m=[0]*9
        rm=[0]*9
        
        b=[0]*len(inlst_line1)
        gp=[[0 for x in range(len(inlst_line1))] for y in range(num)]
        number=[0]*num
        oo2=open(output,"w")
        name=inlst_line1[0].split('\t')[1]
        jj=0
        kk=0
        for i in range(0,len(inlst_line1)):
               if(inlst_line1[i].split('\t')[1]==name):
                    gp[jj][kk]=(inlst_line1[i])
                    kk=kk+1
               else:
                    number[jj]=kk
                    jj=jj+1
                    kk=0
                    name=inlst_line1[i].split('\t')[1]
                    gp[jj][kk]=(inlst_line1[i])
                    kk=kk+1
        number[jj]=kk
        for k in range(0,jj+1):
            a=[0]*(number[k])
            for i in range(0, number[k]):
               # num = int(inlst_line1[i].split('\t')[0])
               x1=float(gp[k][i].split('\t')[6].split('=')[1].split(',')[0])
               y1=float(gp[k][i].split('\t')[6].split('=')[1].split(',')[1])
               alt=float(gp[k][i].split('\t')[5].split('=')[1].split(',')[0])*PI/180
               az=float(gp[k][i].split('\t')[5].split('=')[1].split(',')[1])*PI/180
               phi=float(gp[k][i].split('\t')[5].split('=')[1].split(',')[2])*PI/180
               score1=float(gp[k][i].split('\t')[7].split('=')[1])
               if(a[i]==1):
                    continue
               for j in range(i+1,number[k]):
                    if((i+1)==number[k] or a[j]==1):
                           continue
                    x2=float(gp[k][j].split('\t')[6].split('=')[1].split(',')[0])
                    y2=float(gp[k][j].split('\t')[6].split('=')[1].split(',')[1])
                    dist=abs(math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)))
                    
                    score2=float(gp[k][j].split('\t')[7].split('=')[1])

                    if(dist<center):
                       alt_o=float(gp[k][j].split('\t')[5].split('=')[1].split(',')[0])*PI/180
                       az_o=float(gp[k][j].split('\t')[5].split('=')[1].split(',')[1])*PI/180
                       phi_o=float(gp[k][j].split('\t')[5].split('=')[1].split(',')[2])*PI/180
                    # origin eulerrmrix
                       m_o[0]=float(float(math.cos(phi_o))*float(math.cos(az_o))-float(math.sin(phi_o))*float(math.cos(alt_o))*float(math.sin(az_o)))
	               m_o[1]=math.cos(phi_o)*math.sin(az_o)+math.sin(phi_o)*math.cos(alt_o)*math.cos(az_o)
	               m_o[2]=math.sin(phi_o)*math.sin(alt_o)
	               m_o[3]=(-1)*math.sin(phi_o)*math.cos(az_o)-math.cos(phi_o)*math.cos(alt_o)*math.sin(az_o)
	               m_o[4]=(-1)*math.sin(phi_o)*math.sin(az_o)+math.cos(phi_o)*math.cos(alt_o)*math.cos(az_o)
	               m_o[5]=math.cos(phi_o)*math.sin(alt_o)
	               m_o[6]=math.sin(alt_o)*math.sin(az_o)
	               m_o[7]=(-1)*math.sin(alt_o)*math.cos(az_o)
	               m_o[8]=math.cos(alt_o)
         # calculated eulerrmrix
                       m[0]=math.cos(phi)*math.cos(az)-math.sin(phi)*math.cos(alt)*math.sin(az)
                       m[1]=(-1)*math.sin(phi)*math.cos(az)-math.cos(phi)*math.cos(alt)*math.sin(az)
                       m[2]=math.sin(alt)*math.sin(az)
                       m[3]=math.cos(phi)*math.sin(az)+math.sin(phi)*math.cos(alt)*math.cos(az)
                       m[4]=(-1)*math.sin(phi)*math.sin(az)+math.cos(phi)*math.cos(alt)*math.cos(az)
                       m[5]=(-1)*math.sin(alt)*math.cos(az)
                       m[6]=math.sin(phi)*math.sin(alt)
                       m[7]=math.cos(phi)*math.sin(alt)
                       m[8]=math.cos(alt)
         # relativermrix
                       rm[0]=m[0]*m_o[0]+m[3]*m_o[1]+m[6]*m_o[2]
                       rm[1]=m_o[0]*m[1]+m_o[1]*m[4]+m_o[2]*m[7]
                       rm[2]=m_o[0]*m[2]+m_o[1]*m[5]+m_o[2]*m[8]
                       rm[3]=m_o[3]*m[0]+m_o[4]*m[3]+m_o[5]*m[6]
                       rm[4]=m_o[3]*m[1]+m_o[4]*m[4]+m_o[5]*m[7]
                       rm[5]=m_o[3]*m[2]+m_o[4]*m[5]+m_o[5]*m[8]
                       rm[6]=m_o[6]*m[0]+m_o[7]*m[3]+m_o[8]*m[6]
                       rm[7]=m_o[6]*m[1]+m_o[7]*m[4]+m_o[8]*m[7]
                       rm[8]=m_o[6]*m[2]+m_o[7]*m[5]+m_o[8]*m[8]
         # to quaternion
                       if((rm[0]+rm[4]+rm[8])>0):
                         t=math.sqrt(1+rm[0]+rm[4]+rm[8])*2
                         if(t==0):
                             continue  #case 1
                         q0=0.25*t
                         q1=(rm[7]-rm[5])/t
                         q2=(rm[2]-rm[6])/t
                         q3=(rm[3]-rm[1])/t
                       if((rm[0]+rm[4]+rm[8])<=0):               #case 2
                          if(rm[0]>=rm[4] and rm[0]>=rm[8]):
                               t=math.sqrt(1+rm[0]-rm[4]-rm[8])*2
                               if(t==0):
                                   continue
                               q0=(rm[7]-rm[5])/t
                               q1=t/4
                               q2=(rm[2]+rm[6])/t
                               q3=(rm[1]+rm[3])/t
                          if(rm[4]>rm[8]):
                               if(1-rm[0]+rm[4]-rm[8]<0):
                                   continue
                               t=math.sqrt(1-rm[0]+rm[4]-rm[8])*2
                               if(t==0):
                                   continue
                               q0=(rm[2]-rm[6])/t
                               q1=(rm[1]+rm[3])/t
                               q2=t/4
                               q3=(rm[7]+rm[5])/t
                          else:
                               if((1-rm[0]-rm[4]+rm[8])<0):
                                   continue
                               t=math.sqrt(1-rm[0]-rm[4]+rm[8])*2
                               if(t==0):
                                   continue
                               q0=(rm[3]-rm[1])/t
                               q1=(rm[2]+rm[6])/t
                               q2=(rm[5]+rm[7])/t
                               q3=t/4
                       if(q0>1 or q0<-1):
                               continue
                       euler_dist=abs(math.acos(q0)*2*180/PI)
                       dist2=euler_dist%180
                       if(dist2>euler_thres+1):
                             dist2=180-dist2
                       if(euler_dist<euler_thres or dist2<euler_thres):
                             #print("score1="+str(score1)+"\t"+"score2="+str(score2))
                             if(score1>score2 and a[j]==0):
                                  a[j]=1 
                             if(score1<=score2 and a[i]==0):
                                  a[i]=1  
            for kk in range(0,number[k]):
               if(a[kk]==0):
                    oo2.write(gp[k][kk]) 
            a=[0]*len(inlst_line1)
               
        g.close()
        oo2.close()

	
def parse_command_line():
	usage="%prog <Detection file> <number of windows> <center thres> <euler thres> <output>"
	parser = OptionParser(usage=usage, version="%1")
	
	if len(sys.argv)<6: 
		print "<Detection file> <number of windows> <center thres> <euler thres> <output>"
		sys.exit(-1)
	
	(options, args)=parser.parse_args()
	
	inlst1 = args[0]
	num = int(args[1])
        center = float(args[2])
        euler_thres = float(args[3])
        output = args[4]
	return (inlst1, num, center, euler_thres, output)

if __name__== "__main__":
	main()
