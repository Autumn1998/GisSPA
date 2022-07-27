#!/usr/bin/env python

import math, os, sys, random
try:
	from optparse import OptionParser
except:
	from optik import OptionParser

def main():
	(inlst1, inlst2, scale, clip, apix, output) =  parse_command_line()
	g = open(inlst1, "r")
        f = open(inlst2, "r")
	num=0
	oo=open(output,"w")
	inlst_line1=g.readlines()
        inlst_line2=f.readlines()
        k=3
        mag=50000/apix
        oo.write("# RELION; version 3.0-beta-2\n\ndata_\n\nloop_\n _rlnMicrographName #1 \n_rlnCoordinateX #2 \n_rlnCoordinateY #3 \n_rlnImageName #4 \n_rlnDefocusU #5 \n_rlnDefocusV #6 \n_rlnDefocusAngle #7 \n_rlnVoltage #8 \n_rlnSphericalAberration #9 \n_rlnAmplitudeContrast #10 \n_rlnMagnification #11 \n_rlnDetectorPixelSize #12 \n_rlnCtfFigureOfMerit #13 \n_rlnGroupNumber #14 \n_rlnAngleRot #15 \n_rlnAngleTilt #16 \n_rlnAnglePsi #17 \n_rlnOriginX #18 \n_rlnOriginY #19 \n_rlnClassNumber #20 \n_rlnNormCorrection #21 \n_rlnLogLikeliContribution #22 \n_rlnMaxValueProbDistribution #23 \n_rlnNrOfSignificantSamples #24\n")
	for i in range(0,len(inlst_line1)):
               # num = int(inlst_line1[i].split('\t')[0])
               # m = int(inlst_line1[i].split('\t')[1].split('.')[1])
                #m=9
                if(1):
                    nn = inlst_line1[i].split('\t')[0]
                    s1 = inlst_line1[i].split('\t')[1].split('/')[1].split('.hdf')[0]
                    micro=s1+".mrc"
                    #if(os.path.exists(output)):
                    #      oo=open(output,"a")   
                    
                    df = float(inlst_line1[i].split('\t')[2].split('=')[1])*10000
                    dfdiff = float(inlst_line1[i].split('\t')[3].split('=')[1])*10000
                    dfu=df-dfdiff
                    dfv=df+dfdiff
                    dfang = float(inlst_line1[i].split('\t')[4].split('=')[1])
                    for j in range(1,len(inlst_line2)):
                         nn2 = inlst_line2[j].split('\t')[0]
                         mic = inlst_line2[j].split('\t')[1]
                         if(nn == nn2 and mic == inlst_line1[i].split('\t')[1]):
                              ox = float(inlst_line2[j].split('\t')[5].split('=')[1].split(',')[0])
                              oy = float(inlst_line2[j].split('\t')[5].split('=')[1].split(',')[1])
                              break
                    euler1 = float(inlst_line1[i].split('\t')[5].split('=')[1].split(',')[0])
                    euler2 = float(inlst_line1[i].split('\t')[5].split('=')[1].split(',')[1])
                    euler3 = float(inlst_line1[i].split('\t')[5].split('=')[1].split(',')[2])
                    score = float(inlst_line1[i].split('\t')[7].split('=')[1])
             #   m[i-3] = int(inlst_line1[i].split('\t')[0])
                    s5 = "euler="
                    s6 = "center="
                    cx = (float(inlst_line1[i].split('\t')[6].split('=')[1].split(',')[0]))*scale+(ox-clip*scale/2)
                    cy = (float(inlst_line1[i].split('\t')[6].split('=')[1].split(',')[1]))*scale+(oy-clip*scale/2)
                    oo.write(str(micro)+"\t"+str(cx)+"\t"+str(cy)+"\t000001@"+str(micro)+"\t"+str(dfu)+"\t"+str(dfv)+"\t"+str(dfang)+"\t300.000000\t2.700000\t0.100000\t"+str(mag)+"\t5.000000\t0.140000\t1\t"+str(math.fmod(euler2-90,360))+"\t"+str(euler1)+"\t"+str(math.fmod(euler3+90,360))+"\t0.000000\t0.000000\t1\t0.600000\t1.800000e+05\t1.000000\t1"+"\n")                           
        g.close()
        f.close()
	oo.close()

	
def parse_command_line():
	usage="%prog <Merged lstfile> <original particle file> <scale factor> <scaled window size> <unbinned pixel size> <out star file>"
	parser = OptionParser(usage=usage, version="%1")
	
	if len(sys.argv)<7: 
		print "<Merged lstfile> <original particle file> <scale factor> <scaled window size> <unbinned pixel size> <out star file>"
		sys.exit(-1)
	
	(options, args)=parser.parse_args()
	
	inlst1 = args[0]
	inlst2 = args[1]
        scale = int(args[2])
        clip = int(args[3])
        apix = float(args[4])
        output = args[5]
	return (inlst1, inlst2, scale, clip, apix, output)

if __name__== "__main__":
	main()
