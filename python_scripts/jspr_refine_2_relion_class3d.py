#!/usr/bin/env python

import math, os, sys
try:
	from optparse import OptionParser
except:
	from optik import OptionParser

def main():
	(lstfile,ny,apix,output, meta,debug) =  parse_command_line()
	g = open(lstfile, "r")
	lst_line = g.readlines()
	r=open(output,"w")
        line_len=int(len(lst_line[4].split()))
        for j in range(2,line_len):
            tmp=str(lst_line[4].split('\t')[j].split('=')[0])
            if(tmp=="defocus"):
	        df_index=j
            elif(tmp=="dfdiff"):
		dfdiff_index=j
	    elif(tmp=="dfang"):
		dfang_index=j
	    elif(tmp=="euler"):
		euler_index=j
	    elif(tmp=="center"):
		center_index=j
        r.write("\ndata_")
        r.write("\n\n")
        r.write("loop_\n")
        if(debug==0):
            r.write("_rlnMicrographName #1\n_rlnCoordinateX #2\n_rlnCoordinateY #3\n_rlnImageName #4\n_rlnDefocusU #5\n_rlnDefocusV #6\n_rlnDefocusAngle #7\n_rlnVoltage #8\n_rlnSphericalAberration #9\n_rlnAmplitudeContrast #10\n_rlnMagnification #11\n_rlnDetectorPixelSize #12\n_rlnCtfFigureOfMerit #13\n_rlnGroupNumber #14\n_rlnAngleRot #15\n_rlnAngleTilt #16\n_rlnAnglePsi #17\n_rlnOriginX #18\n_rlnOriginY #19\n_rlnClassNumber #20\n_rlnNormCorrection #21\n_rlnLogLikeliContribution #22\n_rlnMaxValueProbDistribution #23\n_rlnNrOfSignificantSamples #24\n")
        else:
            r.write("_rlnMicrographName #1\n_rlnCoordinateX #2\n_rlnCoordinateY #3\n_rlnImageName #4\n_rlnDefocusU #5\n_rlnDefocusV #6\n_rlnDefocusAngle #7\n_rlnVoltage #8\n_rlnSphericalAberration #9\n_rlnAmplitudeContrast #10\n_rlnMagnification #11\n_rlnDetectorPixelSize #12\n_rlnCtfFigureOfMerit #13\n_rlnGroupNumber #14\n_rlnAngleRot #15\n_rlnAngleTilt #16\n_rlnAnglePsi #17\n_rlnOriginX #18\n_rlnOriginY #19\n_rlnClassNumber #20\n_rlnNormCorrection #21\n_rlnLogLikeliContribution #22\n_rlnMaxValueProbDistribution #23\n_rlnNrOfSignificantSamples #24\n_rlnRandomSubset #25\n")
	for i in range (meta,len(lst_line)):
                num=int(lst_line[i].split()[0])+1
                filename=lst_line[i].split()[1]
		defocus=float((lst_line[i].split()[df_index]).split('=')[1])*10000
		dfdiff=float((lst_line[i].split()[dfdiff_index]).split('=')[1])*10000
		dfang=float((lst_line[i].split()[dfang_index]).split('=')[1])
		dfu=defocus-dfdiff
		dfv=defocus+dfdiff
		euler_tmp=(lst_line[i].split()[euler_index]).split('=')[1]
		euler1=float(euler_tmp.split(',')[0])
		euler2=float(euler_tmp.split(',')[1])
		euler3=float(euler_tmp.split(',')[2])
		center_tmp=(lst_line[i].split()[center_index]).split('=')[1]
		center1=math.floor(float(center_tmp.split(',')[0])+0.5)
		center1_pointer=ny/2-float(center_tmp.split(',')[0])
		center2=math.floor(float(center_tmp.split(',')[1])+0.5)
		center2_pointer=ny/2-float(center_tmp.split(',')[1])
		apix_c=apix*0.0001
		apix_c=5/apix_c
		r.write(filename+"\t")
		r.write(str(center1)+"\t"+str(center2)+"\t")
		r.write("%06d" % int(num))
		r.write("@"+filename+"\t")
		r.write(str(dfu)+"\t"+str(dfv)+"\t"+str(dfang)+"\t300.000000\t2.7000000\t0.100000\t"+str(apix_c)+"\t5.000000\t")
		r.write("0.14000\t1\t")
#		r.write(str(euler1-90)+"\t"+str(math.fmod(euler2+180,360))+"\t"+str(math.fmod(euler3+90,360))+"\t")
		r.write(str(math.fmod(euler2-90,360))+"\t"+str(euler1)+"\t"+str(math.fmod(euler3+90,360))+"\t")
	# 	only for EMAN ICOS to relion I3 symmetry
		r.write(str(center1_pointer)+"\t"+str(center2_pointer)+"\t")
                if(debug==0):
		    r.write("1\t0.6\t1.8e+05\t1.000000\t1\n")
                if(debug==1):
                    r.write("1\t0.6\t1.8e+05\t1.000000\t1\t1\n")
                if(debug==2):
                    r.write("1\t0.6\t1.8e+05\t1.000000\t1\t2\n")
	
	g.close()
	r.close()
	
def parse_command_line():
	usage="%prog <lstfile> <ny> <apix> <output> < metadata line> <1,odd 2,even 0,none>"
	parser = OptionParser(usage=usage, version="%1")
	
	if len(sys.argv)<5: 
		print "<lstfile> <ny> <apix> <output> <metadata line> <1,odd 2,even 0,none>"
		sys.exit(-1)
	
	(options, args)=parser.parse_args()
	lstfile = args[0]
	ny=float(args[1])
	apix=float(args[2])
	output=args[3]
        meta=int(args[4])
        debug=int(args[5])
	return (lstfile,ny,apix,output,meta,debug)

if __name__== "__main__":
	main()


			
