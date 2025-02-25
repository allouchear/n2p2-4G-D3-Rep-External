#!/bin/bash

fileIn=$1
fileOut=$2

fmopacP=${fileIn%.inp}_mop
fmopacIn=${fmopacP}.mop
fmopacOut=${fmopacP}.out

nThreads=1
a0=0.52917721
convkcalHartree=1.59360150e-03
method="AM1"
export LC_NUMERIC="C"

#cat $fileIn
# Build MOPAC input file
awk -v nThreads=$nThreads -v a0=$a0 -v method=$method 'BEGIN{
n=0
type=0
nTv=0
}
{
	if(NF>0)
	{
		if(n==0)
		{
			type=$1
			charge=$2
			multiplicity=$3
			ms=(multiplicity-1)/2
			uhf=" "
			if(multiplicity%2==0) uhf=" UHF "
			if(type==0)
			{
				printf("%s %s 1SCF AUX GNORM=0.01 CHARGE=%s MS=%f THREADS=%d\n", method, uhf, charge, ms, nThreads)
				printf("\n");
				printf("Mopac file generated by nnp-script-pm7.bash \n");
			}
			else
			{
				printf("%s %s 1SCF GRADIENTS AUX GNORM=0.01 CHARGE=%d MS=%f THREADS=%d\n",method, uhf, charge, ms, nThreads)
				printf("\n");
				printf("Mopac file generated by nnp-script-pm7.bash \n");
			}
		}
		else if(NF==4)
		{
			if($1 ~ /Tv/ && NF==4)
			{
				Tv[nTv] = sprintf("%s %20.10e %20.10e %20.10e",$1,$2*a0,$3*a0,$4*a0);
				nTv += 1
			}
			else
			{
				printf("%s %20.10e %20.10e %20.10e\n",$1,$2*a0,$3*a0,$4*a0);
			}
		}
	}
	n +=1;
	
	
}
END{
for(n=0;n<nTv;n++)
	printf("%s\n",Tv[n])
	
}' $fileIn  > $fmopacIn

# run MOPAC
mopac22 $fmopacIn >& /dev/null

# Get energy  from MOPAC output file
#  FINAL HEAT OF FORMATION =      -4657.36884 KCAL/MOL =  -19486.43122 KJ/MOL
energy=$(grep "FINAL HEAT OF FORMATION"  $fmopacOut | awk -v kh=$convkcalHartree '{printf("%20.10e",$6*kh)}')
echo $energy > $fileOut 
type=$(grep -c GRADIENTS $fmopacIn)
if [ $type -gt 0 ] 
then
	# Get energy from MOPAC output file
	grep " CARTESIAN " $fmopacOut | grep ANGSTROM  | awk -v a0=$a0 -v kh=$convkcalHartree 'BEGIN{
		kha0=kh*a0
	}
	{
		printf("%20.10e ", -$7*kha0) # sign - : forces = - gradients
		if($0 ~ / CARTESIAN Z/) printf("\n")
	}' >> $fileOut

fi
#rm ${fmopacP}.*


