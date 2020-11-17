#! /bin/bash

#iteration number
X=1
ladj_ems=1

if [ ladj_ems ]; then
   rm ./BEC_wn13/FLAG_* 
   rm ./BEC_wn13/BINVX_*
fi
cp  GEOSRestart/initial_GEOSChem_rst.2x25_CO2.20171001.nc GEOSRestart/initial_GEOSChem_rst.2x25_CO2.nc
cp FdCo2Out/GEOSChem.SpeciesConc.20171001_0000z.bac.nc4 FdCo2Out/GEOSChem.SpeciesConc.20171001_0000z.nc4
cp  GEOSRestart/initial_GEOSChem_rst.2x25_CO2.nc GEOSRestart/initial_GEOSChem_rst.2x25_CO2_background.nc
MAXITE=25
min_var=0.001
while [ 1 ]
do
  echo "$X" > ITER
  ./geos #&> gcadj.log.$(printf "%02d" $X)
  if [ ladj_ems ]; then
  cp ./BEC_wn13/FLAG_* ./BEC/temp/
  cp ./BEC_wn13/BINVX_* ./BEC/temp/
  rm ./BEC_wn13/FLAG_*
  rm ./BEC_wn13/BINVX_*
  fi
#  if [ $X -eq 1 ];then
#  mv FdCo2Out/GEOSChem.Species* GEOSChemSpe/
#  mv GEOSChemSpe/GEOSChem.SpeciesConc.20160101_0000z.nc4 FdCo2Out/GEOSChem.SpeciesConc.20160101_0000z.nc4
#  fi
#  cp GEOSRestart/initial_GEOSChem_rst.2x25_CO2.nc GEOSRestart/initial_GEOSChem_rst.2x25_CO2.$X.nc
  if [ `cat exit.flag`  -eq 1 ]; then
     break
  fi
#  cost_fir=`printf %f $(cat cfn.$(printf "%02d" 1))`
#  cost_cur=`printf "%f" $(cat cfn.$(printf "%02d" $X)) `
  tot_ite=$(printf "%04d" `cat TOT_ITER`)
#  cost_ite_exit=`echo "scale=10; $cost_fir*1.0" | bc`
 if [ $X -eq 1 ]; then
    cost_ite_diff=100.0
 else
    lastX=$((X-1))
    cost_last=`printf "%f" $(cat cfn.$(printf "%02d" $lastX)) `
    cost_cur=`printf "%f" $(cat cfn.$(printf "%02d" $X)) `
    cost_ite_diff=`echo "scale=10; $cost_last-$cost_cur" | bc`
 fi

 if [ $(echo "$cost_ite_diff <= 0."|bc) -eq 1 ]; then
   cost_ite_diff=$(echo "-1*$cost_ite_diff"|bc)
 fi

  #copy cost function 
  cp cfn.$(printf "%02d" $X) cost_function/cfn.${tot_ite}.$(printf "%02d" $X)
  cp cfnjo.$(printf "%02d" $X) cost_function/cfnjo.${tot_ite}.$(printf "%02d" $X)

#  if [ $(echo "$cost_cur <= $cost_ite_exit"|bc) -eq 1 -a $X -ge 20 ];then
#      echo "nomal"
#      break
#  fi
 
  if [ $(echo "$cost_ite_diff <= $min_var"|bc) -eq 1 ];then
      echo "nomal"
      break
  fi

  if [ $X -ge $MAXITE ];then
     break
  fi
#  if [ $X -eq 4 ]; then
#     break
#  fi
X=$((X+1))
done

