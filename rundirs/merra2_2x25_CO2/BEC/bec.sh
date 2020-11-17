#!/bin/bash

#delete log
rm gcadj.log*

gfortran -fopenmp bec.f90 -o bec
true=1
false=0
lfirst=$true

startymdHM="20171001 0000"
endymdHM="20171001 0600"

startymdHMts=`date -d "$startymdHM" +%s`
endymdHMts=`date -d "$endymdHM" +%s`

echo $startymdHMts
echo $endymdHMts

#cp  GEOSRestart/initial_GEOSChem_rst.2x25_CO2.20171001.nc GEOSRestart/initial_GEOSChem_rst.2x25_CO2.nc
#cp  GEOSRestart/initial_GEOSChem_rst.2x25_CO2.20171001.nc GEOSRestart/initial_GEOSChem_rst.2x25_CO2_background.nc
#cp FdCo2Out/GEOSChem.SpeciesConc.20171001_0000z.bac.nc4 FdCo2Out/GEOSChem.SpeciesConc.20171001_0000z.nc4

MAXITE=30
itenum=1

tot_ite=1
lics=0
bec=1
#wrf4dvar=0
#python_bin='/home/huoxiao/application/anaconda3/envs/esmpy/bin/python'
while [ 1 ]
do
      echo enter
      exit_flag=`(cat EXIT_FLAG)`
      if [ $exit_flag -eq 1 ]; then
         echo exit
         exit

      fi
      while [ 1 ] 
      do
         sf_success=`(cat SF_SUCCESS)`
         if [ $sf_success -eq 1 ]; then
            break 1
         fi
      done
      scp gctm.gc.01 wn12:~/GEOS_CHEM_ADJOINT_V1_INC/rundirs/merra2_2x25_CO2/BEC/
      scp gctm.gc.01 wn13:~/GEOS_CHEM_ADJOINT_V1_INC/rundirs/merra2_2x25_CO2/BEC/
      scp gctm.gc.01 wn14:~/GEOS_CHEM_ADJOINT_V1_INC/rundirs/merra2_2x25_CO2/BEC/
      scp gctm.gc.01 wn12:~/GEOS_CHEM_ADJOINT_V1_INC/rundirs/merra2_2x25_CO2/BEC_wn11/
      scp gctm.gc.cur wn12:~/GEOS_CHEM_ADJOINT_V1_INC/rundirs/merra2_2x25_CO2/BEC/
      scp gctm.gc.cur wn13:~/GEOS_CHEM_ADJOINT_V1_INC/rundirs/merra2_2x25_CO2/BEC/
      scp gctm.gc.cur wn14:~/GEOS_CHEM_ADJOINT_V1_INC/rundirs/merra2_2x25_CO2/BEC/
      scp gctm.gc.cur wn12:~/GEOS_CHEM_ADJOINT_V1_INC/rundirs/merra2_2x25_CO2/BEC_wn11/

      echo 1 > master_f
      scp master_f wn12:~/GEOS_CHEM_ADJOINT_V1_INC/rundirs/merra2_2x25_CO2/BEC/
      scp master_f wn13:~/GEOS_CHEM_ADJOINT_V1_INC/rundirs/merra2_2x25_CO2/BEC/
      scp master_f wn14:~/GEOS_CHEM_ADJOINT_V1_INC/rundirs/merra2_2x25_CO2/BEC/
      scp master_f wn12:~/GEOS_CHEM_ADJOINT_V1_INC/rundirs/merra2_2x25_CO2/BEC_wn11/
      ./bec
#      cp gctm.gc.01 temp/
      
#      rm gctm.gc.01
#      rm gctm.gc.cur
      echo 0 > SF_SUCCESS
#      scp BINVX_01 wn15:~/GEOS_CHEM_ADJOINT_V1_INC/rundirs/merra2_2x25_CO2/BEC/
#      scp FLAG_01 wn15:~/GEOS_CHEM_ADJOINT_V1_INC/rundirs/merra2_2x25_CO2/BEC/
done 
