ROOT             := ..
LIB              := $(ROOT)/lib
MOD              := $(ROOT)/mod
SOURCES          := $(wildcard *.f90) 
TEMP             := $(SOURCES:.f90=.o)
OBJECTS          := $(TEMP:.f90=.o)
MODULES          := $(OBJECTS:.o=.mod)
MODULES          := $(foreach I, $(MODULES), $(MOD)/$(I))
LIBRARY          := libADJOINT.a
include $(ROOT)/Makefile_header.mk
$(warning $(OBJECTS))
all: lib

lib: $(OBJECTS)
	$(AR) crs $(LIBRARY) $(OBJECTS)
	mv $(LIBRARY) $(LIB)

clean: 
	rm -f *.o $(MODULES)

logical_adj_mod.o        :  logical_adj_mod.f90

errcode_adj_mod.o        :  errcode_adj_mod.f90

parameter_adj_mod.o      :  parameter_adj_mod.f90

input_adj_mod.o          :  input_adj_mod.f90

file_ope_mod.o           :  file_ope_mod.f90                                   \
                            inverse_mod.o             state_chm_adj_mod.o      \
                            state_grid_adj_mod.o      state_met_adj_mod.o
             
adj_arrays_mod.o         :  adj_arrays_mod.f90                                 \
                            state_chm_adj_mod.o       state_grid_adj_mod.o     \
                            state_met_adj_mod.o
tpcore_fvdas_adj_mod.o   :  tpcore_fvdas_adj_mod.f90

tpcore_fvdas_adj_mod_d.o :  tpcore_fvdas_adj_mod_d.f90

state_chm_adj_mod.o      :  state_chm_adj_mod.f90                              \
                            errcode_adj_mod.o         state_grid_adj_mod.o

state_grid_adj_mod.o     :  state_grid_adj_mod.f90                             \
                            errcode_adj_mod.o

state_met_adj_mod.o      :  state_met_adj_mod.f90                              \
                            errcode_adj_mod.o         state_grid_adj_mod.o  
 
checkpt_mod.o            :  checkpt_mod.f90                                    \
                            state_grid_adj_mod.o      state_chm_adj_mod.o      \
                            state_met_adj_mod.o

geos_chem_adj_mod.o      :  geos_chem_adj_mod.f90                              \
                            state_grid_adj_mod.o      state_chm_adj_mod.o      \
                            checkpt_mod.o             logical_adj_mod.o        \
                            oco2_co2_mod.o            file_ope_mod.o           \
                            pbl_mix_adj_mod.o         convection_adj_mod.o     \
                            input_adj_mod.o           tpcore_fvdas_adj_mod.o   \
                            transport_adj_mod.o       adj_arrays_mod.o         \
                            tansat_co2_mod.o          tccon_co2_mod.o          \
                            co2_adj_mod.o 

co2_adj_mod.o            :  co2_adj_mod.f90                                    \
                            state_chm_adj_mod.o       state_grid_adj_mod.o     \
                            state_met_adj_mod.o       logical_adj_mod.o        \
                            file_ope_mod.o
                     
convection_adj_mod.o     :  convection_adj_mod.f90                             \
                            state_chm_adj_mod.o       state_grid_adj_mod.o     \
                            state_met_adj_mod.o       errcode_adj_mod.o          

pbl_mix_adj_mod.o        :  pbl_mix_adj_mod.f90                                \
                            state_chm_adj_mod.o       state_grid_adj_mod.o     \
                            state_met_adj_mod.o       parameter_adj_mod.o      \
                            adj_arrays_mod.o 

transport_adj_mod.o      :  transport_adj_mod.f90                              \
                            tpcore_fvdas_adj_mod.o    state_chm_adj_mod.o      \
                            state_grid_adj_mod.o      state_met_adj_mod.o      \
                            parameter_adj_mod.o       input_adj_mod.o          \
                            errcode_adj_mod.o         tpcore_fvdas_adj_mod.o   \
                            tpcore_fvdas_adj_mod_d.o

inverse_mod.o            :  inverse_mod.f90                                    \
                            adj_arrays_mod.o          logical_adj_mod.o        \
                            state_chm_adj_mod.o       state_grid_adj_mod.o     \
                            parameter_adj_mod.o 

oco2_co2_mod.o           :  oco2_co2_mod.f90                                   \
                            state_chm_adj_mod.o       state_grid_adj_mod.o     \
                            state_met_adj_mod.o       adj_arrays_mod.o

tansat_co2_mod.o         :  tansat_co2_mod.f90                                 \
                            state_chm_adj_mod.o       state_grid_adj_mod.o     \
                            state_met_adj_mod.o       adj_arrays_mod.o

tccon_co2_mod.o          :  tccon_co2_mod.f90                                  \
                            state_chm_adj_mod.o       state_grid_adj_mod.o     \
                            state_met_adj_mod.o       adj_arrays_mod.o

geos_chem_d_mod.o        :  geos_chem_d_mod.f90                                \
                            state_grid_adj_mod.o      state_chm_adj_mod.o      \
                            checkpt_mod.o             logical_adj_mod.o        \
                            oco2_co2_mod.o            file_ope_mod.o           \
                            pbl_mix_adj_mod.o         convection_adj_mod.o     \
                            input_adj_mod.o           tpcore_fvdas_adj_mod_d.o \
                            transport_adj_mod.o       adj_arrays_mod.o         \
                            tansat_co2_mod.o          tccon_co2_mod.o

assimilation_structure.o :  assimilation_structure.f90                         \
                            state_chm_adj_mod.o       state_chm_adj_mod.o      \
                            state_met_adj_mod.o       input_adj_mod.o          \
                            logical_adj_mod.o         geos_chem_adj_mod.o      \
                            geos_chem_d_mod.o         tansat_co2_mod.o         \
                            background_error_covariance.o
