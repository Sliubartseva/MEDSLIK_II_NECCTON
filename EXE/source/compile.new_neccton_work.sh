# load modules User must deploy own environment
#module load intel19.5/19.5.281
#module load intel19.5/netcdf/C_4.7.2-F_4.5.2_CXX_4.3.1
#module load intel19.5/hdf5/1.10.5
#module load intel19.5/szip/2.1.1
#module load curl/7.66.0

# set folders
DIR_EXE=$MEDSLIK_EXE
DIR_SRC=$MEDSLIK_EXE/MODEL_SRC
NETCDF=/zeus/opt/intel19.5/netcdf

# 2021.10.16
DIR_SRC=.
DIR_EXE=.

# THe paths below must be modified by user
#NETCDF=/zeus/opt/intel19.5/netcdf/C_4.7.2-F_4.5.2_CXX_4.3.1
#MEDSLIK_DIR=/work/opa/sl35719/MEDSLIK_II_1.01/neccton
#MEDSLIK_FOR=$MEDSLIK_DIR/source/medslik_II_neccton_work.for
#MEDSLIK_EXE=$MEDSLIK_DIR/medslik_II_neccton_work.exe

# compile medslik without checking the bounds
gfortran -I$NETCDF/include -L$NETCDF/lib  $MEDSLIK_FOR -lnetcdf -lnetcdff -o $MEDSLIK_EXE

# compile medslik with checking the bounds
#gfortran -fcheck=bounds -I$NETCDF/C_4.7.2-F_4.5.2_CXX_4.3.1/include -L$NETCDF/C_4.7.2-F_4.5.2_CXX_4.3.1/lib  $DIR_SRC/medslik_II_neccton_test.for -lnetcdf -lnetcdff -o $DIR_EXE/medslik_II_neccton_test.exe

ifort -O3 -o $DIR_EXE/Extract_II_neccton_cmems24_currents.exe $DIR_SRC/Extract_II_neccton_cmems24_currents.for -I$NETCDF/C_4.7.2-F_4.5.2_CXX_4.3.1/include -L/$NETCDF/C_4.7.2-F_4.5.2_CXX_4.3.1/lib -lnetcdff -lnetcdf -fPIC -shared-intel -mcmodel=large
ifort -O3 -o $DIR_EXE/Extract_II_neccton_cmems24_wind.exe $DIR_SRC/Extract_II_neccton_cmems24_wind.for -I$NETCDF/C_4.7.2-F_4.5.2_CXX_4.3.1/include -L/$NETCDF/C_4.7.2-F_4.5.2_CXX_4.3.1/lib -lnetcdff -lnetcdf -fPIC -shared-intel -mcmodel=large

cp $DIR_EXE/Extract_II_neccton_cmems24_currents.exe ../$DIR_EXE/Extract_II_neccton_cmems24_currents.exe
cp $DIR_EXE/Extract_II_neccton_cmems24_wind.exe ../$DIR_EXE/Extract_II_neccton_cmems24_wind.exe

# compile jday.f
gfortran -o $DIR_EXE/jday $DIR_SRC/jday.f 

# compile lat_lon
gfortran -o $DIR_EXE/lat_lon.exe $DIR_SRC/lat_lon.for
