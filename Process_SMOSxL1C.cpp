// "Process_SMOSxL1C": mex function to process L1C SMOS data with Octave/Matlab
//
// Copyright (c) 2016 Pablo Saavedra Garfias
//
// this file is part of 'Process_SMOSxL1C'
//
// 'Process_SMOSxL1C' is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// 'Process_SMOSxL1C' is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with 'Process_SMOSxL1C'. If not, see <http://www.gnu.org/licenses/>.

// ---------------------------------------------------------------------------
// Program to read SMOS SCL1C DBL data files (either dual or full pol)

// version 1.0
// COMPILATION:
// (from linux console):
// >[MATLAB_BIN_PATH]/mex CFLAGS='$CFLAGS -std=gnu++11' readSCL1C_smos_DBL.cpp -lgsl -lgslcblas
// (from MATLAB workspace):
// > mex Process_SMOSxL1C.cpp -lgsl -lgslcblas
//
// USAGE 1:
// > [TSF, SSI] = Process_SMOSxL1C;
// then, a file browser will be open and select a DBL file to open.
// USAGE 2:
// > [TSF, SSI] = Process_SMOSxL1C('fname.DBL',GEO_LIMS,outdir);
// WHERE:
// 'fname.DBL': full path string of DBL file to work with,
// GEOLIMS: a 4 element vector as [LAT1, LON1, LAT2, LON2],
// with elements indicating the latitude and longitude of bottom left
// of box and latitude and longitude of upper rigth of box.
// outdir: (optional) string with directory where filename.mat is stored.
// USAGE 3:
// > [TSF, SSI] = Process_SMOSxL1C('fname.DBL',GEO_LIMS);
// same as usage #2, but an output .mat file will not be generated.
// OUTPUTS:
// TSF: Structure with all variables for [Temp_Swath_Full] dataset.
// SSI: Structure with all variables for [Snapshot_Information]
// dataset (this output variable is optional).
// * Additionally, TSF includes the Brightness Temperatue in HV-pol after
// performing the Faraday rotation from the XY polarization plane.
// * The TSF memeber variables TB_Fixed_IncAngle and Fixed_IncAngle,
// are reprecenting the Brightness Temperatures at homogeneously spaced
// fixed incident angles (normally from 0 to 65 deg with 1 deg step).

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cstring>
#include <cmath>
#include <time.h>
#define __STDC_LIMIT_MACROS
#include <stdint.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#ifdef MATLAB_MEX_FILE
// NOTE: The following three lines-block  [#define char16_t] to [#undef char16_t]
// is a work-around for a problem during compilation error with c++11:
// redeclaration of C++ built-in type ‘char16_t’,  typedef CHAR16_T char16_t;
#define char16_t LIBRARY_char16_t
#include <mex.h>
#undef char16_t      // end of work-around char16_t.
#include "mat.h"
#define MYNAN mxGetNaN()
#else
#include <mex.h>
#include "matio.h"
#define MYNAN NAN
#endif

// -- Definition of global parameters --
#define PI 4.0*atan(1.0)
#define FAK90DEG 90.0/(UINT16_MAX+1)    //  Faktor: 90/2^-16, hint. UINT16_MAX=2^16-1
#define FAK360DEG 360.0/(UINT16_MAX+1)  // Faktor: 360 [deg] for Faraday and Geometric angle
#define FAKRAD 2*PI/(UINT16_MAX+1)      // Faktor: 2PI [rad] for Faraday and Geometric angle
#define TH_SIZE 80      // [km] Spatial resolution size threshold (SMOS specification defult 55 km)
#define TH_ELON 2.2     // Spatial resolution elongation threshold (specification default 1.5)
#define TBXX_MIN 50     // Minimum threshold for brigthness temperature [K]
#define TBXX_MAX 500    // Maximum threshold for brigthness temperature [K]
#define TH_RFI_ST4 50   // Threshould for 4-Stoke parameter Radio Frequency Interference
#define CA_TBS1 5       // parameter A for equation
#define CB_TBS1 4       // parameter B for equation ST4-<ST4> = A + B
#define INCA_INI 0      // initial incidence angle [deg]
#define INCA_NUM 61     // number of incidence angles [deg] (e.g. 66 SMOS)
#define INCA_DEL 1      // steps for  incidence angle [deg]
// -- end of definition for global parameters

using namespace std;

namespace smos{


  // ****** CLASS DEFINITION *********************
  class Be_Discrete {
  private:
    int k, Ninc, Init, Step;
    double *TB_i, *TBave;
    int *N_i;
    double Delta;

  public:
    double *theta;
    // default creator for Be_Discrete class:
    Be_Discrete () : Ninc(INCA_NUM),
      Init(INCA_INI), Step(INCA_DEL) {
      initialize_it();
    }
    // creator for Be_Discrete class with parameters:
    Be_Discrete (int I, int S, int N) : Ninc(N),
      Init(I), Step(S) {
      initialize_it();
    }
    // destroyer for Be_Discrete:
    ~Be_Discrete () {
      delete [] N_i; delete [] TBave; delete [] TB_i; delete [] theta;
    }
    // return the average Brightness temperature within the inc_angle interval:
    double *deliver_it() {
      for(k=0;k<Ninc;++k) TBave[k]=N_i[k]!=0?TB_i[k]/N_i[k]:MYNAN;
      return(TBave);
    }
    // shows the results of the discretization:
    void show_it(){
      for(k=0;k<Ninc;++k)
	cout<<theta[k]<<" "<<TBave[k]<<";"<<endl;
    }
    // take a incidence angle and its TB to include in the discretization:
    void consider_it(double beta,double VARin){
      for(k=0; k<Ninc; ++k)
	if(beta>=(theta[k]-Delta) && beta<(theta[k]+Delta) && VARin==VARin){
	  // NOTE here that VARin==VARin is FALSE only when VARin=NaN
	  TB_i[k] += VARin;
	  N_i[k]++;
	}
      return;
    }

  protected:
    void initialize_it(){
      Delta = double(Step)/2;
      N_i   = new int[Ninc];
      TBave = new double[Ninc];
      TB_i  = new double[Ninc];
      theta = new double[Ninc];
      for(k=0;k<Ninc;++k){
	TBave[k]=0; TB_i[k]=0; N_i[k]=0;
	theta[k]= (double) k*Step+Init;
      };
    }
  };

  // ****** SUBROUTINES DECLARATION **************
  int GetInputDBL_BoxLIM(char *&, char *&, float []);
  void MxTrans_XY2HV(double, double [], double [], int, bool);
  void interp1(double [], double [], int, double [], double [], int);
  float getHDR_fields(const char [], const char []);
  int CreateMAToutput(mxArray *, mxArray *, char [], char []);
#ifndef MATLAB_MEX_FILE
  matvar_t *CloneSTRUCTforMATIO(mxArray *, const char *);
#endif
  int TimeEEC2date(int [], char []);
  void TBHV_at_FixTheta(double [], double [], int, double []);
  void ShowUsage(const char*);
  void ShowGNUPL();

  bool SSI_FLAG=true;
} // end of smos namespace


// ::::::::::: MAIN PROGRAM STARTS ::::::::::::::::::
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

  // Defining the output structure variables:
  mxArray *SnapShotInfo, *TempSwathFull;

  // Following mex variables for the pixel variables:
  const char *field_name[] = {"File_Name","GEOBoxLimits","Start_End_Time",
			      "GridPoint_ID","GridPoint_Latitude","GridPoint_Longitude",
			      "GridPoint_Altitude","GridPoint_Mask","BT_Data_Counter",
			      "flags","BTvalue_re","BTvalue_im","PixelRadiometry_acc",
			      "Incidence_angle","Azimuth_angle","FaradayRotation_angle",
			      "GeometricRotation_angle","Snapshot_ID_Pixels",
			      "Footprint_Axis1","Footprint_Axis2","TBxy","TBhv",
			      "RA","idx_SnapshotID","TB_Fixed_IncAngle","Fixed_IncAngle"};

  const unsigned N_field_name = sizeof(field_name)/sizeof(field_name[0]);

  // Following mex variables for the snapshot variables:
  const char *snapsh_name[] = {"Snapshot_Time","Snapshot_ID","Snapshot_OBET",
			       "XYZ_Position","XYZ_Velocity","Vector_So","Qnn","TEC",
			       "Geomag_F_D_I","Sun_RA_DEC_BT",
			       "Accuracy","Radiometric_Accuracy","XBand","SIAC_Error_flag"};
  const unsigned N_snapsh_name = sizeof(snapsh_name)/sizeof(snapsh_name[0]);


  char *filen, *OUTDIR;
  ifstream fp;
  unsigned int i,j, inn, Snapshot_Counter;
  unsigned int Grid_Point_Counter;
  double RAD_ACC_SC, PIX_FOOTP_SC;
  unsigned STOKES_VEC = 4;      // default: 4 for Stokes vector (optional: 2 for dual-pol)

  // List of Grid Point Data (TempSwathFull variables):
  unsigned int Grid_Point_ID;
  float Grid_Point_Latitude, Grid_Point_Longitude, Grid_Point_Altitude;
  unsigned char Grid_Point_Mask;
  unsigned short int BT_Data_Counter;
  unsigned short int Flags, Pixel_Radiometric_Accuracy;
  float BT_Value_Real, BT_Value_Imag;
  unsigned short int Incidence_Angle, Azimuth_Angle;
  unsigned short int Faraday_Rotation_Angle, Geometric_Rotation_Angle;
  unsigned int Snapshot_ID_of_Pixel;
  unsigned short int Footprint_Axis1, Footprint_Axis2;
  // Others:
  bool INMYBOX; //bool smos::SSI_FLAG=true;
  float BOXLIM[4] = {-90,-180,90,180};   // Default Latitude-Longitude Box
  char smos_date0[]="yyyymmddThh:mi:ss", smos_date1[]="yyyymmddThh:mi:ss";
  // ---------------------------
  // Checking output parameters:
  if(nlhs<1 || nlhs>2) smos::ShowUsage("Need at least one output variable, try again ;)");
  if(nlhs==1) smos::SSI_FLAG=false; else smos::SSI_FLAG=true;

  // Checking the input parameters:
  if (nrhs>=0 && nrhs<4){
    int FileLength;

    switch(nrhs){
    case 3:
      // third input: Output file path
      if (mxIsChar(prhs[2])){
	FileLength = mxGetN(prhs[2])+1;
	OUTDIR = (char *) mxCalloc(FileLength, sizeof(char));
	mxGetString(prhs[2], OUTDIR, FileLength);
	cout<<"% Output dir: "<<OUTDIR<<endl;
      }
      else mexErrMsgTxt("Third input needs to be a string PATH.");

    case 2:
      // Second input: Box to extract [lat_min lon_min lat_max lon_max]
      if (mxIsNumeric(prhs[1]) && mxGetNumberOfElements(prhs[1])==4) {     // Getting the Lat, Lon boundaries:
	// Latitude bottom left
	BOXLIM[0] = (float) *(mxGetPr(prhs[1]));
	// Longitude bottom left
	BOXLIM[1] = (float) *(mxGetPr(prhs[1])+1);
	// Latitude upper rigth
	BOXLIM[2] = (float) *(mxGetPr(prhs[1])+2);
	// Longitude upper rigth
	BOXLIM[3] = (float) *(mxGetPr(prhs[1])+3);
      }
      else mexErrMsgTxt("Second input must be 4x1 numeric, try again ;)");

    case 1:
      // first input: file name:
      if(mxIsChar(prhs[0])){
	FileLength = mxGetN(prhs[0])+1;
	filen = (char *) mxCalloc(FileLength, sizeof(char));
	mxGetString(prhs[0], filen, FileLength);
	mexPrintf("DBL file to open: %s\n",filen);
      }
      else mexErrMsgTxt("First input needs to be a string FILENAME.");
      break;

    case 0:
      //open file brower
      if(smos::GetInputDBL_BoxLIM(filen,OUTDIR,BOXLIM)!=0) mexErrMsgTxt("Wrong input DBL file or (Lat,Lon) min-max limits!");
      break;

    default:
      mexErrMsgTxt("Ups! something is wrong with the input variables!");
    }  // end switch
  }
  else smos::ShowUsage("Wrong number of inputs, try again ;)");

  cout<<"% Position min: ("<<BOXLIM[0]<<","<<BOXLIM[1]<<")"<<endl;
  cout<<"% Position max: ("<<BOXLIM[2]<<","<<BOXLIM[3]<<")"<<endl;
  cout<<__cplusplus<<endl;
  // checking the polarization type from DBL file name:
  string WhatPol = filen;
  int idxPol = WhatPol.find("1C_");
  if(WhatPol.at(idxPol-5)=='_' && WhatPol.at(idxPol-1)=='D') STOKES_VEC = 2; // dual-pol
  else STOKES_VEC = 4;   // 4-Stokes-pol

  // Checking whether file name has extention ".DBL" as SMOS data block binary file:
  if(WhatPol.compare(WhatPol.size()-4,4,".DBL")!=0) mexWarnMsgTxt("Are you sure this is a SMOS 1C DBL binary file?");

  // Opening the BDL file:
  fp.open(filen,ios::binary | ios::in);
  if (!fp.is_open()){
    cout<<filen;
    mexErrMsgTxt(" :DBL file cannot be opened, it might be corrupted!");
    return;
  }

  // Reading the necessary HDR fields:
  RAD_ACC_SC = (double) smos::getHDR_fields(filen,"Radiometric_Accuracy_Scale");
  if(RAD_ACC_SC!=MYNAN) RAD_ACC_SC = RAD_ACC_SC/(UINT16_MAX+1);

  PIX_FOOTP_SC = (double) smos::getHDR_fields(filen,"Pixel_Footprint_Scale");
  if(PIX_FOOTP_SC!=MYNAN) PIX_FOOTP_SC = PIX_FOOTP_SC/(UINT16_MAX+1);

  float CREATOR_VER = smos::getHDR_fields(filen,"Creator_Version");
  if(CREATOR_VER!=MYNAN) CREATOR_VER = 620;  // default value

  // Starting to read the data from DBL:
  fp.read((char *) &Snapshot_Counter,sizeof(Snapshot_Counter));

  int Snapshot_Time[Snapshot_Counter][3]; //= new int*[Snapshot_Counter];
  int Snapshot_ID[Snapshot_Counter];   //= new int[Snapshot_Counter];
  long int Snapshot_OBET;   // = new long int[Snapshot_Counter];
  double X_Position;   // = new double[Snapshot_Counter];
  double Y_Position;   // = new double[Snapshot_Counter];
  double Z_Position;   // = new double[Snapshot_Counter];
  double X_Velocity;   // = new double[Snapshot_Counter];
  double Y_Velocity;   // = new double[Snapshot_Counter];
  double Z_Velocity;   // = new double[Snapshot_Counter];
  double TEC;   //        = new double[Snapshot_Counter];
  double Geomag_F;   //   = new double[Snapshot_Counter];
  double Geomag_D;   //   = new double[Snapshot_Counter];
  double Geomag_I;   //   = new double[Snapshot_Counter];
  float Sun_RA;   //     = new float[Snapshot_Counter];
  float Sun_DEC;   //    = new float[Snapshot_Counter];
  float Sun_BT;   //     = new float[Snapshot_Counter];
  float Accuracy[Snapshot_Counter];   //   = new float[Snapshot_Counter];
  float Radiometric_Accuracy[Snapshot_Counter][2];   // = new float*[Snapshot_Counter];
  unsigned char Software_Error_flag;
  unsigned char Instrument_Error_flag;
  unsigned char ADF_Error_flag;
  unsigned char Calibration_Error_flag[Snapshot_Counter];
  unsigned char Vector_Source, X_Band;
  double Qnn[4];

  // creating the mex variables for snapshot variables:
  mxArray *mex_snaptime=NULL, *mex_snapid=NULL, *mex_snapobet=NULL, *mex_xyzpo=NULL, *mex_xyzve=NULL;
  mxArray *mex_vecsou=NULL, *mex_Qnn=NULL, *mex_tec=NULL, *mex_geomag=NULL, *mex_sunv=NULL;
  mxArray *mex_accu=NULL, *mex_radaccu=NULL, *mex_xband=NULL, *mex_errflag=NULL;
  if (smos::SSI_FLAG){
    mex_snaptime = mxCreateDoubleMatrix(Snapshot_Counter,3,mxREAL);
    mex_snapid   = mxCreateDoubleMatrix(Snapshot_Counter,1,mxREAL);
    mex_snapobet = mxCreateDoubleMatrix(Snapshot_Counter,1,mxREAL);
    mex_xyzpo    = mxCreateDoubleMatrix(Snapshot_Counter,3,mxREAL);
    mex_xyzve    = mxCreateDoubleMatrix(Snapshot_Counter,3,mxREAL);
    mex_vecsou   = mxCreateDoubleMatrix(Snapshot_Counter,1,mxREAL);
    mex_Qnn      = mxCreateDoubleMatrix(4,Snapshot_Counter,mxREAL);
    mex_tec      = mxCreateDoubleMatrix(Snapshot_Counter,1,mxREAL);
    mex_geomag   = mxCreateDoubleMatrix(Snapshot_Counter,3,mxREAL);
    mex_sunv     = mxCreateDoubleMatrix(Snapshot_Counter,3,mxREAL);
    mex_accu     = mxCreateDoubleMatrix(Snapshot_Counter,1,mxREAL);
    mex_radaccu  = mxCreateDoubleMatrix(Snapshot_Counter,2,mxREAL);
    mex_xband    = mxCreateDoubleMatrix(Snapshot_Counter,1,mxREAL);
    mex_errflag  = mxCreateDoubleMatrix(Snapshot_Counter,4,mxREAL);
  }

  for(i=0;i<Snapshot_Counter;++i){
    //Snapshot_Time[i] = new int[3];   // (0:days, 1:secods, 2:microseconds)
    // Radiometric_Accuracy[i] = new float[2];
    fp.read((char *) &Snapshot_Time[i][0], 3*sizeof(int));
    fp.read((char *) &Snapshot_ID[i], sizeof(int));
    fp.read((char *) &Snapshot_OBET, sizeof(long int));
    fp.read((char *) &X_Position, sizeof(double));
    fp.read((char *) &Y_Position, sizeof(double));
    fp.read((char *) &Z_Position, sizeof(double));
    fp.read((char *) &X_Velocity, sizeof(double));
    fp.read((char *) &Y_Velocity, sizeof(double));
    fp.read((char *) &Z_Velocity, sizeof(double));
    fp.read((char *) &Vector_Source,1); // sizeof(Vector_Source));
    fp.read((char *) &Qnn[0], 4*sizeof(double));
    fp.read((char *) &TEC, sizeof(double));
    fp.read((char *) &Geomag_F, sizeof(Geomag_F));
    fp.read((char *) &Geomag_D, sizeof(Geomag_D));
    fp.read((char *) &Geomag_I, sizeof(Geomag_I));
    fp.read((char *) &Sun_RA,  sizeof(Sun_RA));
    fp.read((char *) &Sun_DEC, sizeof(Sun_DEC));
    fp.read((char *) &Sun_BT,  sizeof(Sun_BT));
    fp.read((char *) &Accuracy[i], sizeof(Accuracy[0]));
    fp.read((char *) &Radiometric_Accuracy[i][0], 2*sizeof(float));
    fp.read((char *) &X_Band, 1); //sizeof(X_Band));
    fp.read((char *) &Software_Error_flag, sizeof(Software_Error_flag));
    fp.read((char *) &Instrument_Error_flag, sizeof(Instrument_Error_flag));
    fp.read((char *) &ADF_Error_flag, sizeof(ADF_Error_flag));
    fp.read((char *) &Calibration_Error_flag[i], sizeof(Calibration_Error_flag[0]));

    // filling SnapshotInfo variables to MeX:
    if(smos::SSI_FLAG){
      int k0 = i;
      int k1 = i+Snapshot_Counter;
      int k2 = i+2*Snapshot_Counter;
      int k3 = i+3*Snapshot_Counter;
      // Snapshot Time:
      *(mxGetPr(mex_snaptime)+k0) = Snapshot_Time[i][0];   //days
      *(mxGetPr(mex_snaptime)+k1) = Snapshot_Time[i][1];   //seconds
      *(mxGetPr(mex_snaptime)+k2) = Snapshot_Time[i][2]; //microseconds
      // Snapshot ID:
      *(mxGetPr(mex_snapid)+i) = Snapshot_ID[i];
      *(mxGetPr(mex_snapobet)+i) = Snapshot_OBET;
      // X_Y_Z_Position:
      *(mxGetPr(mex_xyzpo)+k0) = X_Position;
      *(mxGetPr(mex_xyzpo)+k1) = Y_Position;
      *(mxGetPr(mex_xyzpo)+k2) = Z_Position;
      // X_Y_Z_Velocity:
      *(mxGetPr(mex_xyzve)+k0) = X_Velocity;
      *(mxGetPr(mex_xyzve)+k1) = Y_Velocity;
      *(mxGetPr(mex_xyzve)+k2) = Z_Velocity;
      // Vector_Source, Qnn():
      *(mxGetPr(mex_vecsou)+i) = Vector_Source;
      memcpy(mxGetPr(mex_Qnn)+i*4, Qnn, 4*sizeof(double));
      *(mxGetPr(mex_tec)+i) = TEC;
      // Geomag_F_D_I:
      *(mxGetPr(mex_geomag)+k0) = Geomag_F;
      *(mxGetPr(mex_geomag)+k1) = Geomag_D;
      *(mxGetPr(mex_geomag)+k2) = Geomag_I;
      // Sun_Right Ascention, Declination, Brightness Temperature:
      *(mxGetPr(mex_sunv)+k0) = Sun_RA;
      *(mxGetPr(mex_sunv)+k1) = Sun_DEC;
      *(mxGetPr(mex_sunv)+k2) = Sun_BT;
      // Accuracies:
      *(mxGetPr(mex_accu)+i) = Accuracy[i];
      *(mxGetPr(mex_radaccu)+k0) = Radiometric_Accuracy[i][0];
      *(mxGetPr(mex_radaccu)+k1) = Radiometric_Accuracy[i][1];
      // XBand:
      *(mxGetPr(mex_xband)+i) = X_Band;
      // Flags:
      *(mxGetPr(mex_errflag)+k0) = Software_Error_flag;
      *(mxGetPr(mex_errflag)+k1) = Instrument_Error_flag;
      *(mxGetPr(mex_errflag)+k2) = ADF_Error_flag;
      *(mxGetPr(mex_errflag)+k3) = Calibration_Error_flag[i];
    } // end if SSI_FLAG
  }   // end loop over snapshots in orbit.


  // Here start the reading of Snapshots:
  fp.read((char *) &Grid_Point_Counter, sizeof(Grid_Point_Counter));

  // Creating mex variables for TempSwathFull output variable:
  mxArray *mex_gridid = mxCreateDoubleMatrix(Grid_Point_Counter,1,mxREAL);
  mxArray *mex_lat  = mxCreateDoubleMatrix(Grid_Point_Counter,1,mxREAL);
  mxArray *mex_lon  = mxCreateDoubleMatrix(Grid_Point_Counter,1,mxREAL);
  mxArray *mex_alt  = mxCreateDoubleMatrix(Grid_Point_Counter,1,mxREAL);
  mxArray *mex_mask = mxCreateDoubleMatrix(Grid_Point_Counter,1,mxREAL);
  mxArray *mex_BT_Counter = mxCreateDoubleMatrix(Grid_Point_Counter,1,mxREAL);
  mxArray *mex_flag = mxCreateCellMatrix(Grid_Point_Counter,1);
  mxArray *mex_BT_re = mxCreateCellMatrix(Grid_Point_Counter,1);
  mxArray *mex_BT_im = mxCreateCellMatrix(Grid_Point_Counter,1);
  mxArray *mex_pixrad = mxCreateCellMatrix(Grid_Point_Counter,1);
  mxArray *mex_incang = mxCreateCellMatrix(Grid_Point_Counter,1);
  mxArray *mex_aziang = mxCreateCellMatrix(Grid_Point_Counter,1);
  mxArray *mex_fayang = mxCreateCellMatrix(Grid_Point_Counter,1);
  mxArray *mex_geoang = mxCreateCellMatrix(Grid_Point_Counter,1);
  mxArray *mex_idpix = mxCreateCellMatrix(Grid_Point_Counter,1);
  mxArray *mex_foot1 = mxCreateCellMatrix(Grid_Point_Counter,1);
  mxArray *mex_foot2 = mxCreateCellMatrix(Grid_Point_Counter,1);
  mxArray *mex_BTxy  = mxCreateCellMatrix(Grid_Point_Counter,1);
  mxArray *mex_BThv  = mxCreateCellMatrix(Grid_Point_Counter,1);
  mxArray *mex_RA    = mxCreateCellMatrix(Grid_Point_Counter,1);
  mxArray *mex_snapID= mxCreateCellMatrix(Grid_Point_Counter,1);

  double *mat_TBdscr= new double[Grid_Point_Counter*STOKES_VEC*INCA_NUM];
  mxArray *mex_INCAdscr = mxCreateDoubleMatrix(INCA_NUM, 1, mxREAL);

  inn = 0;
  for(i=0;i<Grid_Point_Counter;++i){

    INMYBOX = false;   // it is set up to FALSE for every iteration

    mxArray *vec_flag=NULL, *vec_BTre=NULL, *vec_BTim=NULL, *vec_pixrad=NULL,  *vec_incang=NULL;
    mxArray *vec_aziang=NULL, *vec_fayang=NULL, *vec_geoang=NULL, *vec_idpix=NULL, *vec_foot1=NULL, *vec_foot2=NULL;
    mxArray *vec_BTxy=NULL, *vec_BThv=NULL, *vec_RA=NULL, *vec_snapID=NULL;

    fp.read((char *) &Grid_Point_ID, sizeof(Grid_Point_ID));
    fp.read((char *) &Grid_Point_Latitude, sizeof(float));
    fp.read((char *) &Grid_Point_Longitude, sizeof(float));
    fp.read((char *) &Grid_Point_Altitude , sizeof(float));
    fp.read((char *) &Grid_Point_Mask, sizeof(unsigned char));
    fp.read((char *) &BT_Data_Counter, sizeof(unsigned short int));

    if (Grid_Point_Latitude>=BOXLIM[0] && Grid_Point_Latitude<=BOXLIM[2] &&
    	Grid_Point_Longitude>=BOXLIM[1] && Grid_Point_Longitude<=BOXLIM[3]){
      INMYBOX = true;
      vec_flag   = mxCreateDoubleMatrix(BT_Data_Counter,1,mxREAL);
      vec_RA     = mxCreateDoubleMatrix(BT_Data_Counter,STOKES_VEC,mxREAL);
      vec_BTxy   = mxCreateDoubleMatrix(BT_Data_Counter,STOKES_VEC,mxREAL);
      vec_BThv   = mxCreateDoubleMatrix(BT_Data_Counter,STOKES_VEC,mxREAL);
      vec_BTre   = mxCreateDoubleMatrix(BT_Data_Counter,1,mxREAL);
      vec_BTim   = mxCreateDoubleMatrix(BT_Data_Counter,1,mxREAL);
      vec_pixrad = mxCreateDoubleMatrix(BT_Data_Counter,1,mxREAL);
      vec_incang = mxCreateDoubleMatrix(BT_Data_Counter,1,mxREAL);
      vec_aziang = mxCreateDoubleMatrix(BT_Data_Counter,1,mxREAL);
      vec_fayang = mxCreateDoubleMatrix(BT_Data_Counter,1,mxREAL);
      vec_geoang = mxCreateDoubleMatrix(BT_Data_Counter,1,mxREAL);
      vec_idpix  = mxCreateDoubleMatrix(BT_Data_Counter,1,mxREAL);
      vec_foot1 = mxCreateDoubleMatrix(BT_Data_Counter,1,mxREAL);
      vec_foot2 = mxCreateDoubleMatrix(BT_Data_Counter,1,mxREAL);
      vec_snapID= mxCreateDoubleMatrix(BT_Data_Counter,1,mxREAL);
    }  // enf of INMYBOX 0 block

    // initializing variables for interpolation:
    double Txx_in[BT_Data_Counter];   // = new double[BT_Data_Counter];
    double Tyy_in[BT_Data_Counter];   //    = new double[BT_Data_Counter];
    double Txy_in[BT_Data_Counter];   //    = new double[BT_Data_Counter];
    double BTxx_in[BT_Data_Counter];   //   = new double[BT_Data_Counter];
    double BTyy_in[BT_Data_Counter];   //   = new double[BT_Data_Counter];
    double BTxyRE_in[BT_Data_Counter];   // = new double[BT_Data_Counter];
    double BTxyIM_in[BT_Data_Counter];   // = new double[BT_Data_Counter];
    double T_out[BT_Data_Counter];   //     = new double[BT_Data_Counter];
    double BTxx_out[BT_Data_Counter];   //  = new double[BT_Data_Counter];
    double BTyy_out[BT_Data_Counter];   //  = new double[BT_Data_Counter];
    double BTxyRE_out[BT_Data_Counter];   // = new double[BT_Data_Counter];
    double BTxyIM_out[BT_Data_Counter];   // = new double[BT_Data_Counter];
    double RAxx_in[BT_Data_Counter];   //   = new double[BT_Data_Counter];
    double RAyy_in[BT_Data_Counter];   //   = new double[BT_Data_Counter];
    double RAxy_in[BT_Data_Counter];   //   = new double[BT_Data_Counter];
    double RAxx_out[BT_Data_Counter];   //  = new double[BT_Data_Counter];
    double RAyy_out[BT_Data_Counter];   //  = new double[BT_Data_Counter];
    double RAxy_out[BT_Data_Counter];   //  = new double[BT_Data_Counter];
    int Nxx=0, Nyy=0, Nxy=0;

    double theta_inc[BT_Data_Counter];
    double A_vec[4], X_vec[STOKES_VEC];  // system: M*X = A --> X unknown.
    double alpha[BT_Data_Counter];;  // rotational angle
    double filter1_ax1ax2[BT_Data_Counter], filter2_ax1ax2[BT_Data_Counter], filter3_tbxtby;

    for(j=0;j<BT_Data_Counter;++j){
      fp.read((char *) &Flags, sizeof(unsigned short int));
      fp.read((char *) &BT_Value_Real, sizeof(float));
      fp.read((char *) &BT_Value_Imag, sizeof(float));
      fp.read((char *) &Pixel_Radiometric_Accuracy, sizeof(unsigned short int));
      fp.read((char *) &Incidence_Angle, sizeof(unsigned short int));
      fp.read((char *) &Azimuth_Angle, sizeof(unsigned short int));
      fp.read((char *) &Faraday_Rotation_Angle, sizeof(unsigned short int));
      fp.read((char *) &Geometric_Rotation_Angle, sizeof(unsigned short int));
      fp.read((char *) &Snapshot_ID_of_Pixel, sizeof(unsigned int));
      fp.read((char *) &Footprint_Axis1, sizeof(unsigned short int));
      fp.read((char *) &Footprint_Axis2, sizeof(unsigned short int));
      // Quality Control Checking:
      if (Flags&0x80) continue;   // SUN_POINT flag 8th bit=1, pixel degraded
      // criteria to filter based on footprint spatial resolution requirements,
      filter1_ax1ax2[j] = ((double) Footprint_Axis1)/((double) Footprint_Axis2);
      filter2_ax1ax2[j] = 2*PIX_FOOTP_SC*sqrt((double) (Footprint_Axis1*Footprint_Axis2));
      // RFI filter based on 4th Stoks vector ST4>TH_RFI_ST4 indicates RFI
      bool RFI_ST4 = TH_RFI_ST4<abs(-2*BT_Value_Imag)? true : false;

      if (INMYBOX){

	int idx_snap = distance(Snapshot_ID,
				find(Snapshot_ID,Snapshot_ID+Snapshot_Counter,
				     Snapshot_ID_of_Pixel));
	// getting the time for the snapshot [sec]:
	if(inn==0 && j==0) smos::TimeEEC2date(Snapshot_Time[idx_snap],smos_date0);

	T_out[j] = (double) (Snapshot_Time[idx_snap][1]+
			     Snapshot_Time[idx_snap][2]*1E-6);

	*(mxGetPr(vec_flag)+j) = Flags;
	*(mxGetPr(vec_snapID)+j) = idx_snap+1;

	// Flags' 15 and 16 bit = FALSE -> XX-pol:
	if (((~Flags)&0x01) && ((~Flags)&0x02) && !(Flags&0x40) && !RFI_ST4){
	  Txx_in[Nxx]  = T_out[j];
	  BTxx_in[Nxx] = (double) BT_Value_Real;
	  RAxx_in[Nxx] = (double) (Pixel_Radiometric_Accuracy*RAD_ACC_SC);
	  Nxx++;
	}
	// Flags' 15 bit = FALSE and 16 bit = TRUE -> YY-pol:
	if (( Flags&0x01) && ((~Flags)&0x02) && !(Flags&0x4000)  && !RFI_ST4){
	  Tyy_in[Nyy]  = T_out[j];
	  BTyy_in[Nyy] = (double) BT_Value_Real;
	  RAyy_in[Nyy] = (double) (Pixel_Radiometric_Accuracy*RAD_ACC_SC);
	  Nyy++;
	}
	// Flags' 15 bit = TRUE -> XY-pol:
	if (Flags&0x02  && !RFI_ST4){
	  Txy_in[Nxy]  = T_out[j];
	  BTxyRE_in[Nxy] = (double) BT_Value_Real;
	  BTxyIM_in[Nxy] = (double) BT_Value_Imag;
	  RAxy_in[Nxy]   = (double) (Pixel_Radiometric_Accuracy*RAD_ACC_SC);
	  Nxy++;
	}
	theta_inc[j] = Incidence_Angle*FAK90DEG;

	*(mxGetPr(vec_BTre)+j) = BT_Value_Real;
	*(mxGetPr(vec_BTim)+j) = BT_Value_Imag;
	*(mxGetPr(vec_pixrad)+j) = Pixel_Radiometric_Accuracy*RAD_ACC_SC;
	*(mxGetPr(vec_incang)+j) = theta_inc[j];
	*(mxGetPr(vec_aziang)+j) = Azimuth_Angle*FAK360DEG;
	*(mxGetPr(vec_fayang)+j) = Faraday_Rotation_Angle*FAK360DEG;
	*(mxGetPr(vec_geoang)+j) = Geometric_Rotation_Angle*FAK360DEG;
	*(mxGetPr(vec_idpix)+j)  = Snapshot_ID_of_Pixel;
	*(mxGetPr(vec_foot1)+j)  = Footprint_Axis1*PIX_FOOTP_SC;
	*(mxGetPr(vec_foot2)+j)  = Footprint_Axis2*PIX_FOOTP_SC;
	alpha[j] = (Faraday_Rotation_Angle + Geometric_Rotation_Angle)*FAKRAD;
	if(CREATOR_VER<344) alpha[j] += -2*Faraday_Rotation_Angle*FAKRAD; //old version

      } // end of INMYBOX block
    }   // end of BT_Data_Counter block (j index)

    if (INMYBOX){

      // Here local variables for TB with incident angle discrete:
      // with c++11 TB_dscr may be initialized as {{init,step,Ndcsr},{},...};
      smos::Be_Discrete TB_dscr[STOKES_VEC]; //={{0,1,46},{0,1,46},{0,1,46},{0,1,46}};
      int Ndscr = INCA_NUM;   // default SMOS spec. from 0 to 65 deg (66 angles).

      memcpy(mxGetPr(mex_INCAdscr), TB_dscr[0].theta,Ndscr*sizeof(double));

      // **** Interpolation:
      switch (STOKES_VEC){
      case 4:
	smos::interp1(Txy_in,BTxyRE_in,Nxy,T_out,BTxyRE_out,BT_Data_Counter);
	smos::interp1(Txy_in,BTxyIM_in,Nxy,T_out,BTxyIM_out,BT_Data_Counter);
	smos::interp1(Txy_in,RAxy_in,Nxy,T_out,RAxy_out,BT_Data_Counter);
      case 2:
	smos::interp1(Txx_in,BTxx_in,Nxx,T_out,BTxx_out,BT_Data_Counter);
	smos::interp1(Tyy_in,BTyy_in,Nyy,T_out,BTyy_out,BT_Data_Counter);
	smos::interp1(Txx_in,RAxx_in,Nxx,T_out,RAxx_out,BT_Data_Counter);
	smos::interp1(Tyy_in,RAyy_in,Nyy,T_out,RAyy_out,BT_Data_Counter);
	break;
      default:
	mexErrMsgTxt("Problem with declaration of STOKES VECTOR: need to be 4 or 2!");
      }
      // ****

      // calculating the average of the halved 1st Stokes parameter <ST1>:
      float aveST1;
      int k;
      for(j=0, k=0, aveST1=0;j<BT_Data_Counter;++j){
	if(std::isnan(BTxx_out[j]) || std::isnan(BTyy_out[j])) continue;
	aveST1 += 0.5*(BTxx_out[j]+BTyy_out[j]);
	k++;
      }
      aveST1 /= k;

      // **** Rotation matrix:
      for(j=0;j<BT_Data_Counter;++j){

	// *** Here all filetering criteria ST4_RFI, ST1 and FootPrint:
	bool AX1AX2 = false, TBXTBY = false, RFI_ST1 = false;

	if(filter1_ax1ax2[j]>TH_ELON || filter2_ax1ax2[j]>TH_SIZE) AX1AX2 = true;
	filter3_tbxtby = sqrt(BTxx_out[j]*BTxx_out[j]+BTyy_out[j]*BTyy_out[j]);
	if(filter3_tbxtby<TBXX_MIN || filter3_tbxtby>TBXX_MAX) TBXTBY=true;
	if((0.5*(BTxx_out[j]+BTyy_out[j])-aveST1)>(CA_TBS1+CB_TBS1*RAxx_out[j])) RFI_ST1=true;
	// ***

	if(!TBXTBY && !AX1AX2 && !RFI_ST1){
	  A_vec[0] = BTxx_out[j]; A_vec[1] = BTyy_out[j];
	  A_vec[2] = 2*BTxyRE_out[j]; A_vec[3] = -2*BTxyIM_out[j];

	  smos::MxTrans_XY2HV(alpha[j], A_vec, X_vec, STOKES_VEC, false); // rotation
	}
	else {
	  for(unsigned k1=0; k1<STOKES_VEC; ++k1) X_vec[k1] = MYNAN;
	}

	for(unsigned k1=0; k1<STOKES_VEC; ++k1){
	  // Feeding the "Be_Discrete" class type variables:
	  TB_dscr[k1].consider_it(theta_inc[j], X_vec[k1]);
	  // filling the HV output to mex matrix:
	  *(mxGetPr(vec_BThv)+j+k1*BT_Data_Counter) = X_vec[k1];
	}
	// TB_dscr[1].consider_it(theta_inc[j], X_vec[1]);
	// TB_dscr[2].consider_it(theta_inc[j], X_vec[2]);
	// TB_dscr[3].consider_it(theta_inc[j], X_vec[3]);
	// // filling the HV output to mex matrix:
	// *(mxGetPr(vec_BThv)+j)                   = X_vec[0];
	// *(mxGetPr(vec_BThv)+j+BT_Data_Counter)   = X_vec[1];
	// *(mxGetPr(vec_BThv)+j+2*BT_Data_Counter) = X_vec[2];
	// *(mxGetPr(vec_BThv)+j+3*BT_Data_Counter) = X_vec[3];
	// -----
	// for Radiometric_Accuracy RA:
	if(!TBXTBY && !AX1AX2 && !RFI_ST1){
	  A_vec[0] = fabs(RAxx_out[j]); A_vec[1] = fabs(RAyy_out[j]);
	  A_vec[2] = 2*fabs(RAxy_out[j]); A_vec[3] = 2*fabs(RAxy_out[j]);

	  smos::MxTrans_XY2HV(alpha[j], A_vec, X_vec, STOKES_VEC, true);
	}
	else {
	  for(unsigned k1=0; k1<STOKES_VEC; ++k1) X_vec[k1] = MYNAN;
	}
	// filling the Radiometric accuracy HV to mex matrix:
	for(unsigned k1=0; k1<STOKES_VEC; ++k1)
	  *(mxGetPr(vec_RA)+j+k1*BT_Data_Counter) = X_vec[k1];
	// *(mxGetPr(vec_RA)+j+BT_Data_Counter)   = X_vec[1];
	// *(mxGetPr(vec_RA)+j+2*BT_Data_Counter) = X_vec[2];
	// *(mxGetPr(vec_RA)+j+3*BT_Data_Counter) = X_vec[3];
      }    // end over BT_Data_Counter (number of incident angles)

      // Feeding matrix with the TBs discrete incidence angle:
      // inn: snapshot index, Ndscr: # of angles, STOKES_VEC: dimensions
      unsigned k0  = inn*Ndscr*STOKES_VEC;
      for(unsigned k1 = 0; k1<STOKES_VEC; k0 = ++k1*Ndscr+inn*Ndscr*STOKES_VEC)
	memcpy((void *) (mat_TBdscr+k0), TB_dscr[k1].deliver_it(),Ndscr*sizeof(double));

      // Populating the mex vector variables:
      memcpy(mxGetPr(vec_BTxy), BTxx_out, BT_Data_Counter*sizeof(double));
      memcpy(mxGetPr(vec_BTxy)+BT_Data_Counter, BTyy_out, BT_Data_Counter*sizeof(double));
      if(STOKES_VEC==4){
	memcpy(mxGetPr(vec_BTxy)+2*BT_Data_Counter, BTxyRE_out,
	       BT_Data_Counter*sizeof(double));
	memcpy(mxGetPr(vec_BTxy)+3*BT_Data_Counter, BTxyIM_out,
	       BT_Data_Counter*sizeof(double));
      }
      *(mxGetPr(mex_gridid)+inn) = Grid_Point_ID;
      *(mxGetPr(mex_lat)+inn) = Grid_Point_Latitude;
      *(mxGetPr(mex_lon)+inn) = Grid_Point_Longitude;
      *(mxGetPr(mex_alt)+inn) = Grid_Point_Altitude;
      *(mxGetPr(mex_mask)+inn) = Grid_Point_Mask;
      *(mxGetPr(mex_BT_Counter)+inn) = BT_Data_Counter;

      // Populating the mex cell variables
      mxSetCell(mex_flag, inn, vec_flag);
      mxSetCell(mex_BT_re, inn, vec_BTre);
      mxSetCell(mex_BT_im, inn, vec_BTim);
      mxSetCell(mex_pixrad, inn, vec_pixrad);
      mxSetCell(mex_incang, inn, vec_incang);
      mxSetCell(mex_aziang, inn, vec_aziang);
      mxSetCell(mex_fayang, inn, vec_fayang);
      mxSetCell(mex_geoang, inn, vec_geoang);
      mxSetCell(mex_idpix, inn, vec_idpix);
      mxSetCell(mex_foot1, inn, vec_foot1);
      mxSetCell(mex_foot2, inn, vec_foot2);
      mxSetCell(mex_BTxy, inn, vec_BTxy);
      mxSetCell(mex_BThv, inn, vec_BThv);
      mxSetCell(mex_RA  , inn, vec_RA);
      mxSetCell(mex_snapID, inn, vec_snapID);
      inn++;

    }   // end of INMYBOX block_2

  }
  fp.close();

  // --- Assigning output variables:

  // Setting MATLAB workspace variable mexFile_Name with info for DBL file name:
  mxArray *mexFile_Name = mxCreateString(WhatPol.substr(WhatPol.find_last_of("/")+1).c_str());

  mxArray *mex_boxlim = mxCreateDoubleMatrix(2,2,mxREAL);
  for(i=0;i<4;++i) *(mxGetPr(mex_boxlim)+i) = BOXLIM[i];

  int idx_snap = distance(Snapshot_ID,
			  find(Snapshot_ID,Snapshot_ID+Snapshot_Counter,
			       Snapshot_ID_of_Pixel));
  smos::TimeEEC2date(Snapshot_Time[idx_snap],smos_date1);  // time for last point in BOXLIM

  char *str[] = {smos_date0, smos_date1};   // {Start_Time, End_Time}
  mxArray *mex_startdate = mxCreateCharMatrixFromStrings(2, (const char **)str);

  // --- Setting the actual size of variables according to index: inn
  mxSetM(mex_gridid,inn);
  mxSetM(mex_lat,inn);
  mxSetM(mex_lon,inn);
  mxSetM(mex_alt,inn);
  mxSetM(mex_mask,inn);
  mxSetM(mex_BT_Counter,inn);
  mxSetM(mex_flag,inn);
  mxSetM(mex_BT_re,inn);
  mxSetM(mex_BT_im,inn);
  mxSetM(mex_pixrad,inn);
  mxSetM(mex_incang,inn);
  mxSetM(mex_aziang,inn);
  mxSetM(mex_fayang,inn);
  mxSetM(mex_geoang,inn);
  mxSetM(mex_idpix,inn);
  mxSetM(mex_foot1,inn);
  mxSetM(mex_foot2,inn);
  mxSetM(mex_BTxy,inn);
  mxSetM(mex_BThv,inn);
  mxSetM(mex_RA  ,inn);
  mxSetM(mex_snapID,inn);

  // *NOTE: The following lines get the output as Array [INC_ANGLE,POLARIZATION,NUM_DATA]
  // mwSize dims3d[3]={INCA_NUM,STOKES_VEC,(int) inn};
  // mxArray *mex_TBdscr = mxCreateNumericArray(3,dims3d,mxDOUBLE_CLASS,mxREAL);
  // memcpy(mxGetPr(mex_TBdscr), (const void *) mat_TBdscr,INCA_NUM*inn*STOKES_VEC*sizeof(double));

  // *NOTE: The following lines get the output as Array [NUM_DATA,INC_ANGLE,POLARIZATION]
  // this dimension arrange is more convinient for MATLAB Workspace data treatment.
  mwSize dims3d[3] = {(mwSize) inn, (mwSize) INCA_NUM, (mwSize) STOKES_VEC};
  mxArray *mex_TBdscr = mxCreateNumericArray(3,dims3d,mxDOUBLE_CLASS,mxREAL);

  // Arranging mat_TBdscr array to dimensions {NUM_DATA, INC_ANGLE, POLARIZATION}
  for(unsigned k1=0,idat=0;k1<STOKES_VEC;++k1)
    for(unsigned k2=0;k2<INCA_NUM; ++k2)
      for(unsigned k3=0;k3<inn;k3++,	idat++)
  	*(mxGetPr(mex_TBdscr)+idat) = mat_TBdscr[k1*INCA_NUM+k2+k3*(STOKES_VEC*INCA_NUM)];

  // deleting the temporal variable mat_TBdscr:
  delete[] mat_TBdscr;
  // ---
  // Create the output structure:
  mwSize dims[2]={1,1};
  TempSwathFull = mxCreateStructArray(2, dims, N_field_name, field_name);

  // Setting the variables into the Output Structure:
  mxSetField(TempSwathFull,0,"File_Name", mexFile_Name );
  mxSetField(TempSwathFull,0,"GEOBoxLimits", mex_boxlim);
  mxSetField(TempSwathFull,0,"Start_End_Time", mex_startdate);
  mxSetField(TempSwathFull,0,"GridPoint_ID", mex_gridid);
  mxSetField(TempSwathFull,0,"GridPoint_Latitude",mex_lat);
  mxSetField(TempSwathFull,0,"GridPoint_Longitude",mex_lon);
  mxSetField(TempSwathFull,0,"GridPoint_Altitude",mex_alt);
  mxSetField(TempSwathFull,0,"GridPoint_Mask", mex_mask);
  mxSetField(TempSwathFull,0,"BT_Data_Counter",mex_BT_Counter);
  mxSetField(TempSwathFull,0,"flags", mex_flag);
  mxSetField(TempSwathFull,0,"BTvalue_re",mex_BT_re);
  mxSetField(TempSwathFull,0,"BTvalue_im",mex_BT_im);
  mxSetField(TempSwathFull,0,"PixelRadiometry_acc", mex_pixrad);
  mxSetField(TempSwathFull,0,"Incidence_angle",mex_incang);
  mxSetField(TempSwathFull,0,"Azimuth_angle",mex_aziang);
  mxSetField(TempSwathFull,0,"FaradayRotation_angle",mex_fayang);
  mxSetField(TempSwathFull,0,"GeometricRotation_angle",mex_geoang);
  mxSetField(TempSwathFull,0,"Snapshot_ID_Pixels", mex_idpix);
  mxSetField(TempSwathFull,0,"Footprint_Axis1",mex_foot1);
  mxSetField(TempSwathFull,0,"Footprint_Axis2",mex_foot2);
  mxSetField(TempSwathFull,0,"TBxy",mex_BTxy);
  mxSetField(TempSwathFull,0,"TBhv",mex_BThv);
  mxSetField(TempSwathFull,0,"RA",mex_RA);
  mxSetField(TempSwathFull,0,"idx_SnapshotID",mex_snapID);
  mxSetField(TempSwathFull,0,"TB_Fixed_IncAngle",mex_TBdscr);
  mxSetField(TempSwathFull,0,"Fixed_IncAngle", mex_INCAdscr);

  SnapShotInfo = mxCreateStructArray(2, dims, N_snapsh_name, snapsh_name);
  if (smos::SSI_FLAG){
    mxSetField(SnapShotInfo,0,"Snapshot_Time",mex_snaptime);
    mxSetField(SnapShotInfo,0,"Snapshot_ID",mex_snapid);
    mxSetField(SnapShotInfo,0,"Snapshot_OBET",mex_snapobet);
    mxSetField(SnapShotInfo,0,"XYZ_Position",mex_xyzpo);
    mxSetField(SnapShotInfo,0,"XYZ_Velocity",mex_xyzve);
    mxSetField(SnapShotInfo,0,"Vector_So",mex_vecsou);
    mxSetField(SnapShotInfo,0,"Qnn",mex_Qnn);
    mxSetField(SnapShotInfo,0,"TEC",mex_tec);
    mxSetField(SnapShotInfo,0,"Geomag_F_D_I",mex_geomag);
    mxSetField(SnapShotInfo,0,"Sun_RA_DEC_BT",mex_sunv);
    mxSetField(SnapShotInfo,0,"Accuracy",mex_accu);
    mxSetField(SnapShotInfo,0,"Radiometric_Accuracy",mex_radaccu);
    mxSetField(SnapShotInfo,0,"XBand",mex_xband);
    mxSetField(SnapShotInfo,0,"SIAC_Error_flag",mex_errflag);
  }

  // ***********************************************************
  // IF output directory included, creating MAT file:
  if(nrhs==3 && inn!=0){
    if(smos::CreateMAToutput(TempSwathFull,SnapShotInfo, OUTDIR, filen)==0)
      mexErrMsgTxt("Output MAT file could not be created!.");
  }  // end block creating MAT file

  // ***********************************************
  // Returning Variables to MATLAB workspace
  switch(nlhs){
  case 2:
    plhs[1] = SnapShotInfo;
  case 1:
    plhs[0] = TempSwathFull;
  }

  return;
}
// ******** END OF MAIN PROGRAM ***************
// --------------------------------------------
// ******** AUXILIARY ROUTINES ***************
// *******************************************
//

namespace smos{

  // ----------------------------------------------------------
  // Routine to:
  // 1) open a GUI file browser to select a DBL input file
  // to work with.
  // 2) open a GUI input dialog to write the minimum and maximum
  // Latitude-Longitude pairs to select from the DBL input file.
  // ARGUMENTS:
  // * filen:  input DBL file name (output arg).
  // * OUTDIR: directory where filen is located (output arg).
  // * BOXLIM: 4-element vector (lat_min,lon_min,lat_max,lon_max) (output atg).
  // RETURN: status=0 -> OK, status!=0 -> wrong procedure.
  int GetInputDBL_BoxLIM(char *& filen, char *& OUTDIR, float BOXLIM[4]){
    int strLength, status;
    mxArray *INVAR, *OUTVAR[2];

    smos::ShowGNUPL();      // displaying License, it is free!.
    INVAR = mxCreateString("*.DBL");
    status = mexCallMATLAB(2,OUTVAR,1,&INVAR,"uigetfile");
    if (status!=0) mexErrMsgTxt("File selection not possible!");
    // passing the input and output path
    strLength = mxGetN(OUTVAR[0])+mxGetN(OUTVAR[1])+1;
    filen = (char *) mxCalloc(strLength, sizeof(char));
    mxGetString(OUTVAR[1],filen,strLength);
    // passing the file name
    strLength = mxGetN(OUTVAR[0])+1;
    mxGetString(OUTVAR[0],filen+mxGetN(OUTVAR[1]),strLength);
    strLength = mxGetN(OUTVAR[1])+1;
    OUTDIR = (char *) mxCalloc(strLength,sizeof(char));
    mxGetString(OUTVAR[1], OUTDIR, strLength);
    mexPrintf("DBL file chosen: %s\n",filen);
#ifdef MATLAB_MEX_FILE
    // Getting the Latitude and Longitude Box limits:
    mxArray *OUTbox=mxCreateCellMatrix(1,4);
    mxArray *prompt = mxCreateCellMatrix(1,4);
    const mxArray *name    = mxCreateString("Lat-Lon BOX Limits");
    const mxArray *numlines = mxCreateDoubleMatrix(1,1,mxREAL);
    mxArray *defaultanswer = mxCreateCellMatrix(1,4);
    const mxArray *option = mxCreateString("on");
    const mxArray *INOPTS[5] = {prompt,name,numlines,defaultanswer,option};
    char temp_str[4];

    *mxGetPr(numlines) = 1;
    mxSetCell(prompt,0,mxCreateString("Latitude bottom left [-90:+90]"));
    mxSetCell(prompt,1,mxCreateString("Longitude bottom left [-180:+180]"));
    mxSetCell(prompt,2,mxCreateString("Latitude upper rigth [-90:+90]"));
    mxSetCell(prompt,3,mxCreateString("Longitude upper rigth [-180:+180]"));

    mxSetCell(defaultanswer,0,mxCreateString("-90"));
    mxSetCell(defaultanswer,1,mxCreateString("-180"));
    mxSetCell(defaultanswer,2,mxCreateString("90"));
    mxSetCell(defaultanswer,3,mxCreateString("180"));

    status = mexCallMATLAB(1,&OUTbox,5,(mxArray **) INOPTS,"inputdlg");
    if (status!=0) mexErrMsgTxt("(Lat,Lon) Limits selection wrong!");

    for(int i=0;i<4;++i){
      strLength = mxGetN(mxGetCell(OUTbox,i))+1;
      mxGetString(mxGetCell(OUTbox,i),temp_str,strLength);
      BOXLIM[i]=atof(temp_str);
    }
#else
    cout<<"Please introduce Lat-Lon BOX Limits"<<endl;
    cout<<"Latitude bottom left [-90:+90]   : "; cin>>BOXLIM[0];
    cout<<"Longitude bottom left [-180:+180]: "; cin>>BOXLIM[1];
    cout<<"Latitude upper rigth [-90:+90]   : "; cin>>BOXLIM[2];
    cout<<"Longitude upper rigth [-180:+180]: "; cin>>BOXLIM[3];
#endif
    if(BOXLIM[0]>BOXLIM[2] || BOXLIM[1]>BOXLIM[3] ||
       (BOXLIM[0]<-90 || BOXLIM[0]>90) ||
       (BOXLIM[2]<-90 || BOXLIM[2]>90) ||
       (BOXLIM[1]<-180 || BOXLIM[1]>180) ||
       (BOXLIM[3]<-180 || BOXLIM[3]>180)){
      mexPrintf("Wrong set of min-max parameters (%5.1f,%5.1f) - (%5.1f,%5.1f)\n",
		BOXLIM[0],BOXLIM[1],BOXLIM[2],BOXLIM[3]);
      return(1);}

    return(0);
  }

  // --------------------------------------------------------------
  // Linear Interpolation subroutine:
  // Values at (Xo,Yo) are interpolated based on the (Xi,Yi) pairs.
  // Ni: number of input  data for (Xi[j],Yi[j]), with j=0,...,Ni-1
  // No: number of output data for (Xo[j],Yo[j]), with j=0,...,No-1
  // NOTE 1: In case Xo[j] is out of the Xi lower/upper limits,
  // then Yo[j] is assigned with a non-interpolated value of MYNAN;
  // NOTE 2: In case Ni is less than 6 points, Yo is returned MYNAN
  // (no extrapolation allowed).

  void interp1(double Xi[], double Yi[], int Ni, double Xo[], double Yo[], int No){

    int i,k;
    // Interpolation attempt only when input has more than 6 points!
    if(Ni<7) for(i=0,k=0;i<No;++i) Yo[i] = Xi[k]==Xo[i] ? Yi[k++] : MYNAN;
    else{
      gsl_interp_accel *acc = gsl_interp_accel_alloc ();
      gsl_interp *inter = gsl_interp_alloc (gsl_interp_linear, Ni); //cspline, Ni);
      gsl_interp_init(inter, Xi, Yi, Ni);
      for(i=0;i<No;++i){
	if(Xo[i]>=Xi[0] && Xo[i]<=Xi[Ni-1])
	  Yo[i] = gsl_interp_eval(inter,Xi,Yi, Xo[i], acc);
	else Yo[i]=MYNAN;
      }
      gsl_interp_free(inter);
      gsl_interp_accel_free(acc);
    }
    return;
  }


  // Matrix solver for the system M*X = B, with X the unknown,
  // and M a square matrix of order N. X and B Nx1 arrays.
  // The solver uses LU decompositions and it is applied directly
  // to solve in case of the TBxy. For case of Radiometri Accuracies (RA)
  // matrix M needsto be inverted and then calculate the absolute of M.
  // NOTE: N is either 4 (full-polarization) or 2 (dual-polarization).
  void MxTrans_XY2HV(double alpha, double A[], double Xout[], int N, bool RA){

    // defining rotational Matrix,
    // where alpha is the Geometrical angle + Faraday rotation angle [RADIANS]
    // N: number of elements of vector A. Needs to be either 4 (full-pol) or 2 (dual-pol).
    double MR4[] = {
      pow(cos(alpha),2), pow(sin(alpha),2),-cos(alpha)*sin(alpha), 0, // 1st row
      pow(sin(alpha),2), pow(cos(alpha),2), cos(alpha)*sin(alpha), 0, // 2nd row
      sin(2*alpha), -sin(2*alpha), cos(2*alpha), 0, // 3th row
      0, 0, 0, 1 // 4th row
    };

    double MR2[] = {
      pow(cos(alpha),2), pow(sin(alpha),2), // 1st row
      pow(sin(alpha),2), pow(cos(alpha),2)  // 2nd row
    };

    gsl_matrix_view m = gsl_matrix_view_array (N==4?MR4:MR2, N, N);

    gsl_vector_view b = gsl_vector_view_array (A, N);

    gsl_vector *x = gsl_vector_alloc(N);

    int s;

    gsl_permutation * p = gsl_permutation_alloc (N);

    gsl_linalg_LU_decomp (&m.matrix, p, &s);

    if(RA){  // RA=true for Radiometric Accuracies only
      double MR4_inv[N*N];
      gsl_matrix_view m_inv = gsl_matrix_view_array(MR4_inv,N,N);
      gsl_linalg_LU_invert (&m.matrix, p, &m_inv.matrix);
      for(int k=0;k<N*N;++k) MR4_inv[k] = fabs(MR4_inv[k]);
      // solving system M*X=b as X=abs(inv(M))*B
      gsl_blas_dgemv (CblasNoTrans, 1.0, &m_inv.matrix, &b.vector,0.0, x);
    }
    else
      gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);

    for(int i=0; i<N; ++i) Xout[i] = gsl_vector_get(x,i);

    gsl_permutation_free (p);
    gsl_vector_free (x);
    return;
  }

  // Obtain the specific field's value from the corresponding HDR file
  // given as input the DBL file. Both DBL and HDR files have identical
  // name structure.
  // NOTE: this is not a XML reader, but it only works for fields
  // given as: <Field>value</Field>
  // in the same line within the HDR input file.
  // USAGE: value = getHDR_fields(DBL_filename,Fiels);
  float getHDR_fields(const char filename[], const char subtxt[]){
    ifstream file;
    string str;
    string txt = subtxt;
    string txt_0 = "<"+txt+">";
    string txt_1 = "</"+txt+">";
    size_t idx_txt_0, idx_txt_1;
    float value;

    // Converting input file_name form DBL to HDR:
    string hdr_file=filename;
    hdr_file.replace(hdr_file.length()-3,3,"HDR");
    // Opening the Header file HDL:
    file.open(hdr_file.c_str(),ios::in);
    if(!file.is_open()){
      cerr<<"HDL file not found! "<<hdr_file<<endl;
      return(MYNAN);
    }
    while (getline(file, str))
      {
	// Finding the Field?
	idx_txt_0 = str.find(txt_0);
	if(idx_txt_0 != string::npos) {
	  // cout<<"Found!"<<endl;
	  idx_txt_1 = str.find(txt_1);
	  idx_txt_0 += txt_0.size();
	  idx_txt_1 -= idx_txt_0;

	  // Converting field value from string to intefer:
	  value = atof(str.substr(idx_txt_0,idx_txt_1).c_str());
	  file.close();
	  return(value);
	}
      }
    if(file.eof()) cerr<<"EOF! Field not found: "<<subtxt<<endl;
    file.close();
    return(MYNAN);
  }

  // ************************************************************
  // Routice to create the MAT output file.
  // Inputs are: TempSwathFull (TSF), SnapshotInfo (SSI), output directory (OURDIR)
  // and generic DBL input file (filen)
  // Output is a integer indicating 0 if success or -1 otherwise
  // The output MAT file has the same name as DBL input file but with ".mat" instead of
  //  ".DBL", and is located at OUTDIR. In case no OUTDIR is specified, then no MAT
  // file is created.
  int CreateMAToutput(mxArray *TSF, mxArray *SSI, char OUTDIR[], char filen[]){

    cout<<"Creating MAT output file..."<<endl;
    // making up the output file with outdir
    int FileLength = sizeof(OUTDIR)-2;
    if(OUTDIR[FileLength] != '/') strcat(OUTDIR,"/");

    string mat_file=filen;
    int bslash=mat_file.find_last_of("/");
    mat_file.replace(0,bslash+1,OUTDIR);
    mat_file.replace(mat_file.length()-3,3,"mat");

#ifdef MATLAB_MEX_FILE
    MATFile *pmat;
    pmat = matOpen(mat_file.c_str(), "w");
    if (pmat == NULL) {
      mexErrMsgTxt("Error creating MAT file");
      cout<<mat_file.c_str()<<endl;
      return 0;
    }
    if(matPutVariable(pmat, "TSF", TSF)!=0){
      mexErrMsgTxt("Error storing Temp_Swath_Full");
      return 0;
    }
    if(smos::SSI_FLAG){
      if(matPutVariable(pmat, "SSI", SSI)!=0){
	mexErrMsgTxt("Error storing Snapshot_Information");
	return 0;
      }
    }
    if (matClose(pmat) != 0) {
      mexErrMsgTxt("Error closing MAT file.");
      return 0;
    }
#else  // OCTAVE_MEX_FILE
    matvar_t *TSFStruct;
    mat_t *matfp = Mat_CreateVer(mat_file.c_str(),NULL,MAT_FT_MAT73);
    if(NULL==matfp){
      mexPrintf("%s :", mat_file.c_str());
      mexWarnMsgTxt("Error creating OCTAVE MAT file");
      return 0;
    }
    TSFStruct = CloneSTRUCTforMATIO(TSF,"TSF");
    Mat_VarWrite(matfp,TSFStruct,MAT_COMPRESSION_ZLIB);
    Mat_VarFree(TSFStruct);
    if(smos::SSI_FLAG){
      mexPrintf("inside matio SSI_FLAG\n");
      matvar_t *SSIStruct = CloneSTRUCTforMATIO(SSI,"SSI");
      Mat_VarWrite(matfp,SSIStruct,MAT_COMPRESSION_ZLIB);
      Mat_VarFree(SSIStruct);
    }
    Mat_Close(matfp);

#endif
    return 1;
  }

  // *********************************************************
  // Routine to convert MATLAB structure variables into
  // matvar_t type variable suitable to be used under the
  // MATIO.h library to store MAT files with Octave.
  // Input is a mxArray structure
  // Ouput is a matvar_t structure clone.
#ifndef MATLAB_MEX_FILE
  matvar_t *CloneSTRUCTforMATIO(mxArray *TSF, const char *var_name){

    matvar_t *structTSF;
    size_t struct_dims[2] = {1,1};

    // Getting the elements of Structure input and Field Names:
    unsigned nF = mxGetNumberOfFields(TSF);
    const char *MyNameField[nF];
    for(unsigned k=0;k<nF;++k) MyNameField[k] = mxGetFieldNameByNumber(TSF,k);

    structTSF = Mat_VarCreateStruct(var_name,2,struct_dims,MyNameField,nF);
    if (NULL == structTSF) {
      mexWarnMsgTxt("Error creating struct variable for MATIO\n");
      Mat_VarFree(structTSF);
      return(NULL);
    }

    // the matrices are easy to convert

    for(unsigned k=0; k<nF; ++k){
      matvar_t *field;
      mxArray *temp = mxGetField(TSF,0,MyNameField[k]);
      mxClassID TempClass = mxGetClassID(temp);
      const mwSize *Tempdims = mxGetDimensions(temp);
      mwSize Nodims = mxGetNumberOfDimensions(temp);
      size_t dims[Nodims];

      for(int j=0; j<Nodims; ++j) dims[j] = *(Tempdims+j);

      switch (TempClass){
      case mxCHAR_CLASS:
	field = Mat_VarCreate(NULL,MAT_C_CHAR,MAT_T_UTF8,Nodims,dims,mxGetPr(temp),0);
	if (NULL == field){
	  mexPrintf("Char variable %s : ",MyNameField[k]);
	  mexWarnMsgTxt("cannot be assigned to structure!\n");
	  Mat_VarFree(structTSF);
	  return(NULL);
	}
	break;
      case mxDOUBLE_CLASS:
	field = Mat_VarCreate(NULL,MAT_C_DOUBLE,MAT_T_DOUBLE,Nodims,dims,mxGetPr(temp),0);
	if (NULL == field){
	  mexPrintf("Matrix variable %s : ",MyNameField[k]);
	  mexWarnMsgTxt("cannot be assigned to structure!\n");
	  Mat_VarFree(structTSF);
	  return(NULL);
	}
	break;
      case mxCELL_CLASS:
	field = Mat_VarCreate(NULL,MAT_C_CELL,MAT_T_CELL,Nodims,dims,NULL,0);
	if (NULL == field){
	  mexPrintf("cell variable %s : ",MyNameField[k]);
	  mexWarnMsgTxt("cannot be assigned to structure!\n");
	  Mat_VarFree(structTSF);
	  return(NULL);
	}
	else{
	  for(unsigned j=0; j<dims[0];++j){
	    mxArray *TempCell = mxGetCell(temp,j);
	    size_t Eledims[2];
	    Eledims[0] = mxGetM(TempCell);
	    Eledims[1] = mxGetN(TempCell);
	    matvar_t *cell_element = Mat_VarCreate(NULL,MAT_C_DOUBLE,MAT_T_DOUBLE,2,Eledims,mxGetPr(TempCell),0);
	    Mat_VarSetCell(field,j,cell_element);
	  }
	}
	break;
      default:
	{
	  mexWarnMsgTxt("Structure variable not identified!");
	  field = NULL;
	}
      }  // end of switch
      Mat_VarSetStructFieldByName(structTSF,MyNameField[k],0,field);

    } // end over loop number of variables in structure

    return(structTSF);
  }
#endif

  // *********************************************************
  // Converting the EE CFI time format to a string time format
  // with (ISO 8601)  'yyyymmddTHH:MM:SS'  e.g. 20040509T15:45:17
  // In EE CFI, days count begins on the 1st of January 2000 AD, while
  // seconds are relative to time for current day in UTC.
  int TimeEEC2date(int ee_cfi[], char smos_date[]){

    struct tm y2k;
    long int smos_sec;

    smos_sec = ee_cfi[0]*60*60*24+ee_cfi[1];  //microsec are not considered

    //Seting up the EE CFI initial date:
    y2k.tm_hour = 0; y2k.tm_min = 0; y2k.tm_sec = smos_sec;
    y2k.tm_year = 100;   // 2000 = 1900 + 100
    y2k.tm_mon = 0; y2k.tm_mday = 1; y2k.tm_isdst = -1;
    // mktime(&y2k) is needed, otherwise y2k=20000001T00:00:smos_sec
    // using the fuction timegm in order to convert the time in UTC
    // (maybe not all platforms support this)
    timegm(&y2k);
    if(strftime(smos_date,18*sizeof(smos_date),"%Y%m%dT%T", &y2k)==0){
      cout<<"Error converting EE CFI time format! "<<endl;
      return 0;
    }
    return(1);
  }

  // *********************************************************
  // Funciton to show the USAGE message:
  //
  void ShowUsage(const char *Messages){
    smos::ShowGNUPL();
    mexPrintf("USAGE:\n");
    mexPrintf("\t>> [TSF, SSI]=Process_SMOSxL1C;\n");
    mexPrintf("\t\t File browser will open, select a DBL file and then specify GEOLIMS.\n");
    mexPrintf("\t>> [TSF, SSI]=Process_SMOSxL1C('fname.DBL',GEO_LIMS,OUT_DIR);\n");
    mexPrintf("WHERE:\n");
    mexPrintf("\t\t 'fname.DBL': full path string of DBL file to work with,\n");
    mexPrintf("\t\t GEOLIMS: a 4 elements vector as [LAT1, LON1, LAT2, LON2],\n");
    mexPrintf("\t\t indicating the vertices of a box in a latitude-longitude sapce, being\n");
    mexPrintf("\t\t (LAT1,LON1) & (LAT2,LON2) the latitude and longitude of lower-left and upper-rigth vertex.\n");
    mexPrintf("\t\t OUT_DIR: (optional) string with directory where MAT file 'fname.mat' is stored.\n");
    mexPrintf("\t >> [TSF, SSI]=Process_SMOSxL1C('fname.DBL',GEO_LIMS);\n");
    mexPrintf("\t\t same as usage before, but an output MAT file 'fname.mat' will not be generated.\n");
    mexPrintf("OUTPUTS:\n");
    mexPrintf("\t\t TSF: Structure with all variables for [Temp_Swath_Full] dataset.\n");
    mexPrintf("\t\t SSI: Structure with all variables for [Snapshot_Information]\n");
    mexPrintf("\t\t see documentation for details about variables TSF and SSI.\n");
    mexErrMsgIdAndTxt("SMOSreader:input_output", Messages);
    return;

  }

  // ***************************************************************
  // Funtion to show the GNU Public Lincense notive:
  //
  void ShowGNUPL(){
    mexPrintf("'Process_SMOSxL1C'  Copyright (C) 2016  Pablo Saavedra G.\n");
    mexPrintf("This program comes with ABSOLUTELY NO WARRANTY; for details see GNUPLv3 `LICENCE'.\n");
    mexPrintf("This is free software, and you are welcome to redistribute it\n");
    mexPrintf("under certain conditions; see GNUPLv3 `LICENCE' or <http://www.gnu.org/licenses/> for details.\n");
    return;
  }

} // end of namespace smos

// --------------------------------------
// END OF <Process_SMOSxL1C.cpp>
// ======================================
