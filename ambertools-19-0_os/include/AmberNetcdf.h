#ifndef INC_AMBER_NETCDF_H
#define INC_AMBER_NETCDF_H
/*! \file AmberNetcdf.h 
  * \author Daniel R. Roe
  * \date 2010-12-07
  * \brief A C implementation of routines for reading and writing the Amber 
  * Netcdf trajectory and restart formats.
  *
  * Based on Cpptraj implementation.
  * Original implementation of netcdf in Amber by Jon Mongan.
  */

// NOTE: It would be better to allocate memory for single-precision coords
//       upon loading netcdf traj, but since NAB does not really handle pointers
//       the memory must be allocated during every read/write.
// NOTE: Any changes made to the structure below must also be made to 
//       nab_netcdf.h, however due to nab not recognizing the double
//       type all double vars here must be float vars there. NAB float is
//       equivalent to double (see nab/defreal.h etc), but this is why there 
//       have to be two separate definitions of this structure.
/// Hold info for Amber Netcdf trajectory or restart
struct AmberNetcdf {
  double temp0;       // Temperature of current frame (if TempVID!=-1)
  double restartTime; // Simulation time if Amber restart
  int isNCrestart;    // 0 if trajectory, 1 if restart
  int ncid;           // Netcdf ID of the file when open
  int frameDID;       // ID of frame dimension
  int ncframe;        // Number of frames in the file
  int currentFrame;   // Current frame number 
  int atomDID;        // ID of atom dimension
  int ncatom;         // Number of atoms
  int ncatom3;        // Number of coordinates (ncatom * 3)
  int coordVID;       // ID of coordinates variable
  int velocityVID;    // ID of velocities variable
  int cellAngleVID;   // ID of box angle variable
  int cellLengthVID;  // ID of box length variable
  int spatialDID;
  int labelDID;
  int cell_spatialDID;
  int cell_angularDID;
  int spatialVID;
  int timeVID;
  int cell_spatialVID;
  int cell_angularVID;
  int TempVID;
};

// NOTE: To be NAB-useable, any functions added below must also be referenced in
//       nab/nabcode.h and nab/symbol.c
int netcdfDebug(struct AmberNetcdf*);
int netcdfLoad(struct AmberNetcdf*,char *);
int netcdfClose(struct AmberNetcdf*);
int netcdfWriteRestart(char *, int, double *, double *, double *, double, double);
int netcdfCreate(struct AmberNetcdf*,char *, int, int);
int netcdfGetVelocity(struct AmberNetcdf*, int, double *);
int netcdfGetFrame(struct AmberNetcdf*, int, double*, double*);
int netcdfGetNextFrame(struct AmberNetcdf*,double*,double*);
int netcdfWriteFrame(struct AmberNetcdf*, int, double *, double *);
int netcdfWriteNextFrame(struct AmberNetcdf*, double *, double *);
int netcdfInfo(struct AmberNetcdf*);

#endif
