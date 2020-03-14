// Hold info for Amber Netcdf trajectory or restart
struct AmberNetcdf {
  float temp0;        // Temperature of current frame (if TempVID!=-1)
  float restartTime;  // Simulation time if Amber restart
  int isNCrestart;    // 0 if trajectory, 1 if restart
  int ncid;           // Netcdf ID of the file when open
  int frameDID;       // ID of frame dimension
  int ncframe;        // Number of frames in the file
  int currentFrame;   // Current frame number 
  int atomDID;        // ID of atom dimension
  int ncatom;         // Number of atoms
  int ncatom3;        // Number of coordinates (ncatom * 3)
  int coordVID;       // ID of coordinatess variable
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

