/* File format header used by VisNow for particle codes */
struct ioheader{
  int npart[6];	// The number of particles
  double massarr[6];
  char pad[24];
  int nall[6];
  char fill[256 - 6*4 - 24 - 48 - 6*4];	// Fills to 256 Bytes
};
