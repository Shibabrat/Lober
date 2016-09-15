/************************************************/
/* isInside.h     Pub Interface to -lisInside   */
/* Francois Lekien     [lekien@princeton.edu]   */
/* defines signed_curve and signed_surface that */
/* represent respectively 2D closed polygons    */
/* and 3D polyhedrons.                          */
/* For signed_curve: give the list of points    */
/* (do NOT repeat the first one). For           */
/* signed_surface, list the vertices, then the  */
/* connectivity list (i.e., list of triangles)  */
/************************************************/

char *       __IS_INSIDE__NAME__(void);
unsigned int __IS_INSIDE__VERSION_INTERFACE__(void);
unsigned int __IS_INSIDE__VERSION_REVISION__(void);
unsigned int __IS_INSIDE__VERSION_COMPATIBILITY__(void);
 
class signed_curve {
 public:
  signed_curve(long N, double **x); //constructor from x[2][N]
  /* ARRAY FORMAT:
         x[0][0]=x1
	 x[1][0]=y1
	 x[0][1]=x2
	 x[1][1]=y2
	 ...
	 x[0][N-1]=xN
	 x[1][N-1]=yN
     NOTICE: x1<>xN OR y1<>yN. The curve is closed
             automatically in the constructor.
  ************************************************/
  signed_curve(long N, double *x);  //constructor from x[2*N]
  /* ARRAY FORMAT:
         x[0]=x1
	 x[1]=y1
	 x[2]=x2
	 x[3]=y2
	 ...
	 x[2*N-2]=xN
	 x[2*N-1]=yN
     NOTICE: x1<>xN OR y1<>yN. The curve is closed
             automatically in the constructor.
  ************************************************/
  signed_curve(const char *filename); //constructor from file
  /* FILE FORMAT:
         N
         x1  y1
         x2  y2
         ...
         xN  yN
     NOTICE: x1<>xN OR y1<>yN. The curve is closed
             automatically in the constructor.
  ************************************************/
  ~signed_curve(); //destructor
  int isInside(double x[2]); // 1/0/-1 for inside/boundary/outside
  double d_isInside(double x[2]); // 2PI/PI/0 for inside/boundary/outside
  double area(void); // return area (stored, no computational cost)
  int rotation(void);
 protected:
 private:
  double myArea; //stored area of the surface enclosed (always >=0)
  double area_compute(void); //actual computation of the area (signed)
  long npt; //number of points
  double *pt[2]; 
  int rot; // 1/-1 for counter-clockwise/clockwise
  void make(long npt, double **x);
  void make(long npt, double *x);
};

class signed_surface {
 public:
  signed_surface(long N, long M, double **x, long **tri); //constructor from x[N][2]
  /* ARRAY FORMAT:
         x[0][0]=x1
	 x[0][1]=y1
	 x[1][0]=x2
	 x[1][1]=y2
	 ...
	 x[N-1][0]=xN
	 x[N-1][1]=yN
     NOTICE: x1<>xN OR y1<>yN. The curve is closed
             automatically in the constructor.
  ************************************************/
  signed_surface(long N, long M, double *x, long *tri);  //constructor from x[2*N]
  /* ARRAY FORMAT:
         x[0]=x1
	 x[1]=y1
	 x[2]=x2
	 x[3]=y2
	 ...
	 x[2*N-2]=xN
	 x[2*N-1]=yN
     NOTICE: x1<>xN OR y1<>yN. The curve is closed
             automatically in the constructor.
  ************************************************/
  signed_surface(const char *filename); //constructor from file
  /* FILE FORMAT:
         N   M
         x1  y1
         x2  y2
         ...
         xN  yN
     NOTICE: x1<>xN OR y1<>yN. The curve is closed
             automatically in the constructor.
  ************************************************/
  ~signed_surface(); //destructor
  int isInside(double x[3]); // 1/0/-1 for inside/boundary/outside
  double d_isInside(double x[3]); // 4PI/2PI/0 for inside/boundary/outside
  double volume(void); // return area (stored, no computational cost)
 protected:
 private:
  double myVolume; //stored area of the surface enclosed (always >=0)
  //void signed_surface::normal_vector_tryout(long faceID, int majorPoint, double *v1, double *v2, double *v1xv2_norm, double v1xv2[3], double ve1[3], double ve2[3]);
  void normal_vector_tryout(long faceID, int majorPoint, double *v1, double *v2, double *v1xv2_norm, double v1xv2[3], double ve1[3], double ve2[3]);
  double volume_compute(void); //actual computation of the area (signed)
  double volume_compute_face(long i); //same for a single face
  double d_isInside_face(long i, double x[3], int *bndFlag);//for one face only
  long npt; //number of points
  long ntri; //number of elements
  double *pt[3]; //list of points
  long *tri[3]; //list of elements
  double *normal[3]; //unormalized normal vectors
  double *rotation[3][3]; //rotation matrix (puts normal along 1x
  int *majorPt; //which point (0,1 or 2) on the face is principal
  double *majorV1xV2; // ||v1xv2|| for the major point
  double *V1[2]; //V1 in the new coordinate system  
  double *V2[2]; //V2 in the new coordinate system
  int rot; // 1/-1 for counter-clockwise/clockwise
  void make(long npt, long ntri, double **x, long **tri);
  void make(long npt, long ntri, double *x, long *tri);
  void make_finalize(void);
};
