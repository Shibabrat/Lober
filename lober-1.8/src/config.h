/*********************************************/
/* lober config file                         */
/*********************************************/

///FORMAT FOR OUTPUT IN FILES
#define OUT_FORMAT "%40.30lf"
///FORMAT FOR OUTPUT ON CONSOLE
#define STDOUT_FORMAT "%lf"

///FILENAMES MUST BE SHORTER THAN stringmaxlen
#define stringmaxlen 2048

///PARAMS for the SURFACE GENERATOR
#define lbnx 20
#define lbny 20
#define xmin -2.0
#define xmax -1.0
#define ymin -0.5
#define ymax 0.5
//#define lbnx 500
//#define lbny 500
//#define xmin -1.8
//#define xmax -1.2
//#define ymin -0.25
//#define ymax -0.05

///TURN ON TO REMOVE UNECESSARY PARTS OF THE MANIFOLD
#define cutEnd	1


///DECIMATOR parameters
#define LOBEMINNUMPOINTS -50 //-1 de-activate
#define decimdmin ((xmax-xmin)/(2*lbnx))
