// Translate this file with
//
// g++ -O3 --std=c++11 assignment-2019.c -o assignment
//
// Run it with
//
// ./assignment
//
// There should be a result.pvd file that you can open with Paraview.
// Sometimes, Paraview requires to select the representation "Point Gaussian"
// to see something meaningful.
//
// (C) 2018-2019 Tobias Weinzierl

#include <omp.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <limits>
#include <iomanip>


double t          = 0;
double tFinal     = 0;
double tPlot      = 0;
double tPlotDelta = 0;

int NumberOfBodies = 0;

/**
 * Pointer to pointers. Each pointer in turn points to three coordinates, i.e.
 * each pointer represents one molecule/particle/body.
 * The double asterisk is what makes it a pointer to pointers!
 */
double** x;

/**
 * Equivalent to x storing the velocities.
 * Stores the velocity for each dimension (x,y,z)
 */
double** v;

/**
 * One mass entry per molecule/particle.
 */
double*  mass;

/**
 * Global time step size used.
 */
double   timeStepSize = 0.0;

/**
 * Maximum velocity of all particles.
 */
double   maxV;


/**
 * Minimum distance between two elements.
 */
double   minDx;


/**
 * Set up scenario from the command line.
 *
 * This operation is not to be changed in the assignment.
 */

void setUp(int argc, char** argv) {
  NumberOfBodies = (argc-4) / 7;
  /**
   * Initialising the partciles. x stores coordinates of all particles, 
   * v all the velocities and 
   * mass all the masses
   */
  x    = new double*[NumberOfBodies];
  v    = new double*[NumberOfBodies];
  mass = new double [NumberOfBodies];

  int readArgument = 1;

  tPlotDelta   = std::stof(argv[readArgument]); readArgument++;
  tFinal       = std::stof(argv[readArgument]); readArgument++;
  timeStepSize = std::stof(argv[readArgument]); readArgument++;

  for (int i=0; i<NumberOfBodies; i++) {
    x[i] = new double[3];
    v[i] = new double[3];

    x[i][0] = std::stof(argv[readArgument]); readArgument++;
    x[i][1] = std::stof(argv[readArgument]); readArgument++;
    x[i][2] = std::stof(argv[readArgument]); readArgument++;

    v[i][0] = std::stof(argv[readArgument]); readArgument++;
    v[i][1] = std::stof(argv[readArgument]); readArgument++;
    v[i][2] = std::stof(argv[readArgument]); readArgument++;

    mass[i] = std::stof(argv[readArgument]); readArgument++;

    if (mass[i]<=0.0 ) {
      std::cerr << "invalid mass for body " << i << std::endl;
      exit(-2);
    }
  }

  std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;
  
  if (tPlotDelta<=0.0) {
    std::cout << "plotting switched off" << std::endl;
    tPlot = tFinal + 1.0;
  }
  else {
    std::cout << "plot initial setup plus every " << tPlotDelta << " time units" << std::endl;
    tPlot = 0.0;
  }
}


std::ofstream videoFile;


/**
 * This operation is not to be changed in the assignment.
 */
void openParaviewVideoFile() {
  videoFile.open( "result.pvd" );
  videoFile << "<?xml version=\"1.0\"?>" << std::endl
            << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl
            << "<Collection>";
}


/**
 * This operation is not to be changed in the assignment.
 */
void closeParaviewVideoFile() {
  videoFile << "</Collection>"
            << "</VTKFile>" << std::endl;
}


/**
 * The file format is documented at http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 *
 * This operation is not to be changed in the assignment.
 */
void printParaviewSnapshot() {
  static int counter = -1;
  counter++;
  std::stringstream filename;
  filename << "result-" << counter <<  ".vtp";
  std::ofstream out( filename.str().c_str() );
  out << "<VTKFile type=\"PolyData\" >" << std::endl
      << "<PolyData>" << std::endl
      << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
      << "  <Points>" << std::endl
      << "   <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";


  for (int i=0; i<NumberOfBodies; i++) {
    out << x[i][0]
        << " "
        << x[i][1]
        << " "
        << x[i][2]
        << " ";
  }

  out << "   </DataArray>" << std::endl
      << "  </Points>" << std::endl
      << " </Piece>" << std::endl
      << "</PolyData>" << std::endl
      << "</VTKFile>"  << std::endl;

  videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>" << std::endl;
}


/**
 * This is the only operation you are allowed to change in the assignment.
 */
void updateBody() {

  //  * Instead of creating 'actual' bucket objects and storing values in them, I 
  //  * have instead allocated a bucket value to each particle. This is much more 
  //  * efficient as one can iterate through the particles with ease and with 
  //  * far fewer repititions. 

  // track initial body count for tidying up memory at the end (NumberOfBodies decreases on collision)
  int startBodyCount = NumberOfBodies;

  // double for storing distance
  double distance;

  // diameter of the particles
  double diameter = 0.01;

  // initialising buckets
  int numberOfBuckets = 10;
  int* bucketArray = new int[NumberOfBodies];

  #pragma omp parallel for
  for (int i=0; i<NumberOfBodies; i++){
    // std::cout << "Bucket " << i << " set up on thread: " << omp_get_thread_num() << "\n";
    bucketArray[i] = -1;
  }

  // if first iteration of program, put all particles into first bucket
  if (maxV == 0){
    
    #pragma omp parallel for
    for (int i=0; i<NumberOfBodies; i++){
      // put all bodies in bucket 0
      bucketArray[i] = 0;
    }
  
  }else{


    double vBucket = maxV / numberOfBuckets;

    // loop through all particles    
    #pragma omp parallel for
    for (int i=0; i<NumberOfBodies; i++){
        double particleVel = sqrt(
          (v[i][0] * v[i][0]) +
          (v[i][1] * v[i][1]) +
          (v[i][2] * v[i][2])
        );
      // loop through all buckets, check if particle goes into that
      for (int j=0; j<numberOfBuckets; j++){
      
        if (((particleVel >= j*vBucket) && (particleVel < (j+1)*vBucket)) || ((j == numberOfBuckets -1) && particleVel >= maxV)){
         
            bucketArray[i] = j;
          
        }
        
      }
    }
  }



  // variable to track the minimum distance between any two particles
  minDx  = std::numeric_limits<double>::max();
  
  // matrix that stores the forces acting on all particles in the scene
  double** forceMatrix;
  forceMatrix = new double*[NumberOfBodies];
  
  #pragma parallel for
  for (int i=0; i<NumberOfBodies; i++) {
	  // one entry for each dimension (x,y,z)
	  forceMatrix[i] = new double[3]{0.0, 0.0, 0.0};
  }


  // loop through all particles

  // cannot parallelise first loop - too many dependencies
  for (int i = 0; i < NumberOfBodies; i++) {
    
    // set up looping based on buckets
    int timeStepsToTake = pow(2, (bucketArray[i]));
    int timeStepsRemaining = timeStepsToTake;
    double miniTimeStepSize = timeStepSize / timeStepsToTake;
  
    // loop through all particles that i hasn't been compared to yet, calculate forces
    #pragma omp parallel for reduction (min:minDx)
    for (int j = i + 1;  j < NumberOfBodies; j++) {

      // Calculate the distance from particle i to particle j
      distance = std::sqrt(
            (x[i][0]-x[j][0]) * (x[i][0]-x[j][0]) +
            (x[i][1]-x[j][1]) * (x[i][1]-x[j][1]) +
            (x[i][2]-x[j][2]) * (x[i][2]-x[j][2])
          );
  
      // track minimum distance
      minDx = std::min( minDx,distance );

      
      // apply forces to each dimension k
      for (int k = 0; k < 3; k++) {
      
        double force = (x[j][k]-x[i][k]) * mass[j] * mass[i] / distance / distance / distance ;
        // critical to avoid read/write clashes
        #pragma omp critical
        {
          forceMatrix[i][k] += force;
          forceMatrix[j][k] -= force;
        }
        
      }
    }

    double absV = 0;

    for (int k=0; k<3; k++){
      v[i][k] = v[i][k] + timeStepSize * forceMatrix[i][k] / mass[i];
      absV += v[i][k] * v[i][k];
    }
    maxV = std::max(maxV, std::sqrt(absV));

    // perform mini time steps, update position and check for collision
    while (timeStepsRemaining > 0){
      // update particle position
      for (int k=0; k<3; k++){
        x[i][k] += miniTimeStepSize * v[i][k];
      }

      // check for collision
      // can't parallelise collisions, way too many dependencies
      for (int j=0; j<NumberOfBodies; j++){
        // if i and j are not the same particle
        if (i != j){          
          distance = std::sqrt(
            (x[i][0]-x[j][0]) * (x[i][0]-x[j][0]) +
            (x[i][1]-x[j][1]) * (x[i][1]-x[j][1]) +
            (x[i][2]-x[j][2]) * (x[i][2]-x[j][2])
          );


          // if a collision occurs
          if (distance < diameter){
            // if a collision occurs then allocate ONLY one thread to deal with it
          
            // update new particle's velocity and position halfway between each
            for (int k = 0; k < 3; k++) {
              x[i][k] = (x[i][k] + x[j][k]) / 2;
              v[i][k] = (mass[i] / (mass[i] + mass[j])) * v[i][k] + (mass[j] / (mass[i] + mass[j])) * v[j][k];
            }
            // update new particle's mass
            mass[i] += mass[j];

            // remove j
            // can't parallelise
            for (int k=j; k<NumberOfBodies-1; k++){
              x[k] = x[k+1];
              v[k] = v[k+1];
              mass[k] = mass[k+1];
              forceMatrix[k] = forceMatrix[k+1];
              bucketArray[k] = bucketArray[k+1];
            }
            if (j < i){
              i -= 1;
            }
            j -= 1;
            NumberOfBodies -= 1;
          
          }
        }
      }    
      timeStepsRemaining -= 1;
    }
  }

  // check to see if only one particle left
  if (NumberOfBodies == 1){
    // terminate
    t = tFinal;
    std::cout << "\n\n\nFinal particle coordinates: " << x[0][0] << ", " << x[0][1] << ", " << x[0][2] << "\n\n\n";
  }

  t += timeStepSize;

  for (int i=0; i<startBodyCount; i++){
    delete[] forceMatrix[i];
  }
  delete[] bucketArray;
  delete[] forceMatrix;
  
}

/**
 * Main routine.
 *
 * Not to be changed in assignment.
 */
int main(int argc, char** argv) {

  if (argc==1) {
    std::cerr << "usage: " + std::string(argv[0]) + " snapshot final-time dt objects" << std::endl
              << "  snapshot        interval after how many time units to plot. Use 0 to switch off plotting" << std::endl
              << "  final-time      simulated time (greater 0)" << std::endl
              << "  dt              time step size (greater 0)" << std::endl
              << std::endl
              << "Examples:" << std::endl
              << "0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0 \t One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
              << "0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0     0.0 1.0 0.0  1.0 0.0 0.0  1.0  \t One spiralling around the other one" << std::endl
              << "0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0 \t Three body setup from first lecture" << std::endl
              << "0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0     2.0 1.0 0.0  0.0 0.0 0.0  1.0     2.0 0.0 1.0  0.0 0.0 0.0  1.0 \t Five body setup" << std::endl
              << std::endl
              << "In this naive code, only the first body moves" << std::endl;

    return -1;
  }
  else if ( (argc-4)%7!=0 ) {
    std::cerr << "error in arguments: each planet is given by seven entries (position, velocity, mass)" << std::endl;
    std::cerr << "got " << argc << " arguments (three of them are reserved)" << std::endl;
    std::cerr << "run without arguments for usage instruction" << std::endl;
    return -2;
  }

  std::cout << std::setprecision(15);

  setUp(argc,argv);

  openParaviewVideoFile();

  int snapshotCounter = 0;
  if (t > tPlot) {
    printParaviewSnapshot();
    std::cout << "plotted initial setup" << std::endl;
    tPlot = tPlotDelta;
  }

  int timeStepCounter = 0;
  while (t<=tFinal) {
    updateBody();
    timeStepCounter++;
    if (t >= tPlot) {
      printParaviewSnapshot();
      std::cout << "plot next snapshot"
    		    << ",\t time step=" << timeStepCounter
    		    << ",\t t="         << t
				<< ",\t dt="        << timeStepSize
				<< ",\t v_max="     << maxV
				<< ",\t dx_min="    << minDx
				<< std::endl;

      tPlot += tPlotDelta;
    }
  }

  closeParaviewVideoFile();

  return 0;
}
