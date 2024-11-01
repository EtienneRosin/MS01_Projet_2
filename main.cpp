#include "headers.hpp"

int myRank;
int nbTasks;
// gmsh -2 -part 4 benchmark/mesh.geo
// mpirun -np 4 solver benchmark/mesh.msh
// gmsh benchmark/mesh.msh results/*.msh_*


int main(int argc, char* argv[])
{
  // 1. Initialize MPI ----------------------------------------------
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &nbTasks);


  // 2. Read the mesh, and build lists of nodes for MPI exchanges, local numbering
  if(argc != 2){
    cout << "ERROR: Execution line should be 'solver NameOfMesh.msh'. (" << argc << " input info)." << endl;
    exit(EXIT_FAILURE);
  }
  Mesh mesh;
  readMsh(mesh, argv[1]);
  buildListsNodesMPI(mesh);
  

  // 3. Build problem (vectors and matrices) ------------------------
  double alpha = 1.;
  double a = 2.;
  double b = 1.;

  ScaVector uNum(mesh.nbOfNodes);
  ScaVector uExa(mesh.nbOfNodes);
  ScaVector f(mesh.nbOfNodes);
  for(int i=0; i<mesh.nbOfNodes; ++i){
    double x = mesh.coords(i,0);
    double y = mesh.coords(i,1);
    uNum(i) = 0.;
    uExa(i) = cos(M_PI * x / a) * cos(2. * M_PI * y / b);
    f(i) = (alpha + M_PI * M_PI * (1 / (a * a) + 4 / (b * b))) * uExa(i);
  }

  Problem pbm;
  
  buildProblem(pbm,mesh,alpha,f);
  

  // 4. Solve problem -----------------------------------------------
  double tol = 1e-9; // (Currently useless)
  int maxit = 1000000;
  // jacobi(pbm.A, pbm.b, uNum, mesh, tol, maxit);
  conjgrad(pbm.A, pbm.b, uNum, mesh, tol, maxit);


  // 5. Compute error and export fields -----------------------------
  ScaVector uErr = uNum - uExa;
  exportFieldMsh(uNum, mesh, "solNum", "results/solNum.msh");
  exportFieldMsh(uExa, mesh, "solRef", "results/solExa.msh");
  exportFieldMsh(uErr, mesh, "solErr", "results/solErr.msh");

  double localSquaredError = squaredNormL2(uErr, pbm);
  double globalSquaredError = 0.0;
  MPI_Reduce(
        &localSquaredError,
        &globalSquaredError,
        1,
        MPI_DOUBLE,
        MPI_SUM,
        0, MPI_COMM_WORLD);
  if (myRank == 0){
    std::cout << "error : " << std::sqrt(globalSquaredError) << "\n";
  }
  

  // 6. Finilize MPI ------------------------------------------------
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  
  return 0;
}