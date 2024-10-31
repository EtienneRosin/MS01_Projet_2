#include "headers.hpp"

extern int myRank;
extern int nbTasks;

//================================================================================
// Solution of the system Au=b with Jacobi
//================================================================================

void jacobi(SpMatrix& A, ScaVector& b, ScaVector& u, Mesh& mesh, double tol, int maxit)
{
  if(myRank == 0)
    cout << "== Jacobi (tol " << tol << "; maxit " << maxit << ")" << endl;
  
  // Compute the solver matrices
  int size = A.rows();
  ScaVector Mdiag(size);
  SpMatrix N(size, size);
  for(int k=0; k<A.outerSize(); ++k){
    for(SpMatrix::InnerIterator it(A,k); it; ++it){
      if(it.row() == it.col())
        Mdiag(it.row()) = it.value();
      else
        N.coeffRef(it.row(), it.col()) = -it.value();
    }
  }
  exchangeAddInterfMPI(Mdiag, mesh);
  
  // Jacobi solver
    // Computing the initial residual
  double initialResidu = 0.0;
  double localSquaredResidu = localSquaredNorm2(b, mesh);
  double globalSquaredResidu = 0.0;
  MPI_Allreduce(
        &localSquaredResidu,
        &globalSquaredResidu,
        1,
        MPI_DOUBLE,
        MPI_SUM,
        MPI_COMM_WORLD);
  initialResidu = std::sqrt(globalSquaredResidu);

  if (myRank == 0)
    {
      cout << "   [0] residual: " << initialResidu << endl;
    }

  double residuNorm = 1e2;
  int it = 1;
  while (residuNorm > tol * initialResidu && it < maxit + 1){
    
    // Compute N*u
    ScaVector Nu = N * u;
    exchangeAddInterfMPI(Nu, mesh);

    // Update field
    for(int i=0; i<size; i++){
      u(i) = 1/Mdiag(i) * (Nu(i) + b(i));
    }
    
    // Update residual and iterator
    ScaVector Au = A * u;
    exchangeAddInterfMPI(Au, mesh);
    ScaVector localResiduVector = b - Au;
    double localSquaredResidu = localSquaredNorm2(localResiduVector, mesh);
    double globalSquaredResidu = 0.0;
    MPI_Allreduce(
        &localSquaredResidu,
        &globalSquaredResidu,
        1,
        MPI_DOUBLE,
        MPI_SUM,
        MPI_COMM_WORLD);

    residuNorm = std::sqrt(globalSquaredResidu);
    

    if ((it % 100) == 0)
    {
      if(myRank == 0)
        cout << "   [" << it << "] residual: " << residuNorm << endl;
    }
    it++;
  }
  
  if(myRank == 0){
    cout << "   -> final iteration: " << it - 1 << endl;
    cout << "   -> final residual: " << residuNorm << endl;
  }
}

//================================================================================
// Solution of the system Au=b with Conjugate Gradient
//================================================================================

void conjgrad(SpMatrix& A, ScaVector& b, ScaVector& u, Mesh& mesh, double tol, int maxit)
{

// TO DO

}