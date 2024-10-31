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
  ScaVector Au = A * u;
  exchangeAddInterfMPI(Au, mesh);
  ScaVector initialLocalResidu = b - Au;
  double localSquaredResidu = localSquaredNorm2(initialLocalResidu, mesh);
  double squaredResidu = 0.0;
  MPI_Allreduce(&localSquaredResidu, &squaredResidu, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  initialResidu = std::sqrt(squaredResidu);

  if (myRank == 0)
  {
    cout << "   [0] residual: " << initialResidu << endl;
  }

  double residuNorm = initialResidu;
  int it = 1;
  while (residuNorm > tol * initialResidu && it < maxit + 1){
    
    // Compute N*u
    ScaVector Nu = N * u;
    exchangeAddInterfMPI(Nu, mesh);

    // Update field
    for(int i=0; i<u.size(); i++){
      u(i) = 1/Mdiag(i) * (Nu(i) + b(i));
    }
    
    // Update residual and iterator
    Au = A * u;
    exchangeAddInterfMPI(Au, mesh);
    ScaVector localResiduVector = b - Au;
    double localSquaredResidu = localSquaredNorm2(localResiduVector, mesh);
    double squaredResidu = 0.0;
    MPI_Allreduce(
        &localSquaredResidu,
        &squaredResidu,
        1,
        MPI_DOUBLE,
        MPI_SUM,
        MPI_COMM_WORLD);

    residuNorm = std::sqrt(squaredResidu);
    

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
    if (myRank == 0)
        cout << "== Conjugate Gradient (tol " << tol << "; maxit " << maxit << ")" << endl;

    // Initial residual calculation
    ScaVector Au_0 = A * u;
    exchangeAddInterfMPI(Au_0, mesh);
    ScaVector r_l = b - Au_0;

    double localSquaredResidu = localSquaredNorm2(r_l, mesh);
    double squaredResidu = 0.0;
    MPI_Allreduce(&localSquaredResidu, &squaredResidu, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double initialResidu = std::sqrt(squaredResidu);
    double residuNorm = initialResidu;
    if (std::isnan(initialResidu) || initialResidu == 0.0) {
        cout << "Initial residual is NaN or zero, check matrix A and vector b." << endl;
        return;
    }


    ScaVector p_l = r_l;  // initial direction
    int it = 1;

    ScaVector Ap_l(u.size());
    double alpha_l = 0.0, beta_l = 0.0, A_p_l_dot_p_l = 0.0, local_A_p_l_dot_p_l = 0.0;
    double r_l_dot_p_l = 0.0, local_r_l_dot_p_l = 0.0, local_A_r_l_dot_p_l = 0.0, A_r_l_dot_p_l = 0.0;
    if (myRank == 0)
      cout << "   [0] residual: " << initialResidu << endl;

    // Conjugate Gradient iteration loop
    while (residuNorm > tol * initialResidu && it <= maxit)
    {
        Ap_l = A * p_l;
        exchangeAddInterfMPI(Ap_l, mesh);
        
        local_A_p_l_dot_p_l = localDotProduct2(Ap_l, p_l, mesh);
        MPI_Allreduce(&local_A_p_l_dot_p_l, &A_p_l_dot_p_l, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        // Compute the step
        local_r_l_dot_p_l = localDotProduct2(r_l, p_l, mesh);
        MPI_Allreduce(&local_r_l_dot_p_l, &r_l_dot_p_l, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        alpha_l = r_l_dot_p_l / A_p_l_dot_p_l;


        // Update solution and residual
        u += alpha_l * p_l;
        r_l -= alpha_l * Ap_l;

        // Update direction with beta
        local_A_r_l_dot_p_l = localDotProduct2(A * r_l, p_l, mesh);
        MPI_Allreduce(&local_A_r_l_dot_p_l, &A_r_l_dot_p_l, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        beta_l = -A_r_l_dot_p_l / A_p_l_dot_p_l;
        p_l = r_l + beta_l * p_l;

        // Compute new residual norm
        localSquaredResidu = localSquaredNorm2(r_l, mesh);
        MPI_Allreduce(&localSquaredResidu, &squaredResidu, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        residuNorm = std::sqrt(squaredResidu);

        if (std::isnan(residuNorm)) {
            if (myRank == 0){
              cout << "Residual norm is NaN; stopping iteration." << endl;
            };
            break;
        }

        if ((it % 100) == 0 && myRank == 0)  // Display every 10 iterations
            cout << "   [" << it << "] residual: " << residuNorm << endl;

        it++;
    }
  
    if (myRank == 0) {
        cout << "   -> final iteration: " << it - 1 << endl;
        cout << "   -> final residual: " << residuNorm << endl;
    }
}