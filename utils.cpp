#include "headers.hpp"

//================================================================================
// Compute the local dot product associated to the 2 norm
//================================================================================
double localDotProduct2(const ScaVector& u, const ScaVector& v, const Mesh& mesh){
    double value = 0.0;
    assert(u.size() == v.size() && "Vectors u and v must have the same size.");
    
    for (size_t i = 0; i < u.size(); i++)
    {
        // if ((mesh.nodeMultiplicity(i) !=0) && (mesh.nodeMultiplicity(i) !=1)){
        //     cout << "i = " << mesh.nodeMultiplicity(i) << endl;
        // }
        
        value += u(i) * v(i) / (1 + mesh.nodeMultiplicity(i));
    }
    return value;
}

//================================================================================
// Compute the local squared norm 2 e.g. \sum_i u_i^2
//================================================================================
double localSquaredNorm2(const ScaVector& u, const Mesh& mesh){
    return localDotProduct2(u, u, mesh);
}

//================================================================================
// Compute the local squared norm 2 e.g. u^T M u
//================================================================================
double squaredNormL2(const ScaVector& u, const Problem& pbm){
    return u.dot(pbm.M * u);
}