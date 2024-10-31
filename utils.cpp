#include "headers.hpp"

//================================================================================
// Compute the local squared norm 2 e.g. \sum_i u_i^2
//================================================================================
double localSquaredNorm2(const ScaVector& u, const Mesh& mesh){
    // std::cout << "size : " << mesh.nodeMultiplicity.size() << "\n";
    double value = 0.0;
    for (size_t i = 0; i < u.size(); i++)
    {
        // std::cout << "i = " << i << "\n";
        // if (mesh.nodeMultiplicity(i) != 0){
        //     std::cout << "i = " << i << ", u(i) : " << u(i) << ", Multiplicity : " << mesh.nodeMultiplicity(i) << ", coef : " << u(i) * u(i) / (1 + mesh.nodeMultiplicity(i))<< "\n";
        // }
        // if (mesh.nodeMultiplicity(i) != 0){
        //     std::cout << "i = " << i << ", Multiplicity : " << mesh.nodeMultiplicity(i) << "\n";
        // }
        value += u(i) * u(i) / (1 + mesh.nodeMultiplicity(i));
    }
    return value;
}

//================================================================================
// Compute the local squared norm 2 e.g. u^T M u
//================================================================================
double squaredNormL2(const ScaVector& u, const Problem& pbm){
    double value = u.dot(pbm.M * u);
    return value;
}