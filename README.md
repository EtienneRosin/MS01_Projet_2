# Projet 2 : Problème stationnaire avec grille non-structurée

## Description du Projet
On s'intéresse dans ce projet à la résolution numérique du problème suivant :

### Problème considéré
Soit $\Omega = ]0,a[\times]0,b[$. Pour $f\in L^2(\Omega)$ et $\alpha \in \mathbb{R}$ :

$$
\left|\begin{aligned}
\text{ Trouver } u \in H^1(\Omega) \text{ telle que :} \\
\quad \alpha u - \Delta u &= f,\ & \text{dans } \Omega\\
\quad  \partial_{\boldsymbol{n}} u &= 0,\ & \text{sur } \partial\Omega
\end{aligned}
\right.
$$

Qui est équivalent à la formulation variationnelle :
$$
\left|\begin{aligned}
    &\text{ Trouver } u \in H^1(\Omega) \text{ telle que :} \\
    &\quad \forall v \in H^1(\Omega),\ \alpha\int_\Omega uv\ \text{d} \Omega + \int_\Omega \nabla u \cdot \nabla v\ \text{d}\Omega = \int_\Omega fv\ \text{d} \Omega
    \end{aligned}
    \right.
$$

### Problème discrétisé
On utilise un schéma d'éléments finis $\mathbb{P}_1$ basé sur un maillage de triangle $\mathcal{T}_h$ de $\Omega$ et on obtient le système discrétisé :

$$
\boldsymbol{A} \boldsymbol{u} = \boldsymbol{b}
$$
où : $\boldsymbol{A} = \alpha\boldsymbol{M} + \boldsymbol{K}$ avec $\boldsymbol{M}_{ij} = \int_\Omega w_iw_j\ \text{d} \Omega$, $\boldsymbol{K}_{ij} = \int_\Omega \nabla w_i \cdot \nabla w_j\ \text{d} \Omega$, $\boldsymbol{b} = \boldsymbol{M}\boldsymbol{f}$ avec $\boldsymbol{f}_i = f(M_i)$ et enfin $\boldsymbol{u}_i = u_h(M_i)$.

### Méthodes de résolution
On résoud dans ce projet le système $\boldsymbol{A} \boldsymbol{u} = \boldsymbol{b}$ par méthode de Jacobi et du gradient conjugué dont les implémentations seront parallélisées.

L'énoncé ainsi que les scripts de départ se trouvent aux adresses suivantes : 
- https://ms01.pages.math.cnrs.fr/projet2.pdf
- https://plmlab.math.cnrs.fr/ms01/codes-project-fem-2024/




### Test
```cpp
std::unordered_map<int, int> nodeMultiplicityMap;
for (int nTask = 0; nTask < nbTasks; nTask++) {
    if (nTask != myRank) {
        for (int nExch = 0; nExch < mesh.numNodesToExch(nTask); nExch++) {
            int nodeIndex = mesh.nodesToExch(nTask, nExch);
        nodeMultiplicityMap[nodeIndex]++;
        }
    }
}

  mesh.nodeMultiplicity.resize(mesh.nbOfNodes);
  mesh.nodeMultiplicity.setZero();

for (std::unordered_map<int, int>::const_iterator it = nodeMultiplicityMap.begin(); it != nodeMultiplicityMap.end(); ++it) {
    int nodeIndex = it->first;
    int multiplicity = it->second;
    mesh.nodeMultiplicity(nodeIndex) = multiplicity;
}
```