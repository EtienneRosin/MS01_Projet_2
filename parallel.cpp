#include "headers.hpp"

extern int myRank;
extern int nbTasks;

//================================================================================
// Build the local numbering and list of nodes for MPI communications
//================================================================================

void buildListsNodesMPI(Mesh& mesh)
{
  if(myRank == 0)
    cout << "== build local numbering and list of nodes for MPI communications" << endl;

  //==== Build mask for nodes belonging to each MPI process (i.e. interior + interface)
  
  IntMatrix maskNodesEachProc(nbTasks, mesh.nbOfNodes);
  for(int nTask=0; nTask<nbTasks; nTask++){
    for(int nGlo=0; nGlo<mesh.nbOfNodes; nGlo++){
      maskNodesEachProc(nTask,nGlo) = 0;
    }
  }
  for(int iTriGlo=0; iTriGlo<mesh.nbOfTri; iTriGlo++){
    int nTask = mesh.triPart(iTriGlo);
    int nGlo0 = mesh.triNodes(iTriGlo,0);
    int nGlo1 = mesh.triNodes(iTriGlo,1);
    int nGlo2 = mesh.triNodes(iTriGlo,2);
    maskNodesEachProc(nTask,nGlo0) = 1;
    maskNodesEachProc(nTask,nGlo1) = 1;
    maskNodesEachProc(nTask,nGlo2) = 1;
  }
  
  //==== Build local numbering for nodes belonging to the current MPI process (i.e. interior + interfaces)
  
  IntVector nodesGloToLoc(mesh.nbOfNodes);
  int nLoc = 0;
  for(int nGlo=0; nGlo<mesh.nbOfNodes; nGlo++){
    if(maskNodesEachProc(myRank, nGlo)){
      nodesGloToLoc(nGlo) = nLoc;  // this node belongs to the current MPI process
      nLoc++;
    }
    else{
      nodesGloToLoc(nGlo) = -1;  // this node does not belong to the current MPI process
    }
  }
  
  //==== Build list with nodes to exchange between the current MPI process and each neighboring MPI process (i.e. interfaces)
  
  mesh.numNodesToExch.resize(nbTasks);
  mesh.nodesToExch.resize(nbTasks,mesh.nbOfNodes);
  for(int nTask=0; nTask<nbTasks; nTask++){
    mesh.numNodesToExch(nTask) = 0;
    if(nTask != myRank){
      int count = 0;
      for(int nGlo=0; nGlo<mesh.nbOfNodes; nGlo++){
        if(maskNodesEachProc(myRank,nGlo) && maskNodesEachProc(nTask,nGlo)){
          mesh.nodesToExch(nTask,count) = nodesGloToLoc(nGlo);
          count++;
        }
      }
      mesh.numNodesToExch(nTask) = count;
      // cout << "   -> task " << myRank << " send/recv " << mesh.numNodesToExch(nTask) << " nodes with task " << nTask << endl;
    }
  }
  if(mesh.numNodesToExch.maxCoeff() > 0){
    mesh.nodesToExch.conservativeResize(nbTasks, mesh.numNodesToExch.maxCoeff());
  }
  
  //==== Build local arrays for nodes/triangles
  
  ScaMatrix coordsMyRank(mesh.nbOfNodes,3);
  IntVector triNumMyRank(mesh.nbOfTri);
  IntMatrix triNodesMyRank(mesh.nbOfTri,3);
  

  nLoc = 0;
  int iLinLoc = 0;
  int iTriLoc = 0;
  for(int nGlo=0; nGlo<mesh.nbOfNodes; nGlo++){
    if(nodesGloToLoc(nGlo) >= 0){
      coordsMyRank(nLoc,0) = mesh.coords(nGlo,0);
      coordsMyRank(nLoc,1) = mesh.coords(nGlo,1);
      coordsMyRank(nLoc,2) = mesh.coords(nGlo,2);
      nLoc++;
    }
  }
  for(int iTriGlo=0; iTriGlo<mesh.nbOfTri; iTriGlo++){
    if(mesh.triPart(iTriGlo) == myRank){
      triNumMyRank(iTriLoc) = mesh.triNum(iTriGlo);
      triNodesMyRank(iTriLoc,0) = nodesGloToLoc(mesh.triNodes(iTriGlo,0));
      triNodesMyRank(iTriLoc,1) = nodesGloToLoc(mesh.triNodes(iTriGlo,1));
      triNodesMyRank(iTriLoc,2) = nodesGloToLoc(mesh.triNodes(iTriGlo,2));
      iTriLoc++;
    }
  }
  
  coordsMyRank.conservativeResize(nLoc,3);
  triNumMyRank.conservativeResize(iTriLoc);
  triNodesMyRank.conservativeResize(iTriLoc,3);
  
  mesh.nbOfNodes = nLoc;
  mesh.nbOfTri = iTriLoc;
  
  mesh.coords = coordsMyRank;
  mesh.triNum = triNumMyRank;
  mesh.triNodes = triNodesMyRank;

    // create the nodeMultiplicity vector 
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
}

//================================================================================
// MPI-parallel exchange/add the interface terms
//================================================================================

void exchangeAddInterfMPI(ScaVector& vec, Mesh& mesh)
{
  
  vector<MPI_Request> requestSnd(nbTasks);
  vector<MPI_Request> requestRcv(nbTasks);
  MPI_Status status;
  vector<ScaVector> bufferSnd(nbTasks);
  vector<ScaVector> bufferRcv(nbTasks);

  for(int nTask=0; nTask<nbTasks; nTask++){
    int numToExch = mesh.numNodesToExch(nTask);
    if(numToExch > 0){
      bufferSnd[nTask].resize(numToExch);
      bufferRcv[nTask].resize(numToExch);
      for(int nExch=0; nExch<numToExch; nExch++)
        bufferSnd[nTask][nExch] = vec(mesh.nodesToExch(nTask,nExch));
      MPI_Isend(&bufferSnd[nTask][0], numToExch, MPI_DOUBLE, nTask, 0, MPI_COMM_WORLD, &requestSnd[nTask]);
      MPI_Irecv(&bufferRcv[nTask][0], numToExch, MPI_DOUBLE, nTask, 0, MPI_COMM_WORLD, &requestRcv[nTask]);
    }
  }
  
  for(int nTask=0; nTask<nbTasks; nTask++){
    int numToExch = mesh.numNodesToExch(nTask);
    if(numToExch > 0){
      MPI_Wait(&requestRcv[nTask], &status);
      for(int nExch=0; nExch<numToExch; nExch++)
        vec(mesh.nodesToExch(nTask,nExch)) += bufferRcv[nTask][nExch];
      MPI_Wait(&requestSnd[nTask], &status);
    }
  }
}