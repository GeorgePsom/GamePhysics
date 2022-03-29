#ifndef SCENE_HEADER_FILE
#define SCENE_HEADER_FILE

#include <vector>
#include <queue>
#include <fstream>
#include <igl/bounding_box.h>
#include <igl/readOFF.h>
#include <igl/per_vertex_normals.h>
#include <igl/edge_topology.h>
#include <igl/diag.h>
#include <igl/readMESH.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include "constraints.h"
#include "auxfunctions.h"
#include "ccd.h"

using namespace Eigen;
using namespace std;

void support(const void *_obj, const ccd_vec3_t *_d, ccd_vec3_t *_p);
void stub_dir(const void *obj1, const void *obj2, ccd_vec3_t *dir);
void center(const void *_obj, ccd_vec3_t *dir);

//the class the contains each individual rigid objects and their functionality
class Mesh{
public:
  
  //position
  VectorXd origPositions;     //3|V|x1 original vertex positions in xyzxyz format - never change this!
  VectorXd currPositions;     //3|V|x1 current vertex positions in xyzxyz format
  
  //kinematics
  bool isFixed;               //is the object immobile (infinite mass)
  VectorXd currVelocities;    //3|V|x1 velocities per coordinate in xyzxyz format.
  
  MatrixXi T;                 //|T|x4 tetrahdra
  MatrixXi F;                 //|F|x3 boundary faces
  VectorXd invMasses;         //|V|x1 inverse masses of vertices, computed in the beginning as 1.0/(density * vertex voronoi area)
  VectorXd voronoiVolumes;    //|V|x1 the voronoi volume of vertices
  VectorXd tetVolumes;        //|T|x1 tetrahedra volumes
  int globalOffset;           //the global index offset of the of opositions/velocities/impulses from the beginning of the global coordinates array in the containing scene class
  
  VectorXi boundTets;  //just the boundary tets, for collision
  
  double youngModulus, poissonRatio, density, alpha, beta;
  
  SparseMatrix<double> A, K, M, D;   //The soft-body matrices
  
  SimplicialLLT<SparseMatrix<double>>* ASolver;   //the solver for the left-hand side matrix constructed for FEM
  
  ~Mesh(){if (ASolver!=NULL) delete ASolver;}
  
  //Quick-reject checking collision between mesh bounding boxes.
  bool isBoxCollide(const Mesh& m2){
    RowVector3d XMin1=RowVector3d::Constant(3276700.0);
    RowVector3d XMax1=RowVector3d::Constant(-3276700.0);
    RowVector3d XMin2=RowVector3d::Constant(3276700.0);
    RowVector3d XMax2=RowVector3d::Constant(-3276700.0);
    for (int i=0;i<origPositions.size();i+=3){
      XMin1=XMin1.array().min(currPositions.segment(i,3).array().transpose());
      XMax1=XMax1.array().max(currPositions.segment(i,3).array().transpose());
    }
    for (int i=0;i<m2.origPositions.size();i+=3){
      XMin2=XMin2.array().min(m2.currPositions.segment(i,3).array().transpose());
      XMax2=XMax2.array().max(m2.currPositions.segment(i,3).array().transpose());
    }
    
    /*double rmax1=vertexSphereRadii.maxCoeff();
     double rmax2=m2.vertexSphereRadii.maxCoeff();
     XMin1.array()-=rmax1;
     XMax1.array()+=rmax1;
     XMin2.array()-=rmax2;
     XMax2.array()+=rmax2;*/
    
    //checking all axes for non-intersection of the dimensional interval
    for (int i=0;i<3;i++)
      if ((XMax1(i)<XMin2(i))||(XMax2(i)<XMin1(i)))
        return false;
    
    return true;  //all dimensional intervals are overlapping = possible intersection
  }
  
  bool isNeighborTets(const RowVector4i& tet1, const RowVector4i& tet2){
    for (int i=0;i<4;i++)
      for (int j=0;j<4;j++)
        if (tet1(i)==tet2(j)) //shared vertex
          return true;
    
    return false;
  }
  
  
  //this function creates all collision constraints between vertices of the two meshes
  void createCollisionConstraints(const Mesh& m, const bool sameMesh, const double timeStep, const double CRCoeff, vector<Constraint>& activeConstraints){
    
    
    //collision between bounding boxes
    if (!isBoxCollide(m))
      return;
    
    if ((isFixed && m.isFixed))  //collision does nothing
      return;
    
    //creating tet spheres
    /*MatrixXd c1(T.rows(), 3);
     MatrixXd c2(m.T.rows(), 3);
     VectorXd r1(T.rows());
     VectorXd r2(m.T.rows());*/
    
    MatrixXd maxs1(boundTets.rows(),3);
    MatrixXd mins1(boundTets.rows(),3);
    MatrixXd maxs2(m.boundTets.rows(),3);
    MatrixXd mins2(m.boundTets.rows(),3);
    
    for (int i = 0; i < boundTets.size(); i++) {
      MatrixXd tet1(4, 3); tet1 << currPositions.segment(3 * T(boundTets(i), 0), 3).transpose(),
      currPositions.segment(3 * T(boundTets(i), 1), 3).transpose(),
      currPositions.segment(3 * T(boundTets(i), 2), 3).transpose(),
      currPositions.segment(3 * T(boundTets(i), 3), 3).transpose();
      
      //c1.row(i) = tet1.colwise().mean();
      //r1(i) = ((c1.row(i).replicate(4, 1) - tet1).rowwise().norm()).maxCoeff();
      mins1.row(i)=tet1.colwise().minCoeff();
      maxs1.row(i)=tet1.colwise().maxCoeff();
      
    }
    
    for (int i = 0; i < m.boundTets.size(); i++) {
      
      MatrixXd tet2(4, 3); tet2 << m.currPositions.segment(3 * m.T(m.boundTets(i), 0), 3).transpose(),
      m.currPositions.segment(3 * m.T(m.boundTets(i), 1), 3).transpose(),
      m.currPositions.segment(3 * m.T(m.boundTets(i), 2), 3).transpose(),
      m.currPositions.segment(3 * m.T(m.boundTets(i), 3), 3).transpose();
      
      //c2.row(i) = tet2.colwise().mean();
      //r2(i) = ((c2.row(i).replicate(4, 1) - tet2).rowwise().norm()).maxCoeff();
      mins2.row(i)=tet2.colwise().minCoeff();
      maxs2.row(i)=tet2.colwise().maxCoeff();
      
    }
    
    //checking collision between every tetrahedrons
    std::list<Constraint> collisionConstraints;
    for (int i=0;i<boundTets.size();i++){
      for (int j=0;j<m.boundTets.size();j++){
        
        //not checking for collisions between tetrahedra neighboring to the same vertices
        if (sameMesh)
          if (isNeighborTets(T.row(boundTets(i)), m.T.row(m.boundTets(j))))
            continue;  //not creating collisions between neighboring tets
        
        bool overlap=true;
        for (int k=0;k<3;k++)
          if ((maxs1(i,k)<mins2(j,k))||(maxs2(j,k)<mins1(i,k)))
            overlap=false;
        
        if (!overlap)
          continue;
        
        VectorXi globalCollisionIndices(24);
        VectorXd globalInvMasses(24);
        for (int t=0;t<4;t++){
          globalCollisionIndices.segment(3*t,3)<<globalOffset+3*(T(boundTets(i),t)), globalOffset+3*(T(boundTets(i),t))+1, globalOffset+3*(T(boundTets(i),t))+2;
          globalInvMasses.segment(3*t,3)<<invMasses(T(boundTets(i),t)), invMasses(T(boundTets(i),t)),invMasses(T(boundTets(i),t));
          globalCollisionIndices.segment(12+3*t,3)<<m.globalOffset+3*m.T(m.boundTets(j),t), m.globalOffset+3*m.T(m.boundTets(j),t)+1, m.globalOffset+3*m.T(m.boundTets(j),t)+2;
          globalInvMasses.segment(12+3*t,3)<<m.invMasses(m.T(m.boundTets(j),t)), m.invMasses(m.T(m.boundTets(j),t)),m.invMasses(m.T(m.boundTets(j),t));
        }
        
        ccd_t ccd;
        CCD_INIT(&ccd);
        ccd.support1       = support; // support function for first object
        ccd.support2       = support; // support function for second object
        ccd.center1         =center;
        ccd.center2         =center;
        
        ccd.first_dir       = stub_dir;
        ccd.max_iterations = 100;     // maximal number of iterations
        
        MatrixXd tet1(4, 3); tet1 << currPositions.segment(3 * T(boundTets(i), 0), 3).transpose(),
        currPositions.segment(3 * T(boundTets(i), 1), 3).transpose(),
        currPositions.segment(3 * T(boundTets(i), 2), 3).transpose(),
        currPositions.segment(3 * T(boundTets(i), 3), 3).transpose();
        
        MatrixXd tet2(4, 3); tet2 << m.currPositions.segment(3 * m.T(m.boundTets(j), 0), 3).transpose(),
        m.currPositions.segment(3 * m.T(m.boundTets(j), 1), 3).transpose(),
        m.currPositions.segment(3 * m.T(m.boundTets(j), 2), 3).transpose(),
        m.currPositions.segment(3 * m.T(m.boundTets(j), 3), 3).transpose();
        
        void* obj1=(void*)&tet1;
        void* obj2=(void*)&tet2;
        
        ccd_real_t _depth;
        ccd_vec3_t dir, pos;
        
        int nonintersect = ccdMPRPenetration(obj1, obj2, &ccd, &_depth, &dir, &pos);
        
        if (nonintersect)
          continue;
        
        Vector3d intNormal, intPosition;
        double depth;
        for (int k=0;k<3;k++){
          intNormal(k)=dir.v[k];
          intPosition(k)=pos.v[k];
        }
        
        depth =_depth;
        intPosition-=depth*intNormal/2.0;
        
        Vector3d p1=intPosition+depth*intNormal;
        Vector3d p2=intPosition;
      
        //getting barycentric coordinates of each point
        
        MatrixXd PMat1(4,4); PMat1<<1.0,currPositions.segment(3*T(boundTets(i),0),3).transpose(),
        1.0,currPositions.segment(3*T(boundTets(i),1),3).transpose(),
        1.0,currPositions.segment(3*T(boundTets(i),2),3).transpose(),
        1.0,currPositions.segment(3*T(boundTets(i),3),3).transpose();
        PMat1.transposeInPlace();
        
        Vector4d rhs1; rhs1<<1.0,p1;
        
        Vector4d B1=PMat1.inverse()*rhs1;
        
        MatrixXd PMat2(4,4); PMat2<<1.0,m.currPositions.segment(3*m.T(m.boundTets(j),0),3).transpose(),
        1.0,m.currPositions.segment(3*m.T(m.boundTets(j),1),3).transpose(),
        1.0,m.currPositions.segment(3*m.T(m.boundTets(j),2),3).transpose(),
        1.0,m.currPositions.segment(3*m.T(m.boundTets(j),3),3).transpose();
        PMat2.transposeInPlace();
        
        Vector4d rhs2; rhs2<<1.0,p2;
        
        Vector4d B2=PMat2.inverse()*rhs2;
        
        //cout<<"B1: "<<B1<<endl;
        //cout<<"B2: "<<B2<<endl;
        
        //Matrix that encodes the vector between interpenetration points by the c
        MatrixXd v2cMat1(3,12); v2cMat1.setZero();
        for (int k=0;k<3;k++){
          v2cMat1(k,k)=B1(0);
          v2cMat1(k,3+k)=B1(1);
          v2cMat1(k,6+k)=B1(2);
          v2cMat1(k,9+k)=B1(3);
        }
        
        MatrixXd v2cMat2(3,12); v2cMat2.setZero();
        for (int k=0;k<3;k++){
          v2cMat2(k,k)=B2(0);
          v2cMat2(k,3+k)=B2(1);
          v2cMat2(k,6+k)=B2(2);
          v2cMat2(k,9+k)=B2(3);
        }
        
        MatrixXd v2dMat(3,24); v2dMat<<-v2cMat1, v2cMat2;
        VectorXd constVector=intNormal.transpose()*v2dMat;
        
        //cout<<"intNormal: "<<intNormal<<endl;
        //cout<<"n*(p2-p1): "<<intNormal.dot(p2-p1)<<endl;
        collisionConstraints.push_back(Constraint(COLLISION, INEQUALITY,  globalCollisionIndices, globalInvMasses, constVector, 0, CRCoeff));
        
        //i=10000000;
        //break;
        
      }
    }
    
    activeConstraints.insert(activeConstraints.end(), collisionConstraints.begin(), collisionConstraints.end());
  }
  
  
  
  void createGlobalMatrices(const double timeStep, const double _alpha, const double _beta)
  {
   
    /*************************TODO: create the M, D, K matrices from alpha, beta, poisson ratio, and Young's modulus as learnt in class. Afterward create the matrix "A" with the given timeStep that is the left hand side of the entire system.
     *********/

    //M
    //  M is the 3|V|x3|V| mass matrix, which is a diagonal matrix with mass value per coordinate
    //      in a(t) (that means repeating the mass of each vertex three times).

    int Vnum = invMasses.size();
    
    MatrixXd M = MatrixXd::Zero(Vnum, Vnum);

    for (int i = 0; i < Vnum; i++)
    {
        M(i, i) = 1 / invMasses(i);
    }

    //D is the 3 jV j  3 jV j damping matrix which slows down(damps) the movement.It is an
    //     abstraction of physical damping, which is mostly here for the sake of stabilizing the discrete
    //     simulation.

    //MatrixXd G_e = MatrixXd::Zero(Vnum, Vnum);

    //K

    
    A=M+D*timeStep+K*(timeStep*timeStep);
    
    //Should currently fail since A is empty
    if (ASolver==NULL)
      ASolver=new SimplicialLLT<SparseMatrix<double>>();
    ASolver->analyzePattern(A);
    ASolver->factorize(A);
    
  }
  
  //returns center of mass
  Vector3d initializeVolumesAndMasses()
  {
    //TODO: compute tet volumes and allocate to vertices
    tetVolumes.conservativeResize(T.rows());
    voronoiVolumes.conservativeResize(origPositions.size()/3);
    voronoiVolumes.setZero();
    invMasses.conservativeResize(origPositions.size()/3);
    Vector3d COM; COM.setZero();
    for (int i=0;i<T.rows();i++){
      Vector3d e01=origPositions.segment(3*T(i,1),3)-origPositions.segment(3*T(i,0),3);
      Vector3d e02=origPositions.segment(3*T(i,2),3)-origPositions.segment(3*T(i,0),3);
      Vector3d e03=origPositions.segment(3*T(i,3),3)-origPositions.segment(3*T(i,0),3);
      Vector3d tetCentroid=(origPositions.segment(3*T(i,0),3)+origPositions.segment(3*T(i,1),3)+origPositions.segment(3*T(i,2),3)+origPositions.segment(3*T(i,3),3))/4.0;
      tetVolumes(i)=std::abs(e01.dot(e02.cross(e03)))/6.0;
      for (int j=0;j<4;j++)
        voronoiVolumes(T(i,j))+=tetVolumes(i)/4.0;
      
      COM+=tetVolumes(i)*tetCentroid;
    }
    
    COM.array()/=tetVolumes.sum();
    for (int i=0;i<origPositions.size()/3;i++)
      invMasses(i)=1.0/(voronoiVolumes(i)*density);
    
    return COM;
    
  }
  
  //performing the integration step of the soft body.
  void integrateVelocity(double timeStep){
    
    if (isFixed)
      return;
    
    /****************TODO: construct rhs (right-hand side) and use ASolver->solve(rhs) to solve for velocities********/
    
    VectorXd rhs = VectorXd::Zero(currVelocities.size());  //REMOVE THIS! it's a stub
    currVelocities=ASolver->solve(rhs);
  }
  
  //Update the current position with the integrated velocity
  void integratePosition(double timeStep){
    if (isFixed)
      return;  //a fixed object is immobile
    
    currPositions+=currVelocities*timeStep;
    //cout<<"currPositions: "<<currPositions<<endl;
  }
  
  //the full integration for the time step (velocity + position)
  void integrate(double timeStep){
    integrateVelocity(timeStep);
    integratePosition(timeStep);
  }
  
  
  Mesh(const VectorXd& _origPositions, const MatrixXi& boundF, const MatrixXi& _T, const int _globalOffset, const double _youngModulus, const double _poissonRatio, const double _density, const bool _isFixed, const RowVector3d& userCOM, const RowVector4d& userOrientation){
    origPositions=_origPositions;
    //cout<<"original origPositions: "<<origPositions<<endl;
    T=_T;
    F=boundF;
    isFixed=_isFixed;
    globalOffset=_globalOffset;
    density=_density;
    poissonRatio=_poissonRatio;
    youngModulus=_youngModulus;
    currVelocities=VectorXd::Zero(origPositions.rows());
    
    VectorXd naturalCOM=initializeVolumesAndMasses();
    //cout<<"naturalCOM: "<<naturalCOM<<endl;
    
    
    origPositions-= naturalCOM.replicate(origPositions.rows()/3,1);  //removing the natural COM of the OFF file (natural COM is never used again)
    //cout<<"after natrualCOM origPositions: "<<origPositions<<endl;
    
    for (int i=0;i<origPositions.size();i+=3)
      origPositions.segment(i,3)<<(QRot(origPositions.segment(i,3).transpose(), userOrientation)+userCOM).transpose();
    
    currPositions=origPositions;
    
    if (isFixed)
      invMasses.setZero();
    
    //finding boundary tets
    VectorXi boundVMask(origPositions.rows()/3);
    boundVMask.setZero();
    for (int i=0;i<boundF.rows();i++)
      for (int j=0;j<3;j++)
        boundVMask(boundF(i,j))=1;
    
    cout<<"boundVMask.sum(): "<<boundVMask.sum()<<endl;
    
    vector<int> boundTList;
    for (int i=0;i<T.rows();i++){
      int incidence=0;
      for (int j=0;j<4;j++)
        incidence+=boundVMask(T(i,j));
      if (incidence>2)
        boundTList.push_back(i);
    }
    
    boundTets.resize(boundTList.size());
    for (int i=0;i<boundTets.size();i++)
      boundTets(i)=boundTList[i];
    
    ASolver=NULL;
  }
  
};


//This class contains the entire scene operations, and the engine time loop.
class Scene{
public:
  double currTime;
  
  VectorXd globalPositions;   //3*|V| all positions
  VectorXd globalVelocities;  //3*|V| all velocities
  VectorXd globalInvMasses;   //3*|V| all inverse masses  (NOTE: the invMasses in the Mesh class is |v| (one per vertex)!
  MatrixXi globalT;           //|T|x4 tetraheda in global index
  
  vector<Mesh> meshes;
  
  vector<Constraint> userConstraints;   //provided from the scene
  vector<Constraint> barrierConstraints;  //provided by the platform
  
  //updates from global values back into mesh values
  void global2Mesh(){
    for (int i=0;i<meshes.size();i++){
      meshes[i].currPositions<<globalPositions.segment(meshes[i].globalOffset, meshes[i].currPositions.size());
      meshes[i].currVelocities<<globalVelocities.segment(meshes[i].globalOffset, meshes[i].currVelocities.size());
    }
  }
  
  //update from mesh current values into global values
  void mesh2global(){
    for (int i=0;i<meshes.size();i++){
      globalPositions.segment(meshes[i].globalOffset, meshes[i].currPositions.size())<<meshes[i].currPositions;
      globalVelocities.segment(meshes[i].globalOffset, meshes[i].currVelocities.size())<< meshes[i].currVelocities;
    }
  }
  
  
  //This should be called whenever the timestep changes
  void initScene(double timeStep, const double alpha, const double beta){
    
    for (int i=0;i<meshes.size();i++){
      if (!meshes[i].isFixed)
        meshes[i].createGlobalMatrices(timeStep, alpha, beta);
    }
    
    mesh2global();
    
  }
  
  /*********************************************************************
   This function handles a single time step
   1. Integrating velocities and position from forces and previous impulses
   2. detecting collisions and generating collision constraints
   3. Resolving collision constraints by updating velocities and positions
   *********************************************************************/
  
  void updateScene(double timeStep, double CRCoeff, const double tolerance, const int maxIterations){
    
    /*******************1. Integrating velocity and position from external and internal forces************************************/
    for (int i=0;i<meshes.size();i++)
      meshes[i].integrate(timeStep);
    
    mesh2global();
    
    
    /*******************2. Creating and Aggregating collision constraints************************************/
    
    vector<Constraint> activeConstraints;
            
    //collision constraints
    for (int i=0;i<meshes.size();i++)
      for (int j=i+1;j<meshes.size(); j++)
        meshes[i].createCollisionConstraints(meshes[j], i==j, timeStep, CRCoeff, activeConstraints);
    
    /*************3. Resolving all collision constraints sequentially*******************/
    int currIteration=0;
    for (int currConstIndex = 0;currConstIndex<activeConstraints.size();currConstIndex++){
      
      Constraint currConstraint=activeConstraints[currConstIndex];
      
      VectorXd constraintPositions(currConstraint.globalIndices.size());
      VectorXd constraintVelocities(currConstraint.globalIndices.size());
      for (int i=0;i<currConstraint.globalIndices.size();i++){
        constraintPositions(i)=globalPositions(currConstraint.globalIndices(i));
        constraintVelocities(i)=globalVelocities(currConstraint.globalIndices(i));
      }
      
      //generating impulses
      VectorXd generatedImpulses;
      currConstraint.resolveVelocityConstraint(constraintPositions, constraintVelocities, generatedImpulses, tolerance);
      
      //correcting velocities according to impulses
      for (int i=0;i<currConstraint.globalIndices.size();i++){
        int currIndex=currConstraint.globalIndices(i);
        globalVelocities(currIndex)+=globalInvMasses(currIndex)*generatedImpulses(i);
      }
    }
    
    global2Mesh();
    
    /*******************4. Solving for position drift************************************/
    
    mesh2global();
    
    for (int currConstIndex = 0;currConstIndex<activeConstraints.size();currConstIndex++){
      
      Constraint currConstraint=activeConstraints[currConstIndex];
      
      VectorXd constraintPositions(currConstraint.globalIndices.size());
      VectorXd constraintVelocities(currConstraint.globalIndices.size());
      for (int i=0;i<currConstraint.globalIndices.size();i++){
        constraintPositions(i)=globalPositions(currConstraint.globalIndices(i));
        constraintVelocities(i)=globalVelocities(currConstraint.globalIndices(i));
      }
      
      //generating impulses
      VectorXd generatedPosDiffs;
      currConstraint.resolvePositionConstraint(constraintPositions, constraintVelocities, generatedPosDiffs, tolerance);
      
      //correcting positions according to impulses
      for (int i=0;i<currConstraint.globalIndices.size();i++){
        int currIndex=currConstraint.globalIndices(i);
        globalPositions(currIndex)+=generatedPosDiffs(i);
      }
    }
    
    global2Mesh();
    
    //Extensions: you can add user constraints in here and resolve them sequentially.
    
  }
  
  //adding an object.
  void addMesh(const MatrixXd& V, const MatrixXi& boundF, const MatrixXi& T, const double youngModulus, const double PoissonRatio,  const double density, const bool isFixed, const RowVector3d& userCOM, const RowVector4d userOrientation){
    
    VectorXd Vxyz(3*V.rows());
    for (int i=0;i<V.rows();i++)
      Vxyz.segment(3*i,3)=V.row(i).transpose();
    
    //cout<<"Vxyz: "<<Vxyz<<endl;
    Mesh m(Vxyz,boundF, T, globalPositions.size(), youngModulus, PoissonRatio, density, isFixed, userCOM, userOrientation);
    meshes.push_back(m);
    int oldTsize=globalT.rows();
    globalT.conservativeResize(globalT.rows()+T.rows(),4);
    globalT.block(oldTsize,0,T.rows(),4)=T.array()+globalPositions.size()/3;  //to offset T to global index
    globalPositions.conservativeResize(globalPositions.size()+Vxyz.size());
    globalVelocities.conservativeResize(globalPositions.size());
    int oldIMsize=globalInvMasses.size();
    globalInvMasses.conservativeResize(globalPositions.size());
    for (int i=0;i<m.invMasses.size();i++)
      globalInvMasses.segment(oldIMsize+3*i,3)=Vector3d::Constant(m.invMasses(i));
    
    mesh2global();
  }
  
  //loading a scene from the scene .txt files
  //you do not need to update this function
  bool loadScene(const std::string dataFolder, const std::string sceneFileName, const std::string constraintFileName){
    
    ifstream sceneFileHandle;
    ifstream constraintFileHandle;
    sceneFileHandle.open(dataFolder+std::string("/")+sceneFileName);
    if (!sceneFileHandle.is_open())
      return false;
    
    constraintFileHandle.open(dataFolder+std::string("/")+constraintFileName);
    if (!constraintFileHandle.is_open())
      return false;
    int numofObjects, numofConstraints;
    
    currTime=0;
    sceneFileHandle>>numofObjects;
    for (int i=0;i<numofObjects;i++){
      MatrixXi objT, objF;
      MatrixXd objV;
      std::string MESHFileName;
      bool isFixed;
      double youngModulus, poissonRatio, density;
      RowVector3d userCOM;
      RowVector4d userOrientation;
      sceneFileHandle>>MESHFileName>>density>>youngModulus>>poissonRatio>>isFixed>>userCOM(0)>>userCOM(1)>>userCOM(2)>>userOrientation(0)>>userOrientation(1)>>userOrientation(2)>>userOrientation(3);
      userOrientation.normalize();
      //if the mesh is an OFF file, tetrahedralize it
      if (MESHFileName.find(".off") != std::string::npos){
        MatrixXd VOFF;
        MatrixXi FOFF;
        igl::readOFF(dataFolder+std::string("/")+MESHFileName,VOFF,FOFF);
        RowVectorXd mins=VOFF.colwise().minCoeff();
        RowVectorXd maxs=VOFF.colwise().maxCoeff();
        for (int k=0;k<VOFF.rows();k++)
          VOFF.row(k)<<25.0*(VOFF.row(k)-mins).array()/(maxs-mins).array();
       
        if (!isFixed)
          igl::copyleft::tetgen::tetrahedralize(VOFF,FOFF,"pq1.1", objV,objT,objF);
        else
          igl::copyleft::tetgen::tetrahedralize(VOFF,FOFF,"pq1.414Y", objV,objT,objF);
      } else {
        igl::readMESH(dataFolder+std::string("/")+MESHFileName,objV,objT, objF);
      }
      
      //fixing weird orientation problem
      MatrixXi tempF(objF.rows(),3);
      tempF<<objF.col(2), objF.col(1), objF.col(0);
      objF=tempF;
      
      //cout<<"objF: "<<objF<<endl;
      //cout<<"viewerF: "<<viewerF<<endl;
      addMesh(objV,objF, objT, youngModulus, poissonRatio,  density, isFixed, userCOM, userOrientation);
    }
    
    //reading intra-mesh attachment constraints
    constraintFileHandle>>numofConstraints;
    for (int i=0;i<numofConstraints;i++){
      int attachM1, attachM2, attachV1, attachV2;
      constraintFileHandle>>attachM1>>attachV1>>attachM2>>attachV2;
      
      VectorXi coordIndices(6);
      coordIndices<<meshes[attachM1].globalOffset+3*attachV1,
      meshes[attachM1].globalOffset+3*attachV1+1,
      meshes[attachM1].globalOffset+3*attachV1+2,
      meshes[attachM2].globalOffset+3*attachV2,
      meshes[attachM2].globalOffset+3*attachV2+1,
      meshes[attachM2].globalOffset+3*attachV2+2;
   
      VectorXd constraintInvMasses(6);
      constraintInvMasses<<meshes[attachM1].invMasses(attachV1),
      meshes[attachM1].invMasses(attachV1),
      meshes[attachM1].invMasses(attachV1),
      meshes[attachM2].invMasses(attachV2),
      meshes[attachM2].invMasses(attachV2),
      meshes[attachM2].invMasses(attachV2);
      double refValue=(meshes[attachM1].currPositions.segment(3*attachV1,3)-meshes[attachM2].currPositions.segment(3*attachV2,3)).norm();
      userConstraints.push_back(Constraint(DISTANCE, EQUALITY, coordIndices, constraintInvMasses, VectorXd::Zero(0), refValue, 0.0));
    }
    return true;
  }
  
  
  Scene(){}
  ~Scene(){}
};



/*****************************Auxiliary functions for collision detection. Do not need updating********************************/

/** Support function for libccd*/
void support(const void *_obj, const ccd_vec3_t *_d, ccd_vec3_t *_p)
{
  // assume that obj_t is user-defined structure that holds info about
  // object (in this case box: x, y, z, pos, quat - dimensions of box,
  // position and rotation)
  //std::cout<<"calling support"<<std::endl;
  MatrixXd *obj = (MatrixXd *)_obj;
  RowVector3d p;
  RowVector3d d;
  for (int i=0;i<3;i++)
    d(i)=_d->v[i]; //p(i)=_p->v[i];
  
  
  d.normalize();
  //std::cout<<"d: "<<d<<std::endl;
  
  RowVector3d objCOM=obj->colwise().mean();
  int maxVertex=-1;
  int maxDotProd=-32767.0;
  for (int i=0;i<obj->rows();i++){
    double currDotProd=d.dot(obj->row(i)-objCOM);
    if (maxDotProd < currDotProd){
      maxDotProd=currDotProd;
      //std::cout<<"maxDotProd: "<<maxDotProd<<std::endl;
      maxVertex=i;
    }
    
  }
  //std::cout<<"maxVertex: "<<maxVertex<<std::endl;
  
  for (int i=0;i<3;i++)
    _p->v[i]=(*obj)(maxVertex,i);
  
  //std::cout<<"end support"<<std::endl;
}

void stub_dir(const void *obj1, const void *obj2, ccd_vec3_t *dir)
{
  dir->v[0]=1.0;
  dir->v[1]=0.0;
  dir->v[2]=0.0;
}

void center(const void *_obj,ccd_vec3_t *center)
{
  MatrixXd *obj = (MatrixXd *)_obj;
  RowVector3d objCOM=obj->colwise().mean();
  for (int i=0;i<3;i++)
    center->v[i]=objCOM(i);
}




#endif
