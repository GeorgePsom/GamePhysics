#ifndef SCENE_HEADER_FILE
#define SCENE_HEADER_FILE

#include <vector>
#include <fstream>
#include <igl/bounding_box.h>
#include <igl/readMESH.h>
#include <igl/bounding_box.h>
#include <igl/readOFF.h>
#include <igl/per_vertex_normals.h>
#include <igl/edge_topology.h>
#include <igl/diag.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include "ccd.h"
#include "volInt.h"
#include "auxfunctions.h"
#include "constraints.h"

using namespace Eigen;
using namespace std;


void support(const void *_obj, const ccd_vec3_t *_d, ccd_vec3_t *_p);
void stub_dir(const void *obj1, const void *obj2, ccd_vec3_t *dir);
void center(const void *_obj, ccd_vec3_t *dir);



//Impulse is defined as a pair <position, direction>
typedef std::pair<RowVector3d,RowVector3d> Impulse;


//the class the contains each individual rigid objects and their functionality
class Mesh{
public:
  MatrixXd origV;   //original vertex positions, where COM=(0.0,0.0,0.0) - never change this!
  MatrixXd currV;   //current vertex position
  MatrixXi F;   //faces of the tet mesh
  MatrixXi T;   //Tets in the tet mesh

  VectorXi boundTets;  //indices (from T) of just the boundary tets, for collision

  //position of object in space. We must always have that currV = QRot(origV, orientation)+ COM
  RowVector4d orientation; //current orientation
  RowVector3d COM;  //current center of mass
  RowVector3d COMprevious;
  Matrix3d invIT;  //Original *inverse* inertia tensor around the COM, defined in the rest state to the object (so to the canonical world system)

  VectorXd tetVolumes;    //|T|x1 tetrahedra volumes
  VectorXd invMasses;     //|T|x1 tetrahedra *inverse* masses

  //kinematics
  bool isFixed;  //is the object immobile
  double totalMass;  //sum(1/invMass)
  double totalVolume;
  RowVector3d comVelocity;  //the linear velocity of the center of mass
  RowVector3d angVelocity;  //the angular velocity of the object.

  //checking collision between bounding boxes, and consequently the boundary tets if succeeds.
  //you do not need to update these functions (isBoxCollide and isCollide) unless you are doing a different collision

  bool isBoxCollide(const Mesh& m){
    RowVector3d VMin1=currV.colwise().minCoeff();
    RowVector3d VMax1=currV.colwise().maxCoeff();
    RowVector3d VMin2=m.currV.colwise().minCoeff();
    RowVector3d VMax2=m.currV.colwise().maxCoeff();

    //checking all axes for non-intersection of the dimensional interval
    for (int i=0;i<3;i++)
      if ((VMax1(i)<VMin2(i))||(VMax2(i)<VMin1(i)))
        return false;

    return true;  //all dimensional intervals are overlapping = intersection

  }

  bool isCollide(const Mesh& m, double& depth, RowVector3d& intNormal, RowVector3d& intPosition){


    if ((isFixed && m.isFixed))  //collision does nothing
      return false;

    //collision between bounding boxes
    if (!isBoxCollide(m))
      return false;

    //otherwise, full test
    ccd_t ccd;
    CCD_INIT(&ccd);
    ccd.support1       = support; // support function for first object
    ccd.support2       = support; // support function for second object
    ccd.center1         =center;
    ccd.center2         =center;

    ccd.first_dir       = stub_dir;
    ccd.max_iterations = 100;     // maximal number of iterations


    void* obj1=(void*)this;
    void* obj2=(void*)&m;

    ccd_real_t _depth;
    ccd_vec3_t dir, pos;

    int nonintersect = ccdMPRPenetration(obj1, obj2, &ccd, &_depth, &dir, &pos);

    if (nonintersect)
      return false;

    for (int k=0;k<3;k++){
      intNormal(k)=dir.v[k];
      intPosition(k)=pos.v[k];
    }

    depth =_depth;
    intPosition-=depth*intNormal/2.0;

    //Vector3d p1=intPosition+depth*intNormal;
    //Vector3d p2=intPosition;
    //std::cout<<"intPosition: "<<intPosition<<std::endl;

    //std::cout<<"depth: "<<depth<<std::endl;
    //std::cout<<"After ccdGJKIntersect"<<std::endl;

    //return !nonintersect;

    return true;

  }


  //return the current inverted inertia tensor around the current COM. Update it by applying the orientation
  Matrix3d getCurrInvInertiaTensor(){
   /********
    TODO: complete from Practical 1

    *******/

      Matrix3d R = Q2RotMatrix(orientation);



      /***************
       TODO
       ***************/
      Matrix3d newI = R.transpose() * invIT * R;

      return newI;
  }


  //Update the current position and orientation by integrating the linear and angular velocities, and update currV accordingly
  //You need to modify this according to its purpose
  void updatePosition(double timeStep){
    //just forward Euler now
    if (isFixed)
      return;  //a fixed object is immobile

    /********
     TODO: complete from Practical 1
     *******/
    COMprevious = COM;
    COM += comVelocity * timeStep;

    RowVector3d normAng = angVelocity;

    //std::cout << "Angular: " << angVelocity << std::endl;
    //std::cout << "Norm angular: " << normAng << std::endl;

    RowVector4d q(0, normAng.x(), normAng.y(), normAng.z());

    orientation += 0.5 * timeStep * QMult(q, orientation);
    orientation.normalize();

    for (int i = 0; i < currV.rows(); i++)
    {
        currV.row(i) << QRot(origV.row(i), orientation) + COM;
    }

  }


  RowVector3d initStaticProperties(const double density)
  {
    //TODO: compute tet volumes and allocate to vertices
    tetVolumes.conservativeResize(T.rows());

    RowVector3d naturalCOM; naturalCOM.setZero();
    Matrix3d IT; IT.setZero();
    for (int i=0;i<T.rows();i++){
      Vector3d e01=origV.row(T(i,1))-origV.row(T(i,0));
      Vector3d e02=origV.row(T(i,2))-origV.row(T(i,0));
      Vector3d e03=origV.row(T(i,3))-origV.row(T(i,0));
      Vector3d tetCentroid=(origV.row(T(i,0))+origV.row(T(i,1))+origV.row(T(i,2))+origV.row(T(i,3)))/4.0;
      tetVolumes(i)=std::abs(e01.dot(e02.cross(e03)))/6.0;

      naturalCOM+=tetVolumes(i)*tetCentroid;

    }

    totalVolume=tetVolumes.sum();
    totalMass=density*totalVolume;
    naturalCOM.array()/=totalVolume;

    //computing inertia tensor
    for (int i=0;i<T.rows();i++){
      RowVector4d xvec; xvec<<origV(T(i,0),0)-naturalCOM(0),origV(T(i,1),0)-naturalCOM(0),origV(T(i,2),0)-naturalCOM(0),origV(T(i,3),0)-naturalCOM(0);
      RowVector4d yvec; yvec<<origV(T(i,0),1)-naturalCOM(1),origV(T(i,1),1)-naturalCOM(1),origV(T(i,2),1)-naturalCOM(1),origV(T(i,3),1)-naturalCOM(1);
      RowVector4d zvec; zvec<<origV(T(i,0),2)-naturalCOM(2),origV(T(i,1),2)-naturalCOM(2),origV(T(i,2),2)-naturalCOM(2),origV(T(i,3),2)-naturalCOM(2);

      double I00, I11, I22, I12, I21, I01, I10, I02, I20;
      Matrix4d sumMat=Matrix4d::Constant(1.0)+Matrix4d::Identity();
      I00 = density*6*tetVolumes(i)*(yvec*sumMat*yvec.transpose()+zvec*sumMat*zvec.transpose()).sum()/120.0;
      I11 = density*6*tetVolumes(i)*(xvec*sumMat*xvec.transpose()+zvec*sumMat*zvec.transpose()).sum()/120.0;
      I22 = density*6*tetVolumes(i)*(xvec*sumMat*xvec.transpose()+yvec*sumMat*yvec.transpose()).sum()/120.0;
      I12 = I21 = -density*6*tetVolumes(i)*(yvec*sumMat*zvec.transpose()).sum()/120.0;
      I10 = I01 = -density*6*tetVolumes(i)*(xvec*sumMat*zvec.transpose()).sum()/120.0;
      I20 = I02 = -density*6*tetVolumes(i)*(xvec*sumMat*yvec.transpose()).sum()/120.0;

      Matrix3d currIT; currIT<<I00, I01, I02,
      I10, I11, I12,
      I20, I21, I22;

      IT+=currIT;

    }
    invIT=IT.inverse();
    if (isFixed)
      invIT.setZero();  //infinite resistance to rotation

    return naturalCOM;

  }


  //Updating the linear and angular velocities of the object
  //You need to modify this to integrate from acceleration in the field (basically gravity)
  void updateVelocity(double timeStep){

    if (isFixed)
      return;

    /********
     TODO: complete from Practical 1
     *******/

    Vector3d gravity; gravity << 0, -9.8, 0.0;
    comVelocity += (gravity * timeStep);
    //std::cout << totalMass << std::endl;
    //comVelocity -= comVelocity * dragCoeff * timeStep;
    //angVelocity -= (angVelocity * timeStep * dragCoeff);
  }


  //the full integration for the time step (velocity + position)
  //You need to modify this if you are changing the integration
  void integrate(double timeStep){
    updateVelocity(timeStep);
    updatePosition(timeStep);
  }


  Mesh(const MatrixXd& _V, const MatrixXi& _F, const MatrixXi& _T, const double density, const bool _isFixed, const RowVector3d& _COM, const RowVector4d& _orientation){
    origV=_V;
    F=_F;
    T=_T;
    isFixed=_isFixed;
    COM=_COM;
    orientation=_orientation;
    comVelocity.setZero();
    angVelocity.setZero();

    COMprevious = COM;

    RowVector3d naturalCOM;  //by the geometry of the object

    //initializes the original geometric properties (COM + IT) of the object
    naturalCOM = initStaticProperties(density);

    origV.rowwise()-=naturalCOM;  //removing the natural COM of the OFF file (natural COM is never used again)

    currV.resize(origV.rows(), origV.cols());
    for (int i=0;i<currV.rows();i++){
        if (!isFixed)
        {
            origV.row(i) *= 0.25;
        }

        if (density == 111) {
            origV.row(i) *= 4.0;
        }

        currV.row(i)<<QRot(origV.row(i), orientation)+COM;
    }


    VectorXi boundVMask(origV.rows());
    boundVMask.setZero();
    for (int i=0;i<F.rows();i++)
      for (int j=0;j<3;j++)
        boundVMask(F(i,j))=1;

    //cout<<"boundVMask.sum(): "<<boundVMask.sum()<<endl;

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
  }

  ~Mesh(){}
};

//This class contains the entire scene operations, and the engine time loop.
class Scene{
public:
  double currTime = 0.0;
  int numFullV, numFullT;
  std::vector<Mesh> meshes;
  const double EPSILON = 1.0e-4f;

  double poly6(RowVector3d r, double h)
  {
      double rBar = r.norm();

      if (rBar < EPSILON || rBar > h) {
          return 0.0f;
      }

      // (315 / (64 * PI * h^9)) * (h^2 - |r|^2)^3
      double h9 = (h * h * h * h * h * h * h * h * h);
      if (h9 < EPSILON) {
          return 0.0f;
      }
      double A = 1.566681471061f / h9;
      double B = (h * h) - (rBar * rBar);

      return A * (B * B * B);
  }

  RowVector3d spiky(RowVector3d r, double h)
  {
      double rBar = r.norm();

      if (rBar < EPSILON || rBar > h) {
          return RowVector3d::Zero();
      }

      // (45 / (PI * h^6)) * (h - |r|)^2 * (r / |r|)
      double h6 = (h * h * h * h * h * h);
      if (h6 < EPSILON) {
          return RowVector3d::Zero();
      }
      double A = 14.323944878271f / h6;
      double B = (h - rBar);
      RowVector3d out = A * (B * B) * (r / (rBar + EPSILON));
      //out[3] = 0.0f;
      return out;
  }

  //Practical 2
  vector<Constraint> constraints;   //The (user) constraints of the scene

  //adding an objects. You do not need to update this generally
  void addMesh(const MatrixXd& V, const MatrixXi& F, const MatrixXi& T, const double density, const bool isFixed, const RowVector3d& COM, const RowVector4d& orientation){

    Mesh m(V,F, T, density, isFixed, COM, orientation);
    meshes.push_back(m);
  }



  /*********************************************************************
   This function handles collision constraints between objects m1 and m2 when found
   Input: meshes m1, m2
   depth: the depth of penetration
   contactNormal: the normal of the conact measured m1->m2
   penPosition: a point on m2 such that if m2 <= m2 + depth*contactNormal, then penPosition+depth*contactNormal is the common contact point
   CRCoeff: the coefficient of restitution

   You should create a "Constraint" class, and use its resolveVelocityConstraint() and resolvePositionConstraint() *alone* to resolve the constraint.
   You are not allowed to use practical 1 collision handling
   *********************************************************************/
  void handleCollision(Mesh& m1, Mesh& m2,const double& depth, const RowVector3d& contactNormal,const RowVector3d& penPosition, const double CRCoeff, const double tolerance){


    //std::cout<<"Collision: " << std::endl;
    //std::cout<<"penPosition: "<<penPosition<<std::endl;

    double invMass1 = (m1.isFixed ? 0.0 : 1.0/m1.totalMass);  //fixed meshes have infinite mass
    double invMass2 = (m2.isFixed ? 0.0 : 1.0/m2.totalMass);

    RowVector3d com1 = m1.COM;
    RowVector3d com2 = m2.COM;

    RowVector3d contactPosition;
    if (m1.isFixed) {

        com2 = m2.COM + depth * contactNormal.normalized();
        contactPosition = penPosition + depth * contactNormal.normalized();
    }
    else if (m2.isFixed) {

        com1 = m1.COM - depth * contactNormal.normalized();
        contactPosition = penPosition + depth * contactNormal.normalized();

    }
    else { //inverse mass weighting

        com2 = m2.COM + ((depth * m2.totalMass) / (m1.totalMass + m2.totalMass)) * contactNormal;
        com1 = m1.COM - ((depth * m1.totalMass) / (m1.totalMass + m2.totalMass)) * contactNormal;

        contactPosition = penPosition + ((depth * m2.totalMass) / (m1.totalMass + m2.totalMass)) * contactNormal;
    }

    //m1.COM = com1;
    //m2.COM = com2;
    double vel1 = (m1.comVelocity.y() > 0 ? 1 : -1) * sqrt(abs(m1.comVelocity.y() * m1.comVelocity.y() - 2 * 9.8 * depth));
    double vel2 = (m2.comVelocity.y() > 0 ? 1 : -1) * sqrt(abs(m2.comVelocity.y() * m2.comVelocity.y() - 2 * 9.8 * depth));
     // collision arm
    RowVector3d r1 = contactPosition - com1;
    RowVector3d r2 = contactPosition - com2;
    m1.comVelocity << m1.comVelocity.x(), m1.isFixed ? m1.comVelocity.y() : vel1, m1.comVelocity.z();
    m2.comVelocity << m2.comVelocity.x(), m2.isFixed ? m2.comVelocity.y() : vel2, m2.comVelocity.z();

    /***************
     TODO: practical 2
     update m(1,2) comVelocity, angVelocity and COM variables by using a Constraint class of type COLLISION
     ***********************/

    Constraint collisionConstraint(COLLISION, INEQUALITY, 0, 0, 0, 0, invMass1, invMass2, -contactNormal, depth, CRCoeff);
    Matrix<double, 2, 3> currCOMPosition;
    currCOMPosition(0,0) = m1.COM.x();
    currCOMPosition(0, 1) = m1.COM.y();
    currCOMPosition(0, 2) = m1.COM.z();
    currCOMPosition(1, 0) = m2.COM.x();
    currCOMPosition(1, 1) = m2.COM.y();
    currCOMPosition(1, 2) = m2.COM.z();

    //std::cout << m2.COM.y() << std::endl;
    Matrix<double, 2, 3> currVertexPositions;
    currVertexPositions(0, 0) = r1.x();
    currVertexPositions(0, 1) = r1.y();
    currVertexPositions(0, 2) = r1.z();
    currVertexPositions(1, 0) = r2.x();
    currVertexPositions(1, 1) = r2.y();
    currVertexPositions(1, 2) = r2.z();


    Matrix<double, 2, 3> currCOMVelocities;
    currCOMVelocities(0, 0) = m1.comVelocity.x();
    currCOMVelocities(0, 1) = m1.comVelocity.y();
    currCOMVelocities(0, 2) = m1.comVelocity.z();
    currCOMVelocities(1, 0) = m2.comVelocity.x();
    currCOMVelocities(1, 1) = m2.comVelocity.y();
    currCOMVelocities(1, 2) = m2.comVelocity.z();


    Matrix<double, 2, 3> currAngularVelocities;
    currAngularVelocities(0, 0) = m1.angVelocity.x();
    currAngularVelocities(0, 1) = m1.angVelocity.y();
    currAngularVelocities(0, 2) = m1.angVelocity.z();
    currAngularVelocities(1, 0) = m2.angVelocity.x();
    currAngularVelocities(1, 1) = m2.angVelocity.y();
    currAngularVelocities(1, 2) = m2.angVelocity.z();

    MatrixXd correctedCOMVelocities =MatrixXd::Zero(2, 3);
    MatrixXd correctedAngularVelocities = MatrixXd::Zero(2, 3);
    MatrixXd correctedCOMPositions = MatrixXd::Zero(2, 3);
    collisionConstraint.resolveVelocityConstraint(currCOMPosition, currVertexPositions, currCOMVelocities,
        currAngularVelocities, m1.getCurrInvInertiaTensor(), m2.getCurrInvInertiaTensor(), correctedCOMVelocities,
        correctedAngularVelocities, 0.0001);


    if (!m1.isFixed)
    {
        m1.comVelocity = correctedCOMVelocities.row(0);
        m1.angVelocity = correctedAngularVelocities.row(0);
    }
    if (!m2.isFixed)
    {
        m2.comVelocity = correctedCOMVelocities.row(1);
        m2.angVelocity = correctedAngularVelocities.row(1);
    }

    Matrix<double, 2, 3> currVertexPositionsN;
    currVertexPositionsN(0, 0) = penPosition.x();
    currVertexPositionsN(0, 1) = penPosition.y();
    currVertexPositionsN(0, 2) = penPosition.z();
    currVertexPositionsN(1, 0) = contactPosition.x();
    currVertexPositionsN(1, 1) = contactPosition.y();
    currVertexPositionsN(1, 2) = contactPosition.z();

    collisionConstraint.resolvePositionConstraint(currCOMPosition, currVertexPositionsN, correctedCOMPositions, 0.01f);
    if (!m1.isFixed)
        m1.COM = correctedCOMPositions.row(0);
    if (!m2.isFixed)
        m2.COM = correctedCOMPositions.row(1);

  }

  /*********************************************************************
   This function handles a single time step by:
   1. Integrating velocities, positions, and orientations by the timeStep
   2. (Practical 2) Detecting collisions and encoding them as constraints
   3. (Practical 2) Iteratively resolved positional and velocity constraints

   You do not need to update this function in Practical 2
   *********************************************************************/
  void updateScene(const double timeStep, const double CRCoeff, const double tolerance, const int maxIterations, vector<vector<RowVector3d>>& dataArray){

    //integrating velocity, position and orientation from forces and previous states
    //omp_set_num_threads(omp_get_max_threads());
#pragma omp parallel for
    for (int i=0;i<meshes.size();i++)
      meshes[i].integrate(timeStep);
#pragma omp barrier

    //detecting and handling collisions when found
    //This is done exhaustively: checking every two objects in the scene.
    double depth;
    RowVector3d contactNormal, penPosition;
    //omp_set_num_threads(omp_get_max_threads());


    for (int i=0;i<meshes.size();i++)
      for (int j=i+1;j<meshes.size();j++)
        if (meshes[i].isCollide(meshes[j],depth, contactNormal, penPosition))
          handleCollision(meshes[i], meshes[j],depth, contactNormal, penPosition,CRCoeff, tolerance);



    const int N = 3;
    const double REST_DENSITY = 37.76;
    const double INV_REST_DENSITY = 1.0f / REST_DENSITY;
    const double particlesRadius = 2.34375;
    const double smoothingRadius = 5.0;
    //const double particlesRadius = 9.5;
    //const double smoothingRadius = 20.0;
    const double relaxation = 0.0033;
    const double artificialPressureK = 0.1;
    const double artificalPressureN = 5;
    const double vorticity = 0.1;
    const double viscocity = 0.01;


    vector<double> Densities(meshes.size());
    vector<double> Lambdas(meshes.size());
    vector<RowVector3d> DeltaPs(meshes.size());
    vector<RowVector3d> Curls(meshes.size());
    vector<RowVector3d> Vorticities(meshes.size());
    vector<RowVector3d> Viscocities(meshes.size());
    for (int iter = 0; iter < N; iter++)
    {
        //estimate density
omp_set_num_threads(omp_get_max_threads());
#pragma omp parallel for
        for (int i = 0; i < meshes.size(); i++)
        {
            Densities[i] = 0;
            if (meshes[i].isFixed)
                continue;

            for (int j = 0; j < meshes.size(); j++)
            {
                if (meshes[j].isFixed || i == j)
                    continue;
                Densities[i] += poly6(meshes[i].COM - meshes[j].COM, smoothingRadius);
            }

        }
#pragma omp barrier


#pragma omp parallel for
        // Compute Lambda
        for (int i = 0; i < meshes.size(); i++)
        {
            if (meshes[i].isFixed)
                continue;
            double C_i = Densities[i] * INV_REST_DENSITY - 1.0;
            double gradientS = 0.0;

            RowVector3d gv_i;
            gv_i << 0.0, 0.0, 0.0;
            for (int j = 0; j < meshes.size(); j++)
            {
                if (meshes[j].isFixed || i == j)
                    continue;
                gv_i += spiky(meshes[i].COM - meshes[j].COM, smoothingRadius);
            }
            double gv_iLen = INV_REST_DENSITY * gv_i.norm();
            gradientS += gv_iLen * gv_iLen;

            double gv_sLengths = 0.0;

            for (int j = 0; j < meshes.size(); j++)
            {
                if (meshes[j].isFixed || i == j)\
                    continue;
                RowVector3d gradVector = -spiky(meshes[i].COM - meshes[j].COM, smoothingRadius) * INV_REST_DENSITY;
                double gradVectorL = gradVector.norm();
                gv_sLengths += gradVectorL;
            }
            gradientS += gv_sLengths;
            if (gradientS == 0.0)
                gradientS == EPSILON;
            Lambdas[i] = -C_i / (gradientS + relaxation);

        }
#pragma omp barrier

        // Compute Curls
#pragma omp parallel for
        for (int i = 0; i < meshes.size(); i++)
        {
            if (meshes[i].isFixed)
                continue;
            Curls[i] = RowVector3d::Zero();
            for (int j = 0; j < meshes.size(); j++)
            {
                if (i == j || meshes[j].isFixed)
                    continue;
                RowVector3d Vel12 = meshes[i].comVelocity - meshes[j].comVelocity;
                RowVector3d gradient = spiky(meshes[i].COM - meshes[j].COM, smoothingRadius);
                Curls[i] += Vel12.cross(gradient);
            }
        }
#pragma omp barrier

#pragma omp parallel for
        // Compute Dp
        for (int i = 0; i < meshes.size(); i++)
        {
            if (meshes[i].isFixed)
                continue;
            DeltaPs[i] = RowVector3d::Zero();
            Vorticities[i] = RowVector3d::Zero();
            Viscocities[i] = RowVector3d::Zero();
            double lambda_i = Lambdas[i];
            for (int j = 0; j < meshes.size(); j++)
            {
                if (i == j || meshes[j].isFixed)
                    continue;
                // Artificial Pressure Corrector
                double h = smoothingRadius;
                RowVector3d r = meshes[i].COM - meshes[j].COM;
                RowVector3d gradient = spiky(r, h);
                double n = poly6(r, h);

                double offset = (0.3f * h);
                RowVector3d deltaQ = meshes[i].COM + RowVector3d(offset, offset, offset);
                double d = poly6(deltaQ, h);
                double nd = abs(d) <= EPSILON ? 0.0f : n / d;

                double nd2 = nd * nd;
                double s_corr = -artificialPressureK * pow(nd, artificialPressureK);
                DeltaPs[i] += (lambda_i + Lambdas[j] + s_corr) * gradient;

                RowVector3d curl_ij = Curls[i] - Curls[j];
                double omegaBar = curl_ij.norm();
                Vorticities[i] += RowVector3d(omegaBar / r.x(), omegaBar / r.y(), omegaBar / r.z());

                RowVector3d Vel12 = meshes[i].comVelocity - meshes[j].comVelocity;
                double W_ij = poly6(r, smoothingRadius);
                Viscocities[i] += W_ij * Vel12;

            }

            double n = Vorticities[i].norm();
            RowVector3d N = RowVector3d::Zero();
            if (n > EPSILON)
                N = Vorticities[i].normalized();
            RowVector3d f_curl = vorticity * N.cross(Curls[i]);
            Vorticities[i] = timeStep * f_curl;
            Viscocities[i] *= viscocity;



        }
#pragma omp barrier

#pragma omp parallel for
        for (int i = 0; i < meshes.size(); i++)
        {
            if (meshes[i].isFixed)
                continue;
            meshes[i].COM += DeltaPs[i];
        }
#pragma omp barrier
    }
//
#pragma omp parallel for
    for (int i = 0; i < meshes.size(); i++)
    {
        if (meshes[i].isFixed)
            continue;
        meshes[i].comVelocity = 1.0 / 0.02 * (meshes[i].COM - meshes[i].COMprevious) + Viscocities[i] + Vorticities[i];

        for (int j = 0; j < meshes[i].currV.rows(); j++)
        {
            meshes[i].currV.row(j) += DeltaPs[i];
        }
    }
#pragma omp barrier

    ////Resolving user constraints iteratively until either:
    ////1. Positions or velocities are valid up to tolerance (a full streak of validity in the iteration)
    ////2. maxIterations has run out
    //
    //
    ////Resolving velocity
    //int currIteration=0;
    //int zeroStreak=0;  //how many consecutive constraints are already below tolerance without any change; the algorithm stops if all are.currVertexPositions
    //int currConstIndex=0;
    //while ((zeroStreak<constraints.size())&&(currIteration*constraints.size()<maxIterations)){
    //
    //  Constraint currConstraint=constraints[currConstIndex];
    //
    //  RowVector3d origConstPos1=meshes[currConstraint.m1].origV.row(currConstraint.v1);
    //  RowVector3d origConstPos2=meshes[currConstraint.m2].origV.row(currConstraint.v2);
    //
    //  RowVector3d currConstPos1 = QRot(origConstPos1, meshes[currConstraint.m1].orientation)+meshes[currConstraint.m1].COM;
    //  RowVector3d currConstPos2 = QRot(origConstPos2, meshes[currConstraint.m2].orientation)+meshes[currConstraint.m2].COM;
    //  //cout<<"(currConstPos1-currConstPos2).norm(): "<<(currConstPos1-currConstPos2).norm()<<endl;
    //  //cout<<"(meshes[currConstraint.m1].currV.row(currConstraint.v1)-meshes[currConstraint.m2].currV.row(currConstraint.v2)).norm(): "<<(meshes[currConstraint.m1].currV.row(currConstraint.v1)-meshes[currConstraint.m2].currV.row(currConstraint.v2)).norm()<<endl;
    //  MatrixXd currCOMPositions(2,3); currCOMPositions<<meshes[currConstraint.m1].COM, meshes[currConstraint.m2].COM;
    //  MatrixXd currConstPositions(2,3); currConstPositions<<currConstPos1, currConstPos2;
    //  MatrixXd currCOMVelocities(2,3); currCOMVelocities<<meshes[currConstraint.m1].comVelocity, meshes[currConstraint.m2].comVelocity;
    //  MatrixXd currAngVelocities(2,3); currAngVelocities<<meshes[currConstraint.m1].angVelocity, meshes[currConstraint.m2].angVelocity;
    //
    //  Matrix3d invInertiaTensor1=meshes[currConstraint.m1].getCurrInvInertiaTensor();
    //  Matrix3d invInertiaTensor2=meshes[currConstraint.m2].getCurrInvInertiaTensor();
    //  MatrixXd correctedCOMVelocities = MatrixXd::Zero(2, 3);
    //  MatrixXd correctedAngVelocities = MatrixXd::Zero(2, 3);
    //
    //
    //  bool velocityWasValid=currConstraint.resolveVelocityConstraint(currCOMPositions, currConstPositions, currCOMVelocities, currAngVelocities, invInertiaTensor1, invInertiaTensor2, correctedCOMVelocities,correctedAngVelocities, tolerance);
    //
    //  if (velocityWasValid){
    //    zeroStreak++;
    //  }else{
    //    //only update the COM and angular velocity, don't both updating all currV because it might change again during this loop!
    //    zeroStreak=0;
    //    meshes[currConstraint.m1].comVelocity =correctedCOMVelocities.row(0);
    //    meshes[currConstraint.m2].comVelocity =correctedCOMVelocities.row(1);
    //
    //    meshes[currConstraint.m1].angVelocity =correctedAngVelocities.row(0);
    //    meshes[currConstraint.m2].angVelocity =correctedAngVelocities.row(1);
    //
    //  }
    //
    //  currIteration++;
    //  currConstIndex=(currConstIndex+1)%(constraints.size());
    //}
    //
    //if (currIteration*constraints.size()>=maxIterations)
    //  cout<<"Velocity Constraint resolution reached maxIterations without resolving!"<<endl;
    //
    //
    ////Resolving position
    //currIteration=0;
    //zeroStreak=0;  //how many consecutive constraints are already below tolerance without any change; the algorithm stops if all are.
    //currConstIndex=0;
    //while ((zeroStreak<constraints.size())&&(currIteration*constraints.size()<maxIterations)){
    //
    //  Constraint currConstraint=constraints[currConstIndex];
    //
    //  RowVector3d origConstPos1=meshes[currConstraint.m1].origV.row(currConstraint.v1);
    //  RowVector3d origConstPos2=meshes[currConstraint.m2].origV.row(currConstraint.v2);
    //
    //  RowVector3d currConstPos1 = QRot(origConstPos1, meshes[currConstraint.m1].orientation)+meshes[currConstraint.m1].COM;
    //  RowVector3d currConstPos2 = QRot(origConstPos2, meshes[currConstraint.m2].orientation)+meshes[currConstraint.m2].COM;

    //  MatrixXd currCOMPositions(2,3); currCOMPositions<<meshes[currConstraint.m1].COM, meshes[currConstraint.m2].COM;
    //  MatrixXd currConstPositions(2,3); currConstPositions<<currConstPos1, currConstPos2;
    //
    //  MatrixXd correctedCOMPositions = MatrixXd::Zero(2, 3);
    //
    //  bool positionWasValid=currConstraint.resolvePositionConstraint(currCOMPositions, currConstPositions,correctedCOMPositions, tolerance);
    //
    //  if (positionWasValid){
    //    zeroStreak++;
    //  }else{
    //    //only update the COM and angular velocity, don't both updating all currV because it might change again during this loop!
    //    zeroStreak=0;

    //    meshes[currConstraint.m1].COM =correctedCOMPositions.row(0);
    //    meshes[currConstraint.m2].COM =correctedCOMPositions.row(1);
    //
    //  }
    //
    //  currIteration++;
    //  currConstIndex=(currConstIndex+1)%(constraints.size());
    //}
    //
    //if (currIteration*constraints.size()>=maxIterations)
    //  cout<<"Position Constraint resolution reached maxIterations without resolving!"<<endl;
    //
    //
    ////updating currV according to corrected COM
    //for (int i=0;i<meshes.size();i++)
    //  for (int j=0;j<meshes[i].currV.rows();j++)
    //    meshes[i].currV.row(j)<<QRot(meshes[i].origV.row(j), meshes[i].orientation)+meshes[i].COM;

    int index = currTime / timeStep;

    if (currTime == 0) index = 0;

    for (int i = 0; i < meshes.size(); i++) {
        dataArray[index][i] = meshes[i].COM;
    }

    currTime+=timeStep;
  }

  //loading a scene from the scene .txt files
  //you do not need to update this function
  bool loadScene(const std::string dataFolder, const std::string sceneFileName, const std::string constraintFileName, bool stretch = false){

    ifstream sceneFileHandle, constraintFileHandle;
    sceneFileHandle.open(dataFolder+std::string("/")+sceneFileName);
    if (!sceneFileHandle.is_open())
      return false;
    int numofObjects;

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
      if (MESHFileName.find(".off") != std::string::npos) {
          MatrixXd VOFF;
          MatrixXi FOFF;
          igl::readOFF(dataFolder + std::string("/") + MESHFileName, VOFF, FOFF);
          RowVectorXd mins = VOFF.colwise().minCoeff();
          RowVectorXd maxs = VOFF.colwise().maxCoeff();
          for (int k = 0; k < VOFF.rows(); k++)
              VOFF.row(k) << 25.0 * (VOFF.row(k) - mins).array() / (maxs - mins).array();

          if (!isFixed)
              igl::copyleft::tetgen::tetrahedralize(VOFF, FOFF, "pq1.1", objV, objT, objF);
          else
              igl::copyleft::tetgen::tetrahedralize(VOFF, FOFF, "pq1.414Y", objV, objT, objF);
      }
      else {
          igl::readMESH(dataFolder + std::string("/") + MESHFileName, objV, objT, objF);
      }
      //igl::readMESH(dataFolder+std::string("/")+MESHFileName,objV,objT, objF);

      //fixing weird orientation problem
      MatrixXi tempF(objF.rows(),3);
      tempF<<objF.col(2), objF.col(1), objF.col(0);
      objF=tempF;


      //addMesh(objV,objF, objT,density, isFixed, userCOM, userOrientation);

      if (i == 1) {
          for (int k = 0; k < 10; k++)
          {
              for (int j = 0; j < 5 - k / 2; j++)
              {
                  for (int i = 0; i < 5 - k / 2; i++)
                  {
                      RowVector3d com;
                      com << userCOM.x() + i * 8.0, userCOM.y() + 8.0 * k, userCOM.z() + 8.0 * j;
                      //com << userCOM.x() + i * 25.0, userCOM.y() + 25.0 * k, userCOM.z() + 25.0 * j;
                      addMesh(objV, objF, objT, density, isFixed, com, userOrientation);
                  }
              }
          }
      }
      else if (i==0) {
          RowVector3d com;
          com << userCOM.x(), userCOM.y()- 160.0 + 160 /20 + 30, userCOM.z();
          addMesh(objV, objF, objT, density, isFixed, com, userOrientation);
      }

      cout << "COM: " << userCOM <<endl;
      cout << "orientation: " << userOrientation <<endl;
    }

    //Practical 2 change
    //reading intra-mesh attachment constraints
    int numofConstraints;
    constraintFileHandle.open(dataFolder+std::string("/")+constraintFileName);
    if (!constraintFileHandle.is_open())
      return false;
    constraintFileHandle>>numofConstraints;
    for (int i=0;i<numofConstraints;i++){
      int attachM1, attachM2, attachV1, attachV2;
      constraintFileHandle>>attachM1>>attachV1>>attachM2>>attachV2;

      double initDist=(meshes[attachM1].currV.row(attachV1)-meshes[attachM2].currV.row(attachV2)).norm();
      cout<<"initDist: "<<initDist<<endl;
      double invMass1 = (meshes[attachM1].isFixed ? 0.0 : 1.0/meshes[attachM1].totalMass);  //fixed meshes have infinite mass
      double invMass2 = (meshes[attachM2].isFixed ? 0.0 : 1.0/meshes[attachM2].totalMass);
      if (stretch)
      {
          constraints.push_back(Constraint(STRETCH, EQUALITY, attachM1, attachV1, attachM2, attachV2, invMass1, invMass2, RowVector3d::Zero(), initDist, 0.0));
        //  constraints.push_back(Constraint(COMPRESS, INEQUALITY, attachM1, attachV1, attachM2, attachV2, invMass1, invMass2, RowVector3d::Zero(), initDist - 0.5 * initDist, 0.0));

      }
      else
          constraints.push_back(Constraint(DISTANCE, EQUALITY, attachM1, attachV1, attachM2, attachV2, invMass1, invMass2, RowVector3d::Zero(), initDist, 0.0));


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
  Mesh *obj = (Mesh *)_obj;
  RowVector3d p;
  RowVector3d d;
  for (int i=0;i<3;i++)
    d(i)=_d->v[i]; //p(i)=_p->v[i];

  d.normalize();
  //std::cout<<"d: "<<d<<std::endl;

  int maxVertex=-1;
  int maxDotProd=-32767.0;
  for (int i=0;i<obj->currV.rows();i++){
    double currDotProd=d.dot(obj->currV.row(i)-obj->COM);
    if (maxDotProd < currDotProd){
      maxDotProd=currDotProd;
      //std::cout<<"maxDotProd: "<<maxDotProd<<std::endl;
      maxVertex=i;
    }

  }
  //std::cout<<"maxVertex: "<<maxVertex<<std::endl;

  for (int i=0;i<3;i++)
    _p->v[i]=obj->currV(maxVertex,i);

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
  Mesh *obj = (Mesh *)_obj;
  for (int i=0;i<3;i++)
    center->v[i]=obj->COM(i);
}








#endif
