#ifndef CONSTRAINTS_HEADER_FILE
#define CONSTRAINTS_HEADER_FILE

using namespace Eigen;
using namespace std;

typedef enum ConstraintType{DISTANCE, COLLISION} ConstraintType;   //You can expand it for more constraints
typedef enum ConstraintEqualityType{EQUALITY, INEQUALITY} ConstraintEqualityType;

//there is such constraints per two variables that are equal. That is, for every attached vertex there are three such constraints for (x,y,z);
class Constraint{
public:
  
  int m1, m2;                     //Two participating meshes (can be the same)  - auxiliary data for users (constraint class shouldn't use that)
  int v1, v2;                     //Two vertices from the respective meshes - auxiliary data for users (constraint class shouldn't use that)
  double invMass1, invMass2;       //inverse masses of two bodies
  double refValue;                //Reference values to use in the constraint, when needed (like distance)
  RowVector3d refVector;             //Reference vector when needed (like vector)
  double CRCoeff;                 //extra velocity bias
  ConstraintType constraintType;  //The type of the constraint, and will affect the value and the gradient. This SHOULD NOT change after initialization!
  ConstraintEqualityType constraintEqualityType;  //whether the constraint is an equality or an inequality
  
  Constraint(const ConstraintType _constraintType, const ConstraintEqualityType _constraintEqualityType, const int& _m1, const int& _v1, const int& _m2, const int& _v2, const double& _invMass1, const double& _invMass2, const RowVector3d& _refVector, const double& _refValue, const double& _CRCoeff):constraintType(_constraintType), constraintEqualityType(_constraintEqualityType), m1(_m1), v1(_v1), m2(_m2), v2(_v2), invMass1(_invMass1), invMass2(_invMass2),  refValue(_refValue), CRCoeff(_CRCoeff){
    refVector=_refVector;
  }
  
  ~Constraint(){}
  
  
  
  //computes the impulse needed for all particles to resolve the velocity constraint, and corrects the velocities accordingly.
  //The velocities are a vector (vCOM1, w1, vCOM2, w2) in both input and output.
  //returns true if constraint was already valid with "currVelocities", and false otherwise (false means there was a correction done)
  //currCOMPositions is a 2x3 matrix, where each row is per one of the sides of the constraints; the rest of the relevant variables are similar, and so should the outputs be resized.
  bool resolveVelocityConstraint(const MatrixXd& currCOMPositions, const MatrixXd& currVertexPositions, 
      const MatrixXd& currCOMVelocities, const MatrixXd& currAngularVelocities, const Matrix3d& invInertiaTensor1,
      const Matrix3d& invInertiaTensor2, MatrixXd& correctedCOMVelocities, MatrixXd& correctedAngularVelocities,
      double tolerance){
      //std::cout << "Resolve velocity" << std::endl;
    MatrixXd invMassMatrix=MatrixXd::Zero(12,12);
    RowVectorXd constGradient(12);
    RowVector3d v1 = currVertexPositions.row(0);
    RowVector3d v2 = currVertexPositions.row(1);

    Vector3d r1N = v1.cross(refVector);
    Vector3d r2N = -v2.cross(refVector);
    if (constraintType == COLLISION)
    {
        constGradient << refVector.x(), refVector.y(), refVector.z(),
            r1N.x(), r1N.y(), r1N.z(),
            -refVector.x(), -refVector.y(), -refVector.z(),
            r2N.x(), r2N.y(), r2N.z();
    }
    else
    {
        RowVector3d d12 = currVertexPositions.row(0) - currVertexPositions.row(1);
        float Cp = d12.norm() - refValue;
        RowVector3d r1 = currVertexPositions.row(0) - currCOMPositions.row(0);
        RowVector3d r2 = currVertexPositions.row(0) - currCOMPositions.row(0);
        d12.normalize();

        Vector3d r1N = v1.cross(d12);
        Vector3d r2N = -v2.cross(d12);
        constGradient << d12.x(), d12.y(), d12.z(),
            r1N.x(), r1N.y(), r1N.z(),
            -d12.x(), -d12.y(), -d12.z(),
            r2N.x(), r2N.y(), r2N.z();


      
    }
       

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if(i == j)
                invMassMatrix(i, j) = invMass1;
        }
    }

    for (int i = 6; i < 9; i++)
    {
        for (int j = 6; j < 9; j++)
        {
            if(i == j)
                invMassMatrix(i, j) = invMass2;
        }
    }


    //invMassMatrix(0, 0) = invMass1;
    //invMassMatrix(1, 1) = invMass1;
    //invMassMatrix(2, 2) = invMass1;
    //invMassMatrix(6, 6) = invMass2;
    //invMassMatrix(7, 7) = invMass2;
    //invMassMatrix(8, 8) = invMass2;

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {

            invMassMatrix(3 + i, 3 + j) = invInertiaTensor1(i, j);
        }
    }

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            invMassMatrix(9 + i, 9 + j) = invInertiaTensor2(i, j);
        }
    }

    RowVectorXd velocities(12);
    velocities << currCOMVelocities.row(0), currAngularVelocities.row(0), currCOMVelocities.row(1), currAngularVelocities.row(1);
    
    double Ju = constGradient * velocities.transpose();
    if ((abs(Ju) <= tolerance && constraintEqualityType == EQUALITY)
        || ((Ju) >= 0 && constraintEqualityType == INEQUALITY))
    {
        correctedCOMVelocities = currCOMVelocities;
        correctedAngularVelocities = currAngularVelocities;
        return true;
    }

    double lamda = -(1 + CRCoeff) * Ju / (constGradient * invMassMatrix * constGradient.transpose());

    velocities += invMassMatrix * constGradient.transpose() * lamda;

    correctedCOMVelocities.row(0) << velocities(0), velocities(1), velocities(2);
    correctedAngularVelocities.row(0) << velocities(3), velocities(4), velocities(5);
    correctedCOMVelocities.row(1) << velocities(6), velocities(7), velocities(8);
    correctedAngularVelocities.row(1) << velocities(9), velocities(10), velocities(11);


    return false;
    
    
    /**************
     TODO: write velocity correction procedure:
     1. If the velocity Constraint is satisfied up to tolerate ("abs(Jv)<=tolerance"), set corrected values to original ones and return true
     
     2. Otherwise, correct linear and angular velocities as learnt in class.
     
     Note to differentiate between different constraint types; for inequality constraints you don't do anything unless it's unsatisfied.
     ***************/
    
     //Stub code: remove upon implementation
    
     //end of stub code
  }
  
  //projects the position unto the constraint
  //returns true if constraint was already valid with "currPositions"
  bool resolvePositionConstraint(const MatrixXd& currCOMPositions, const MatrixXd& currConstPositions, MatrixXd& correctedCOMPositions, double tolerance){
    
      //std::cout << "Resolve position" << std::endl;

    MatrixXd invMassMatrix=MatrixXd::Zero(6,6);
    invMassMatrix(0, 0) = invMass1;
    invMassMatrix(1, 1) = invMass1;
    invMassMatrix(2, 2) = invMass1;
    invMassMatrix(3, 3) = invMass2;
    invMassMatrix(4, 4) = invMass2;
    invMassMatrix(5, 5) = invMass2;

    
    RowVectorXd constGradient(6);
    double Cp;
    if (constraintType == COLLISION)
    {
        constGradient << refVector, -refVector;
        Cp = (currCOMPositions.row(0) - currCOMPositions.row(1)).dot(refVector);
    }  
    else
    {
        RowVector3d d12 = currConstPositions.row(0) - currConstPositions.row(1);
        Cp = d12.norm() - refValue;
        d12.normalize();
        constGradient << d12, -d12;
    }
        
   
    if ((constraintEqualityType == INEQUALITY && Cp >= 0) ||
        constraintEqualityType == EQUALITY && abs(Cp) <= tolerance)
    {
        correctedCOMPositions = currCOMPositions;
        return true;
    }


    double lamda = -Cp / (constGradient * invMassMatrix * constGradient.transpose());
    RowVectorXd pos(6);
    pos << currCOMPositions.row(0), currCOMPositions.row(1);
    pos += (lamda * invMassMatrix * constGradient.transpose()).transpose();
    correctedCOMPositions.row(0) << pos(0), pos(1), pos(2);
    correctedCOMPositions.row(1) << pos(3), pos(4), pos(5);
    return false;




    
    // (com1 - com2)* n - d >= 0
     //dC/dCom1 = n
    // ... = -n
    /**************
     TODO: write position correction procedure:
     1. If the position Constraint is satisfied up to tolerate ("abs(C(p)<=tolerance"), set corrected values to original ones 
     and return true
     
     2. Otherwise, correct COM position as learnt in class.
     Note that since this is a linear correction, correcting COM position == correcting all positions the same offset. 
     the currConstPositions are used to measure the constraint, and the COM values are corrected accordingly to create the effect.
     
     Note to differentiate between different constraint types; for inequality constraints you don't do anything unless 
     it's unsatisfied.
     ***************/
    
  }
};



#endif /* constraints_h */
