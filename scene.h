#ifndef SCENE_HEADER_FILE
#define SCENE_HEADER_FILE

#include <vector>
#include <fstream>
#include <igl/bounding_box.h>
#include <igl/readMESH.h>
#include "ccd.h"
#include "volInt.h"
#include "auxfunctions.h"

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
  Matrix3d invIT;  //Original *inverse* inertia tensor around the COM, defined in the rest state to the object (so to the canonical world system)
  
  VectorXd tetVolumes;    //|T|x1 tetrahedra volumes
  VectorXd invMasses;     //|T|x1 tetrahedra *inverse* masses
  
  //kinematics
  bool isFixed;  //is the object immobile
  double totalMass;  //sum(1/invMass)
  double totalVolume;
  RowVector3d comVelocity;  //the linear velocity of the center of mass
  RowVector3d angVelocity;  //the angular velocity of the object.
  
  //dynamics
  std::vector<Impulse> currImpulses;  //current list of impulses, updated by collision handling
  
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
  //HAS TODO
  Matrix3d getCurrInvInertiaTensor(){
    Matrix3d R=Q2RotMatrix(orientation);
    Matrix3d I_world = R * invIT * R.transpose();

    return I_world; 
  }
  
  
  //Update the current position and orientation by integrating the linear and angular velocities, and update currV accordingly
  //You need to modify this according to its purpose
  //HAS TODO
  void updatePosition(double timeStep){   
    //just forward Euler now
    if (isFixed)
      return;  //a fixed object is immobile
    
    /***************
    TODO
     ***************/

    //Vector3d gravity; gravity << 0, -9.8, 0.0;
    //COM += comVelocity * timeStep;
    //cout << COM << "\n";
    COM += comVelocity * timeStep;
    //cout << COM;

    //// Convert orientation to quaternion
    Quaterniond q(orientation[0], orientation[1], orientation[2], orientation[3]);

    //// Calculate the quaternion representing the rotation due to angular velocity
    Quaterniond deltaQ;
    Vector3d omega = angVelocity * timeStep;
    double theta = omega.norm();
    Vector3d rotation_axis = omega.normalized();
    
    //deltaQ = Quaterniond(1, angVelocity[0] * timeStep / 2, angVelocity[1] * timeStep / 2, angVelocity[2] * timeStep / 2);
    deltaQ = Quaterniond(cos(theta / 2), sin(theta / 2) * rotation_axis.x(), sin(theta / 2) * rotation_axis.y(), sin(theta / 2) * rotation_axis.z());
    
    //// Update orientation by quaternion multiplication
    q = q * deltaQ;

    //// Normalize the quaternion
    q.normalize();

    //// Extract the updated orientation as a 4D vector
    orientation << q.w(), q.x(), q.y(), q.z();
    
    for (int i=0;i<currV.rows();i++)
      currV.row(i)<<QRot(origV.row(i), orientation)+COM;
  }
  
  
  //Updating velocity *instantaneously*. i.e., not integration from acceleration, but as a result of a collision impulse from the "impulses" list
  //You need to modify this for that purpose.
  //HAS TODO
  void updateImpulseVelocities(){
    
    if (isFixed){
      comVelocity.setZero();
      currImpulses.clear();
      angVelocity.setZero();
      return;
    }
    RowVector3d impulses = RowVector3d::Zero();
    RowVector3d POC, R;
    cout << "comVelocity before: " << comVelocity << endl;
    for (int i = 0; i < currImpulses.size(); i++) {
        impulses = currImpulses[i].second;
        comVelocity += (impulses / totalMass);
        cout << "comVelocity after: " << comVelocity << endl;

        POC = currImpulses[i].first;
        //R = POC - COM;
        //norm_R = R;
        //norm_R.normalize();
        //Vector3d R_crossed = R.cross(norm_R);
        R = POC - COM;
        RowVector3d torque = R.cross(impulses);
        cout << "angVelocity before: " << angVelocity << endl;

        angVelocity += getCurrInvInertiaTensor() * torque.transpose();
        //angVelocity += (impulses * getCurrInvInertiaTensor()) * R_crossed; //WHY IS THIS NOT WORKING AGHHHHH
        cout << "angVelocity after: " << angVelocity << endl;
    }

    currImpulses.clear();
    //update linear and angular velocity according to all impulses
    /***************
     TODO
     ***************/
  }
  

  //HAS TODO???
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
  
    return naturalCOM;
    
  }
  
  
  //Integrating the linear and angular velocities of the object
  //You need to modify this to integrate from acceleration in the field (basically gravity)
  void updateVelocity(double timeStep, double dragCoeff){
    if (isFixed)
      return;
    //integrating external forces (only gravity)
    Vector3d gravity; gravity<<0,-9.8,0.0;
    comVelocity+=gravity*timeStep;
    comVelocity -= dragCoeff * comVelocity;
    angVelocity -= dragCoeff * angVelocity;
  }
  
  
  //the full integration for the time step (velocity + position)
  //You need to modify this if you are changing the integration
  void integrate(double timeStep, double dragCoeff){
    updateVelocity(timeStep, dragCoeff);
    updatePosition(timeStep);
  }

  void applyExternalForce(RowVector3d vec) {
      currImpulses.push_back(Impulse(COM, vec));
      updateImpulseVelocities();
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
    
    RowVector3d naturalCOM;  //by the geometry of the object
    
    //initializes the original geometric properties (COM + IT) of the object
    naturalCOM = initStaticProperties(density);
    
    origV.rowwise()-=naturalCOM;  //removing the natural COM of the OFF file (natural COM is never used again)
    
    currV.resize(origV.rows(), origV.cols());
    for (int i=0;i<currV.rows();i++)
      currV.row(i)<<QRot(origV.row(i), orientation)+COM;
    
    
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
  double currTime;
  int numFullV, numFullT;
  std::vector<Mesh> meshes;
  int selectedObj;
  
  //adding an objects. You do not need to update this generally
  void addMesh(const MatrixXd& V, const MatrixXi& F, const MatrixXi& T, const double density, const bool isFixed, const RowVector3d& COM, const RowVector4d& orientation){
    
    Mesh m(V,F, T, density, isFixed, COM, orientation);
    meshes.push_back(m);
  }
  
  /*********************************************************************
   This function handles a collision between objects ro1 and ro2 when found, by assigning impulses to both objects.
   Input: RigidObjects m1, m2
   depth: the depth of penetration
   contactNormal: the normal of the conact measured m1->m2
   penPosition: a point on m2 such that if m2 <= m2 + depth*contactNormal, then penPosition+depth*contactNormal is the common contact point
   CRCoeff: the coefficient of restitution
   *********************************************************************/
  void handleCollision(Mesh& m1, Mesh& m2,const double& depth, const RowVector3d& contactNormal,const RowVector3d& penPosition, const double CRCoeff){
    
    
    std::cout<<"contactNormal: "<<contactNormal<<std::endl;
    std::cout<<"penPosition: "<<penPosition<<std::endl;
    //std::cout<<"handleCollision begin"<<std::endl;
    
    
    //Interpretation resolution: move each object by inverse mass weighting, unless either is fixed, and then move the other. Remember to respect the direction of contactNormal and update penPosition accordingly.
    RowVector3d contactPosition = penPosition + depth * contactNormal;
    std::cout << "contactPosition: " << contactPosition << std::endl;


    double tot_inv_mass;
    RowVector3d moveperINVmass;
    if (m1.isFixed) {
        m2.COM += depth * contactNormal;
        tot_inv_mass = 1 / m2.totalMass;
        // moveperINVmass = depth * contactNormal / tot_inv_mass;
        // RowVector3d correction = moveperINVmass / m2.totalMass;
        // m2.COM += correction;
        

    } else if (m2.isFixed){
        m1.COM -= depth * contactNormal;
        tot_inv_mass = 1 / m1.totalMass;
        //moveperINVmass = depth * contactNormal / tot_inv_mass;
        //RowVector3d correction = -moveperINVmass / m1.totalMass;
        //m1.COM += correction;

    } else { //inverse mass weighting
        tot_inv_mass = (1 / m1.totalMass) + (1 / m2.totalMass);
        double m1Ratio = (1 / m1.totalMass) / tot_inv_mass;
        //moveperINVmass = depth * contactNormal / tot_inv_mass;
        //RowVector3d correctionm1 = -moveperINVmass / m1.totalMass;
        //RowVector3d correctionm2 = moveperINVmass / m2.totalMass;
        m1.COM -= depth*contactNormal * m1Ratio;
        m2.COM += depth * contactNormal * (1- m1Ratio);
    }
    
    std::cout << "m1.COM: " << m1.COM << std::endl;
    std::cout << "m2.COM: " << m2.COM << std::endl;

    double j, upper, lower;
    /*RowVector3d r_a, r_b;
    r_a = m1.COM - contactPosition;
    r_b = m2.COM - contactPosition;

    Vector3d b_without_transpose = r_b.cross(contactNormal);
    RowVector3d b_transposed = b_without_transpose.transpose();
    Matrix3d b_inertia = m2.getCurrInvInertiaTensor();
    Vector3d b_not_transposed_inertia = b_inertia * b_without_transpose;
    double b_dot = b_transposed * b_not_transposed_inertia;

    Vector3d a_without_transpose = r_a.cross(contactNormal);
    RowVector3d a_transposed = a_without_transpose.transpose();
    Matrix3d a_inertia = m1.getCurrInvInertiaTensor();
    Vector3d a_not_transposed_inertia = a_inertia * a_without_transpose;
    double a_dot = a_transposed * a_not_transposed_inertia;
    std::cout << "a_dot: " << a_dot << std::endl;
    std::cout << "b_dot: " << b_dot << std::endl;
    //std::cout << b_dot << std::endl;


    if (m1.isFixed) {
        upper = (1 + CRCoeff) * (-m2.comVelocity.dot(contactNormal));
        lower = (1 / m2.totalMass) + b_dot;
        j = upper / lower;

    }
    else if (m2.isFixed) {
        upper = (1 + CRCoeff) * (m1.comVelocity.dot(contactNormal));
        lower = (1 / m1.totalMass) + a_dot;
        j = -upper / lower;

    }
    else { //inverse mass weighting
    upper = (1 + CRCoeff) * ((m1.comVelocity - m2.comVelocity).dot(contactNormal));
    lower = (1 / m1.totalMass) + (1 / m2.totalMass) + a_dot + b_dot;
    j = upper / lower;
    //}
    */
    RowVector3d r1 = contactPosition - m1.COM;
    RowVector3d r2 = contactPosition - m2.COM;
    RowVector3d relvel = (m2.comVelocity + m2.angVelocity.cross(r2)) - (m1.comVelocity + m1.angVelocity.cross(r1));
    cout << "rad1: " << r1 << endl;
    cout << "rad2: " << r2 << endl;

    double velnormal = relvel.dot(contactNormal);

    //Not sure which raCrossN to use for better results

    if (velnormal > 0) return;
    cout << "closeVelocity: " << velnormal << endl;
    RowVector3d raTrans = r1.cross(contactNormal).transpose();
    //double raCrossN = r1.cross(contactNormal).dot((r1.cross(contactNormal) * m1.getCurrInvInertiaTensor()));
    double raCrossN = r1.cross(contactNormal).dot((raTrans * m1.getCurrInvInertiaTensor()));

    RowVector3d rbTrans = r2.cross(contactNormal).transpose();
    //double rbCrossN = r2.cross(contactNormal).dot((r2.cross(contactNormal) * m2.getCurrInvInertiaTensor()));
    double rbCrossN = r2.cross(contactNormal).dot((rbTrans * m2.getCurrInvInertiaTensor()));
    //double invEffectiveMass = tot_inv_mass + a_dot + b_dot;
    
    double invEffectiveMass = tot_inv_mass + raCrossN + rbCrossN;
    j = -(1 + CRCoeff) * velnormal / invEffectiveMass;
    std::cout << "invEffectiveMass: " << 1 / invEffectiveMass << std::endl;
    //Create impulse and push them into m1.impulses and m2.impulses.
    
    RowVector3d impulse= j * contactNormal;  //changed this to my result
    
    std::cout<<"impulse: "<<impulse<<std::endl;
    if (impulse.norm()>10e-6){
      m1.currImpulses.push_back(Impulse(contactPosition, -impulse));
      m2.currImpulses.push_back(Impulse(contactPosition, impulse));
    }
    
    //std::cout<<"handleCollision end"<<std::endl;
    
    //updating velocities according to impulses
    m1.updateImpulseVelocities();
    m2.updateImpulseVelocities();
  }

  void applyExternalForce(int i,double x, double y, double z) {
      meshes[i].applyExternalForce(RowVector3d(x,y,z));
  }
  /*********************************************************************
   This function handles a single time step by:
   1. Integrating velocities, positions, and orientations by the timeStep
   2. detecting and handling collisions with the coefficient of restitutation CRCoeff
   3. updating the visual scene in fullV and fullT
   *********************************************************************/
  void updateScene(double timeStep, double CRCoeff, double dragCoeff){
    
    //integrating velocity, position and orientation from forces and previous states
    for (int i=0;i<meshes.size();i++)
        meshes[i].integrate(timeStep, dragCoeff);
    
    //detecting and handling collisions when found
    //This is done exhaustively: checking every two objects in the scene.
    double depth;
    RowVector3d contactNormal, penPosition;
    for (int i=0;i<meshes.size();i++)
      for (int j=i+1;j<meshes.size();j++)
        if (meshes[i].isCollide(meshes[j],depth, contactNormal, penPosition))
          handleCollision(meshes[i], meshes[j],depth, contactNormal, penPosition,CRCoeff);
    
    currTime+=timeStep;
  }
  
  //loading a scene from the scene .txt files
  //you do not need to update this function
  bool loadScene(const std::string dataFolder, const std::string sceneFileName){
    
    ifstream sceneFileHandle;
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
      igl::readMESH(dataFolder+std::string("/")+MESHFileName,objV,objT, objF);
      
      //fixing weird orientation problem
      MatrixXi tempF(objF.rows(),3);
      tempF<<objF.col(2), objF.col(1), objF.col(0);
      objF=tempF;
      
      addMesh(objV,objF, objT,density, isFixed, userCOM, userOrientation);
      cout << "COM: " << userCOM <<endl;
      cout << "orientation: " << userOrientation <<endl;
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
