#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <iostream>
#include "scene.h"
#include "line_cylinders.h"


Eigen::MatrixXd V;
Eigen::MatrixXi F;
igl::opengl::glfw::Viewer mgpViewer;

float currTime = 0;
bool animationHack;  //fixing the weird camera bug in libigl
//initial values
float timeStep = 0.02;
float CRCoeff= 0.1;
bool strecth = false;
int frameNum = 1000;

bool withObject = true;

double tolerance = 10e-3;
int maxIterations=10000;

bool offline = true;

bool platVisibility = true;

Scene scene;

Eigen::MatrixXd platV;
Eigen::MatrixXi platF;
Eigen::MatrixXi platT;
Eigen::RowVector3d platCOM;
Eigen::RowVector4d platOrientation;

Eigen::MatrixXd platV1;
Eigen::MatrixXi platF1;
Eigen::MatrixXi platT1;
Eigen::RowVector3d platCOM1;
Eigen::RowVector4d platOrientation1;

Eigen::MatrixXd platV2;
Eigen::MatrixXi platF2;
Eigen::MatrixXi platT2;
Eigen::RowVector3d platCOM2;
Eigen::RowVector4d platOrientation2;

Eigen::MatrixXd platV3;
Eigen::MatrixXi platF3;
Eigen::MatrixXi platT3;
Eigen::RowVector3d platCOM3;
Eigen::RowVector4d platOrientation3;

Eigen::MatrixXd platV4;
Eigen::MatrixXi platF4;
Eigen::MatrixXi platT4;
Eigen::RowVector3d platCOM4;
Eigen::RowVector4d platOrientation4;

Eigen::MatrixXd platV5;
Eigen::MatrixXi platF5;
Eigen::MatrixXi platT5;
Eigen::RowVector3d platCOM5;
Eigen::RowVector4d platOrientation5;

vector<vector<RowVector3d>> dataArray;

void createPlatform()
{
  double platWidth=270.0;
  platCOM<<0.0, -platWidth+platWidth/20.0, 0.0;
  platV.resize(9,3);
  platF.resize(12,3);
  platT.resize(12,4);
  platV<<-platWidth,0.0,-platWidth,
  -platWidth,0.0,platWidth,
  platWidth,0.0,platWidth,
  platWidth,0.0, -platWidth,
  -platWidth,-platWidth/10.0,-platWidth,
  -platWidth,-platWidth/10.0,platWidth,
  platWidth,-platWidth/10.0,platWidth,
  platWidth,-platWidth/10.0, -platWidth,
  0.0,-platWidth/20.0, 0.0;
  platF<<0,1,2,
  2,3,0,
  6,5,4,
  4,7,6,
  1,0,5,
  0,4,5,
  2,1,6,
  1,5,6,
  3,2,7,
  2,6,7,
  0,3,4,
  3,7,4;
  
  platOrientation<<1.0,0.0,0.0,0.0;
  
  platT<<platF, VectorXi::Constant(12,8);


  platCOM1 << platWidth- platWidth / 20.0, platWidth / 20.0, 0.0;
  //platCOM1 << platWidth, platWidth / 20.0, 0.0;
  platV1.resize(9, 3);
  platF1.resize(12, 3);
  platT1.resize(12, 4);
  platV1 << -platWidth, 0.0, -platWidth,
      -platWidth, 0.0, platWidth,
      platWidth, 0.0, platWidth,
      platWidth, 0.0, -platWidth,
      -platWidth, -platWidth / 10.0, -platWidth,
      -platWidth, -platWidth / 10.0, platWidth,
      platWidth, -platWidth / 10.0, platWidth,
      platWidth, -platWidth / 10.0, -platWidth,
      0.0, -platWidth / 20.0, 0.0;
  platF1 << 0, 1, 2,
      2, 3, 0,
      6, 5, 4,
      4, 7, 6,
      1, 0, 5,
      0, 4, 5,
      2, 1, 6,
      1, 5, 6,
      3, 2, 7,
      2, 6, 7,
      0, 3, 4,
      3, 7, 4;

  platOrientation1 << 0.0, 1.0, 1.0, 0.0;

  platT1 << platF1, VectorXi::Constant(12, 8);



  platCOM2 << 0.0, platWidth / 20.0, platWidth- platWidth / 20.0;
  //platCOM2 << 0.0, platWidth / 20.0, platWidth;
  platV2.resize(9, 3);
  platF2.resize(12, 3);
  platT2.resize(12, 4);
  platV2 << -platWidth, 0.0, -platWidth,
      -platWidth, 0.0, platWidth,
      platWidth, 0.0, platWidth,
      platWidth, 0.0, -platWidth,
      -platWidth, -platWidth / 10.0, -platWidth,
      -platWidth, -platWidth / 10.0, platWidth,
      platWidth, -platWidth / 10.0, platWidth,
      platWidth, -platWidth / 10.0, -platWidth,
      0.0, -platWidth / 20.0, 0.0;
  platF2 << 0, 1, 2,
      2, 3, 0,
      6, 5, 4,
      4, 7, 6,
      1, 0, 5,
      0, 4, 5,
      2, 1, 6,
      1, 5, 6,
      3, 2, 7,
      2, 6, 7,
      0, 3, 4,
      3, 7, 4;

  platOrientation2 << 0.0, 0.0, 1.0, 1.0;

  platT2 << platF2, VectorXi::Constant(12, 8);



  platCOM3 << -platWidth+ platWidth / 20.0, platWidth / 20.0, 0.0;
  //platCOM3 << -platWidth, platWidth / 20.0, 0.0;
  platV3.resize(9, 3);
  platF3.resize(12, 3);
  platT3.resize(12, 4);
  platV3 << -platWidth, 0.0, -platWidth,
      -platWidth, 0.0, platWidth,
      platWidth, 0.0, platWidth,
      platWidth, 0.0, -platWidth,
      -platWidth, -platWidth / 10.0, -platWidth,
      -platWidth, -platWidth / 10.0, platWidth,
      platWidth, -platWidth / 10.0, platWidth,
      platWidth, -platWidth / 10.0, -platWidth,
      0.0, -platWidth / 20.0, 0.0;
  platF3 << 0, 1, 2,
      2, 3, 0,
      6, 5, 4,
      4, 7, 6,
      1, 0, 5,
      0, 4, 5,
      2, 1, 6,
      1, 5, 6,
      3, 2, 7,
      2, 6, 7,
      0, 3, 4,
      3, 7, 4;

  platOrientation3 << 0.0, 1.0, 1.0, 0.0;

  platT3 << platF3, VectorXi::Constant(12, 8);

  platCOM4 << 0.0, platWidth / 20.0, -platWidth + platWidth / 20.0;
  //platCOM4 << 0.0, platWidth / 20.0, -platWidth;
  platV4.resize(9, 3);
  platF4.resize(12, 3);
  platT4.resize(12, 4);
  platV4 << -platWidth, 0.0, -platWidth,
      -platWidth, 0.0, platWidth,
      platWidth, 0.0, platWidth,
      platWidth, 0.0, -platWidth,
      -platWidth, -platWidth / 10.0, -platWidth,
      -platWidth, -platWidth / 10.0, platWidth,
      platWidth, -platWidth / 10.0, platWidth,
      platWidth, -platWidth / 10.0, -platWidth,
      0.0, -platWidth / 20.0, 0.0;
  platF4 << 0, 1, 2,
      2, 3, 0,
      6, 5, 4,
      4, 7, 6,
      1, 0, 5,
      0, 4, 5,
      2, 1, 6,
      1, 5, 6,
      3, 2, 7,
      2, 6, 7,
      0, 3, 4,
      3, 7, 4;

  platOrientation4 << 0.0, 0.0, 1.0, 1.0;

  platT4 << platF4, VectorXi::Constant(12, 8);


  int cubeWidth = 70;
  platCOM5 << -50.0,-platWidth + platWidth/20 + 60, -100.0;
  //platCOM4 << 0.0, platWidth / 20.0, -platWidth;
  platV5.resize(9, 3);
  platF5.resize(12, 3);
  platT5.resize(12, 4);
  platV5 << -cubeWidth, 0.0, -cubeWidth,
      -cubeWidth, 0.0, cubeWidth,
      cubeWidth, 0.0, cubeWidth,
      cubeWidth, 0.0, -cubeWidth,
      -cubeWidth, -cubeWidth , -cubeWidth,
      -cubeWidth, -cubeWidth , cubeWidth,
      cubeWidth, -cubeWidth , cubeWidth,
      cubeWidth, -cubeWidth , -cubeWidth,
      0.0, -cubeWidth , 0.0;
  platF5 << 0, 1, 2,
      2, 3, 0,
      6, 5, 4,
      4, 7, 6,
      1, 0, 5,
      0, 4, 5,
      2, 1, 6,
      1, 5, 6,
      3, 2, 7,
      2, 6, 7,
      0, 3, 4,
      3, 7, 4;

  platOrientation5 << 0.0, 0.0, 1.8, 0.6;

  platT5 << platF5, VectorXi::Constant(12, 8);
  
  
}

void updateMeshes(igl::opengl::glfw::Viewer &viewer)
{
  RowVector3d platColor; platColor<<0.8,0.8,0.8;
  RowVector3d meshColor; meshColor<<0.8,0.2,0.2;
  RowVector3d bluemeshColor; bluemeshColor << 0.2, 0.2, 0.8;

  //RowVector3d transparentColor; transparentColor << 0.8, 0.2, 0.2, 0.0;
  viewer.core().align_camera_center(scene.meshes[0].currV);
  for (int i=0;i<scene.meshes.size();i++){
    viewer.data_list[i].clear();
    viewer.data_list[i].set_mesh(scene.meshes[i].currV, scene.meshes[i].F);
    viewer.data_list[i].set_face_based(true);
    if (i > 5) {
        viewer.data_list[i].set_colors(bluemeshColor);
    }
    else {
        viewer.data_list[i].set_colors(meshColor);
    }
    viewer.data_list[i].show_lines=false;
  }
  viewer.data_list[0].show_lines=false;
  viewer.data_list[0].set_colors(platColor.replicate(scene.meshes[0].F.rows(),1));
  //viewer.data_list[0].set_colors(transparentColor);
  viewer.data_list[0].set_face_based(true);
  viewer.data_list[0].set_visible(platVisibility);

  viewer.data_list[1].show_lines = false;
  viewer.data_list[1].set_colors(platColor.replicate(scene.meshes[1].F.rows(), 1));
  viewer.data_list[1].set_face_based(true);
  viewer.data_list[1].set_visible(platVisibility);

  viewer.data_list[2].show_lines = false;
  viewer.data_list[2].set_colors(platColor.replicate(scene.meshes[2].F.rows(), 1));
  viewer.data_list[2].set_face_based(true);
  viewer.data_list[2].set_visible(platVisibility);

  viewer.data_list[3].show_lines = false;
  viewer.data_list[3].set_colors(platColor.replicate(scene.meshes[3].F.rows(), 1));
  viewer.data_list[3].set_face_based(true);
  viewer.data_list[3].set_visible(false);

  viewer.data_list[4].show_lines = false;
  viewer.data_list[4].set_colors(platColor.replicate(scene.meshes[3].F.rows(), 1));
  viewer.data_list[4].set_face_based(true);
  viewer.data_list[4].set_visible(false);

  if (withObject == true) {
      viewer.data_list[5].show_lines = false;
      viewer.data_list[5].set_colors(meshColor);
      viewer.data_list[5].set_face_based(true);
      viewer.data_list[5].set_visible(true);
  }
  //viewer.core.align_camera_center(scene.meshes[0].currV);
  
  //updating constraint viewing
  MatrixXi constF;
  MatrixXd constV, constC;
  
  MatrixXd P1(scene.constraints.size(),3);
  MatrixXd P2(scene.constraints.size(),3);
  for (int i=0;i<scene.constraints.size();i++){
    P1.row(i)=scene.meshes[scene.constraints[i].m1].currV.row(scene.constraints[i].v1);
    P2.row(i)=scene.meshes[scene.constraints[i].m2].currV.row(scene.constraints[i].v2);
  }
  
  MatrixXd cyndColors=RowVector3d(1.0,1.0,0.0).replicate(P1.size(),1);
  
  double radius = 0.5;
  hedra::line_cylinders(P1,P2,radius,cyndColors,8,constV,constF,constC);
  viewer.data_list[scene.meshes.size()].set_mesh(constV, constF);
  viewer.data_list[scene.meshes.size()].set_face_based(true);
  viewer.data_list[scene.meshes.size()].set_colors(constC);
  viewer.data_list[scene.meshes.size()].show_lines=false;
  
  
}


void updateMeshesOffline(igl::opengl::glfw::Viewer& viewer)
{
    RowVector3d platColor; platColor << 0.8, 0.8, 0.8;
    RowVector3d meshColor; meshColor << 0.8, 0.2, 0.2;
    RowVector3d bluemeshColor; bluemeshColor << 0.2, 0.2, 0.8;

    omp_set_num_threads(omp_get_max_threads());

    int index = currTime / timeStep;

#pragma omp parallel for
    for (int i = 5; i < scene.meshes.size(); i++) {
        for (int j = 0; j < scene.meshes[i].currV.rows(); j++)
        {
            scene.meshes[i].currV.row(j) << scene.meshes[i].origV.row(j) + dataArray[index][i];

            /*if (scene.meshes[i].density > 110 && scene.meshes[i].density < 112) {
                RowVector3d normAng = {0,0.01,0};
                RowVector4d q(0, normAng.x(), normAng.y(), normAng.z());
                scene.meshes[i].orientation += 0.5 * timeStep*4 * QMult(q, scene.meshes[i].orientation);

                scene.meshes[i].currV.row(j) << QRot(scene.meshes[i].origV.row(j), scene.meshes[i].orientation) + dataArray[index][i];
            }*/
        }
    }
#pragma omp barrier

    int bound;

    if (withObject == true) {
        bound = 5;
    }
    else {
        bound = 4;
    }

    viewer.core().align_camera_center(scene.meshes[0].currV);
    for (int i = 0; i < scene.meshes.size(); i++) {
        viewer.data_list[i].clear();
        viewer.data_list[i].set_mesh(scene.meshes[i].currV, scene.meshes[i].F);
        viewer.data_list[i].set_mesh(scene.meshes[i].currV, scene.meshes[i].F);
        viewer.data_list[i].set_face_based(true);

        if (i > bound) {
            viewer.data_list[i].set_colors(bluemeshColor);
        }
        else {
            viewer.data_list[i].set_colors(meshColor);
        }
        viewer.data_list[i].show_lines = false;
    }
    viewer.data_list[0].show_lines = false;
    viewer.data_list[0].set_colors(platColor.replicate(scene.meshes[0].F.rows(), 1));
    //viewer.data_list[0].set_colors(transparentColor);
    viewer.data_list[0].set_face_based(true);
    viewer.data_list[0].set_visible(platVisibility);

    viewer.data_list[1].show_lines = false;
    viewer.data_list[1].set_colors(platColor.replicate(scene.meshes[1].F.rows(), 1));
    viewer.data_list[1].set_face_based(true);
    viewer.data_list[1].set_visible(platVisibility);

    viewer.data_list[2].show_lines = false;
    viewer.data_list[2].set_colors(platColor.replicate(scene.meshes[2].F.rows(), 1));
    viewer.data_list[2].set_face_based(true);
    viewer.data_list[2].set_visible(platVisibility);

    viewer.data_list[3].show_lines = false;
    viewer.data_list[3].set_colors(platColor.replicate(scene.meshes[3].F.rows(), 1));
    viewer.data_list[3].set_face_based(true);
    viewer.data_list[3].set_visible(false);

    viewer.data_list[4].show_lines = false;
    viewer.data_list[4].set_colors(platColor.replicate(scene.meshes[4].F.rows(), 1));
    viewer.data_list[4].set_face_based(true);
    viewer.data_list[4].set_visible(false);

    if (withObject == true) {
        viewer.data_list[5].show_lines = false;
        viewer.data_list[5].set_colors(meshColor);
        viewer.data_list[5].set_face_based(true);
        viewer.data_list[5].set_visible(true);
    }
    //viewer.core.align_camera_center(scene.meshes[0].currV);

    //updating constraint viewing
    MatrixXi constF;
    MatrixXd constV, constC;

    MatrixXd P1(scene.constraints.size(), 3);
    MatrixXd P2(scene.constraints.size(), 3);
    for (int i = 0; i < scene.constraints.size(); i++) {
        P1.row(i) = scene.meshes[scene.constraints[i].m1].currV.row(scene.constraints[i].v1);
        P2.row(i) = scene.meshes[scene.constraints[i].m2].currV.row(scene.constraints[i].v2);
    }

    MatrixXd cyndColors = RowVector3d(1.0, 1.0, 0.0).replicate(P1.size(), 1);

    double radius = 0.5;
    hedra::line_cylinders(P1, P2, radius, cyndColors, 8, constV, constF, constC);
    viewer.data_list[scene.meshes.size()].set_mesh(constV, constF);
    viewer.data_list[scene.meshes.size()].set_face_based(true);
    viewer.data_list[scene.meshes.size()].set_colors(constC);
    viewer.data_list[scene.meshes.size()].show_lines = false;


}


bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier)
{
  if (key == ' ')
  {
    viewer.core().is_animating = !viewer.core().is_animating;
    if (viewer.core().is_animating)
      cout<<"Simulation running"<<endl;
    else
      cout<<"Simulation paused"<<endl;
    return true;
  }
  
  if (key == 'S')
  {
    if (!viewer.core().is_animating){
      scene.updateScene(timeStep, CRCoeff, tolerance, maxIterations, dataArray);
      currTime+=timeStep;
      updateMeshes(viewer);
      std::cout <<"currTime: "<<currTime<<std::endl;
      return true;
    }
  }

  return false;
}


bool pre_draw(igl::opengl::glfw::Viewer &viewer)
{
  using namespace Eigen;
  using namespace std;
  
  if (offline == false) {
      if (viewer.core().is_animating) {
          if (!animationHack)
              scene.updateScene(timeStep, CRCoeff, tolerance, maxIterations, dataArray);
          else
              viewer.core().is_animating = false;
          animationHack = false;
          currTime += timeStep;
          //cout <<"currTime: "<<currTime<<endl;
          updateMeshes(viewer);
      }
  }
  else {
      if (viewer.core().is_animating) {
          if ( currTime / timeStep >= frameNum-6) {
              currTime = 0;
          }
          else {
              currTime += timeStep;
              //currTime += 4 * timeStep;
              //cout << currTime << endl;
          }
          //cout <<"currTime: "<<currTime<<endl;
          updateMeshesOffline(viewer);
      }
  }
 
  
  return false;
}

class CustomMenu : public igl::opengl::glfw::imgui::ImGuiMenu
{
  
  virtual void draw_viewer_menu() override
  {
    // Draw parent menu
    //ImGuiMenu::draw_viewer_menu();
    
    // Add new group
    if (ImGui::CollapsingHeader("Algorithm Options", ImGuiTreeNodeFlags_DefaultOpen))
    {
      ImGui::InputFloat("CR Coeff",&CRCoeff,0,0,"%.2f");
      
      
      if (ImGui::InputFloat("Time Step", &timeStep)) {
        mgpViewer.core().animation_max_fps = (((int)1.0/timeStep));
      }
    }
  }
};



int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  
  
  // Load scene
  if (argc<4){
    cout<<"Please provide path (argument 1), name of scene file (argument 2), and name of constraints file (argument 3)!"<<endl;
    return 0;
  }
  cout<<"scene file: "<<std::string(argv[2])<<endl;
  //create platform
  createPlatform();
  scene.addMesh(platV, platF, platT, 10000.0, true, platCOM, platOrientation);
  scene.addMesh(platV1, platF1, platT1, 10000.0, true, platCOM1, platOrientation1);
  scene.addMesh(platV2, platF2, platT2, 10000.0, true, platCOM2, platOrientation2);
  scene.addMesh(platV3, platF3, platT3, 10000.0, true, platCOM3, platOrientation3);
  scene.addMesh(platV4, platF4, platT4, 10000.0, true, platCOM4, platOrientation4);

  if(withObject == true)
    scene.addMesh(platV5, platF5, platT5, 111.0, true, platCOM5, platOrientation5);
  
  //load scene from file
  scene.loadScene(std::string(argv[1]),std::string(argv[2]),std::string(argv[3]), strecth);

  

  dataArray = vector<vector<RowVector3d>>(frameNum, vector<RowVector3d>(scene.meshes.size(), { 0,0,0 }));

  scene.updateScene(0.0, CRCoeff, tolerance, maxIterations, dataArray);

  if (offline == true) {
      for (int i = 0; i < frameNum; i++) {
          scene.updateScene(timeStep, CRCoeff, tolerance, maxIterations, dataArray);
          currTime += timeStep;
          cout << currTime << endl;
      }

      currTime = 0;
  }

  // Viewer Settings
  for (int i=0;i<scene.meshes.size();i++){
    if (i!=0)
      mgpViewer.append_mesh();
    //mgpViewer.data_list[i].set_mesh(scene.meshes[i].currV, scene.meshes[i].F);
  }
  //mgpViewer.core.align_camera_center(scene.meshes[0].currV);
  
  //constraints mesh (for lines)
  mgpViewer.append_mesh();
  mgpViewer.callback_pre_draw = &pre_draw;
  mgpViewer.callback_key_down = &key_down;
  mgpViewer.core().is_animating = true;
  animationHack = true;
  mgpViewer.core().animation_max_fps = 50.;
  
  CustomMenu menu;
  igl::opengl::glfw::imgui::ImGuiPlugin plugin;
  mgpViewer.plugins.push_back(&plugin);
  plugin.widgets.push_back(&menu);
  
  cout<<"Press [space] to toggle continuous simulation" << endl;
  cout<<"Press 'S' to advance time step-by-step"<<endl;
  
  updateMeshes(mgpViewer);
  mgpViewer.launch();
 
}
