#include <string>
#include <vector>
#include <algorithm>

#include "drake/multibody/rigid_body_tree.h"
#include "drake/multibody/constraint/rigid_body_constraint.h"
#include "parse_json.h"


namespace drake {
namespace examples {
namespace kuka_iiwa_arm {
namespace {

std::vector<std::unique_ptr<RigidBodyConstraint>> parse_constraints(const std::string&, RigidBodyTree<double>*);
std::unique_ptr<RigidBodyConstraint> parse_constraint(const std::string&, RigidBodyTree<double>*);
std::unique_ptr<PostureConstraint> parse_posture_constraint(const std::string&, RigidBodyTree<double>*);
std::unique_ptr<WorldPositionConstraint> parse_world_position_constraint(const std::string&, RigidBodyTree<double>*);
std::unique_ptr<WorldQuatConstraint> parse_world_quat_constraint(const std::string&, RigidBodyTree<double>*);


std::vector<std::unique_ptr<RigidBodyConstraint>> parse_constraints(const std::string &data, RigidBodyTree<double>* tree){
  auto json_constraint_list = parse_json_list(data);
  std::vector<std::unique_ptr<RigidBodyConstraint>> constraint_vec;
  for (int i=0; i<json_constraint_list.size(); i++){
    // wrap in a try-catch block for some of the floating base constraints
    try{
      constraint_vec.push_back(parse_constraint(json_constraint_list[i],tree));  
    }catch(std::runtime_error &e){
      std::cout << "Warning! Error thrown when constructing constraint: ";
      std::cout << e.what() << std::endl;
    }
  }
  return constraint_vec;
}

std::unique_ptr<RigidBodyConstraint> parse_constraint(const std::string &data, RigidBodyTree<double>* tree){
  
  if (contains(data, "PostureConstraint")){
    return parse_posture_constraint(data, tree);
  }else if (contains(data, "PositionConstraint")){
    return parse_world_position_constraint(data, tree);
  }else if (contains(data, "QuatConstraint")){
    return parse_world_quat_constraint(data, tree);
  }
  std::cout << "Got an unknown constraint type: " << data << std::endl; 
  throw("Constraint not implemented");
}

std::unique_ptr<PostureConstraint> parse_posture_constraint(const std::string &data, RigidBodyTree<double>* tree){
  
  auto constraint = parse_json_object(data);
  
  // initialize defaults
  // tspan
  auto tspan = DrakeRigidBodyConstraint::default_tspan;
  auto tspan_it = constraint.find("tspan");
  if (tspan_it != constraint.end())
    tspan = json_to_eigen_d<2>(tspan_it->second);

  // joint inidices
  Eigen::VectorXi joint_idx;
  auto joint_name_it = constraint.find("joints");
  if (joint_name_it != constraint.end()){
    auto joint_names = parse_json_list(joint_name_it->second);
    for (int i=0; i<joint_names.size(); i++){
      joint_idx << tree->FindIndexOfChildBodyOfJoint(joint_names[i]);
    }
  }else{
    // bad constraint: no joint names
    throw(1);
  }

  // TODO: finish constructing the posture constraint
  
  // for debugging
  // std::cout << "Creating PostureConstraint: " << data << std::endl; 
  // std::cout << "Timespan: \n" << extract_array(data,"tspan") << std::endl;
  // std::cout << "Joint Indices: \n" << joint_idx << std::endl;
  
  auto pc1 = std::make_unique<PostureConstraint>(tree, tspan);
  return pc1;
}

std::unique_ptr<WorldPositionConstraint> parse_world_position_constraint(const std::string &data, RigidBodyTree<double>* tree){
  
  auto constraint = parse_json_object(data);

  // initialize defaults
  // tspan
  auto tspan = DrakeRigidBodyConstraint::default_tspan;
  auto tspan_it = constraint.find("tspan");
  if (tspan_it != constraint.end())
    tspan = json_to_eigen_d<2>(tspan_it->second);

  std::cout << "Decoding the reference frame" << std::endl;
  // referenceFrame
  Eigen::Vector3d reference_pos = Eigen::Vector3d::Zero();
  auto reference_frame_it = constraint.find("referenceFrame");
  if (reference_frame_it != constraint.end()){
    std::cout << reference_frame_it->second << std::endl;
    auto reference_frame = parse_json_object(reference_frame_it->second);
    auto pos_it = reference_frame.find("position");
    if (pos_it != reference_frame.end()){
      std::cout << "Found position: " << pos_it->second << std::endl;
      reference_pos = json_to_eigen_d<3>(pos_it->second);
    }
  } 

  // bounds
  Eigen::Vector3d pos_lb = Eigen::Vector3d::Zero();
  Eigen::Vector3d pos_ub = Eigen::Vector3d::Zero();
  auto pos_lb_it = constraint.find("lowerBound");
  auto pos_ub_it = constraint.find("upperBound");
  if (pos_lb_it != constraint.end())
    pos_lb = json_to_eigen_d<3>(pos_lb_it->second);
  if (pos_ub_it != constraint.end())
    pos_ub = json_to_eigen_d<3>(pos_ub_it->second);
  pos_lb += reference_pos;
  pos_ub += reference_pos;

  // link_name
  std::string link_name = "default_link";
  auto link_name_it = constraint.find("linkName");
  if (link_name_it != constraint.end())
    link_name = link_name_it->second;
  
  // print everything for debugging
  std::cout << "Creating WorldPositionConstraint: " << data << std::endl;
  std::cout << "Timespan: \n" << tspan << std::endl;
  std::cout << "Lower Bound: \n" << pos_lb << std::endl;
  std::cout << "Upper Bound: \n" << pos_ub << std::endl;
  std::cout << "Link Name: \n" << link_name << std::endl;

  // create world position constraint
  auto wpc = std::make_unique<WorldPositionConstraint>(tree, tree->FindBodyIndex(link_name),
                              Eigen::Vector3d::Zero(), pos_lb, pos_ub, tspan);
  return wpc;
}

std::unique_ptr<WorldQuatConstraint> parse_world_quat_constraint(const std::string &data, RigidBodyTree<double>* tree){
  auto constraint = parse_json_object(data);

  // initialize defaults
  // tspan
  auto tspan = DrakeRigidBodyConstraint::default_tspan;
  auto tspan_it = constraint.find("tspan");
  if (tspan_it != constraint.end())
    tspan = json_to_eigen_d<2>(tspan_it->second);

  // quat_des
  Eigen::Vector4d quat_des = Eigen::Vector4d::Zero();
  auto quat_it = constraint.find("quaternion");
  if (quat_it != constraint.end()){
    auto quat = parse_json_object(quat_it->second);
    auto quat_des_it = quat.find("quaternion");
    if(quat_des_it != quat.end())
      quat_des = json_to_eigen_d<4>(quat_des_it->second);
  }

  // tol
  double tol = 0.1;

  // link_name
  std::string link_name = "default_link";
  auto link_name_it = constraint.find("linkName");
  if (link_name_it != constraint.end())
    link_name = link_name_it->second;

  // print everything for debugging
  std::cout << "Creating WorldQuatConstraint: " << data << std::endl;
  std::cout << "Timespan: \n" << tspan << std::endl;
  std::cout << "Tolerance: \n" << tol << std::endl;
  std::cout << "Quaternion: \n" << quat_des << std::endl;
  std::cout << "Link Name: \n" << link_name << std::endl;
  
  // create world quat constraint
  auto wqc = std::make_unique<WorldQuatConstraint>(tree, tree->FindBodyIndex(link_name),
                              quat_des, tol, Eigen::Vector2d(1, 3));
  return wqc;
}

} // anonymous
} // kuka_iiwa_arm
} // examples
} // drake
