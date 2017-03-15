#add_library_with_exports(LIB_NAME drakeKukaArmDDPTraj SOURCE_FILES
#    kuka_arm.cpp
#    cost_function_kuka_arm.cpp
#    ilqrsolver.cpp
#    udpsolver.cpp)
#target_link_libraries(drakeKukaArmDDPTraj
#    drakeCommon
#    drakeIK
#    drakeRBM
#    drakeMultibodyParsers
#    drakeTrajectories)
#drake_install_libraries(drakeKukaArmDDPTraj)
drake_install_headers(
    kuka_arm.h
    cost_function_kuka_arm.h
    ilqrsolver.h
    udpsolver.h)
#drake_install_pkg_config_file(drake-kuka-iiwa-arm-ddp-traj
#    TARGET drakeKukaArmDDPTraj
#    LIBS -drakeKukaArmDDPTraj
#REQUIRES drake-lcm-system drake-lcmtypes-cpp drake-rbm)