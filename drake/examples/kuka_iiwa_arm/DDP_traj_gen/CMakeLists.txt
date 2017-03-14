#add_library_with_exports(LIB_NAME drakeKukaArmDDPTraj SOURCE_FILES
#    cart_pole.cpp
#    cost_function_cart_pole.cpp
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
    cart_pole.h
    cost_function_cart_pole.h
    ilqrsolver.h
    udpsolver.h)
#drake_install_pkg_config_file(drake-kuka-iiwa-arm-ddp-traj
#    TARGET drakeKukaArmDDPTraj
#    LIBS -drakeKukaArmDDPTraj
#REQUIRES drake-lcm-system drake-lcmtypes-cpp drake-rbm)