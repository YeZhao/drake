#ifndef CONFIG_H
#define CONFIG_H

#include <Eigen/Dense>
#include <Eigen/StdVector>

EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::MatrixXd)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::VectorXd)

#define stateSize 4
#define commandSize 1
#define fullstatecommandSize 5
#define TimeHorizon 0.5
#define TimeStep 0.01

// typedef for stateSize types
typedef Eigen::Matrix<double,stateSize,1> stateVec_t;                       // stateSize x 1
typedef Eigen::Matrix<double,1,stateSize> stateVecTrans_t;                  // 1 x stateSize
typedef Eigen::Matrix<double,stateSize,stateSize> stateMat_t;               // stateSize x stateSize
typedef Eigen::Matrix<double,stateSize,stateSize> stateTens_t[stateSize];   // stateSize x stateSize x stateSize

// typedef for commandSize types
typedef Eigen::Matrix<double,commandSize,1> commandVec_t;                           // commandSize x 1
typedef Eigen::Matrix<double,1,commandSize> commandVecTrans_t;                      // 1 x commandSize
typedef Eigen::Matrix<double,commandSize,commandSize> commandMat_t;                 // commandSize x commandSize
typedef Eigen::Matrix<double,commandSize,commandSize> commandTens_t[commandSize];   // stateSize x commandSize x commandSize

// typedef for mixed stateSize and commandSize types
typedef Eigen::Matrix<double,stateSize,commandSize> stateR_commandC_t;                          // stateSize x commandSize
typedef Eigen::Matrix<double,stateSize,commandSize> stateR_commandC_stateD_t[stateSize];        // stateSize x commandSize x stateSize
typedef Eigen::Matrix<double,stateSize,commandSize> stateR_commandC_commandD_t[commandSize];    // stateSize x commandSize x commandSize
typedef Eigen::Matrix<double,commandSize,stateSize> commandR_stateC_t;                          // commandSize x stateSize
typedef Eigen::Matrix<double,commandSize,stateSize> commandR_stateC_stateD_t[stateSize];        // commandSize x stateSize x stateSize
typedef Eigen::Matrix<double,commandSize,stateSize> commandR_stateC_commandD_t[commandSize];    // commandSize x stateSize x commandSize
typedef Eigen::Matrix<double,stateSize,stateSize> stateR_stateC_commandD_t[commandSize];        // stateSize x stateSize x commandSize
typedef Eigen::Matrix<double,commandSize,commandSize> commandR_commandC_stateD_t[stateSize];    // commandSize x commandSize x stateSize
typedef Eigen::Matrix<double,stateSize+commandSize,1> stateAug_t;                          // stateSize + commandSize x 1

// typedef for half commandSize and stateSize types
typedef Eigen::Matrix<double,stateSize/2,1> stateVec_half_t;                                    // stateSize/2 x 1
typedef Eigen::Matrix<double,stateSize/2,stateSize/2> stateMat_half_t;                          // stateSize/2 x stateSize/2
typedef Eigen::Matrix<double,stateSize/2,1> stateVec_half_t;                                    // stateSize/2 x 1
typedef Eigen::Matrix<double,stateSize/2,commandSize> stateR_half_commandC_t;                   // stateSize/2 x commandSize

// typedef for vectorized state and command matrix (over the horizon)
typedef std::vector<stateVec_t> stateVecTab_t;
typedef std::vector<double> costVecTab_t;
typedef std::vector<commandVec_t> commandVecTab_t;
typedef std::vector<stateMat_t> stateMatTab_t;
typedef std::vector<commandMat_t> commandMatTab_t;
typedef std::vector<stateR_commandC_t> stateR_commandC_tab_t;
typedef std::vector<commandR_stateC_t> commandR_stateC_tab_t;

//typedef std::vector<stateTens_t> stateTensTab_t;
typedef std::vector<std::vector<stateMat_t> > stateTensTab_t;
typedef std::vector<std::vector<stateR_commandC_t> > stateR_commandC_Tens_t;

#endif // CONFIG_H