#include "matrixUtil.h"

bool isNan(commandVec_t& xVec){
   for (int row(0); row < xVec.rows(); ++row){
        if (isnan(xVec(row))){
            return true;
        }
    }
   return false;
}

bool isNan(const commandVec_t& xVec){
   for (int row(0); row < xVec.rows(); ++row){
        if (isnan(xVec(row))){
            return true;
        }
    }
   return false;
}