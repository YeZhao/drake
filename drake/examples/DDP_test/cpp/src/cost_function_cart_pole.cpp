#include "cost_function_cart_pole.h"

CostFunctionCartPole::CostFunctionCartPole()
{
    Q << .1,0.0,0.0,0.0,
         0.0,0.1,0.0,0.0,
         0.0,0.0,0.1,0.0,
         0.0,0.0,0.0,0.1;
    Qf << 1000.0,0.0,0.0,0.0,
          0.0,1000,0.0,0.0,
          0.0,0.0,1000,0.0,
          0.0,0.0,0.0,1000;
    R << 0.01;
}