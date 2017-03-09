#include "dynamicmodel.h"

// DynamicModel::DynamicModel()
// {
// }
unsigned int DynamicModel::getStateNb()
{
    return stateNb;
}

unsigned int DynamicModel::getCommandNb()
{
    return commandNb;
}

commandVec_t& DynamicModel::getLowerCommandBounds()
{
    return lowerCommandBounds;
}

commandVec_t& DynamicModel::getUpperCommandBounds()
{
    return upperCommandBounds;
}

stateMatTab_t& DynamicModel::getfxList()
{
    return fxList;
}

stateR_commandC_tab_t& DynamicModel::getfuList()
{
    return fuList;
}