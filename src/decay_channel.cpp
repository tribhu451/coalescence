#include "decay_channel.h"


decay_channel::decay_channel(double a_branch_ratio, int a_daughter_1, int a_daughter_2, 
int a_daughter_3, int a_daughter_4, int a_daughter_5 ): branching_ratio(a_branch_ratio), daughter_1(a_daughter_1), 
daughter_2(a_daughter_2), daughter_3(a_daughter_3), daughter_4(a_daughter_4), daughter_5(a_daughter_5) {
}

decay_channel::~decay_channel(){
}
