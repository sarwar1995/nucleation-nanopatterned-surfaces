//
//  RectPatch.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 9/25/20.
//  Copyright © 2020 Sarwar Hussain. All rights reserved.
//

#ifndef RectPatch_hpp
#define RectPatch_hpp


/* This can be used when more patche functionality is added instead of just rectangular patch which is currently what the patch class is implementing*/

#include <stdio.h>
#include "Patch.hpp"
class RectPatch:public Patch{
public:
    RectPatch();
    ~RectPatch();

private:
    std::vector<double> dimensions; //0: X-dim 1:Y-dim
};

#endif /* RectPatch_hpp */
