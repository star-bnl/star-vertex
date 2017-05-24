#ifndef VertexTrackFilter_h
#define VertexTrackFilter_h

#include "tables/St_VertexCuts_Table.h"


namespace star_vertex {


class VertexTrackFilter
{
public:

   VertexTrackFilter(VertexCuts_st vertexCuts);

private:

   VertexCuts_st  mVertexCuts;
};


}

#endif
