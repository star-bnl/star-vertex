#ifndef TDecayVertexContainers_h
#define TDecayVertexContainers_h

#include <vector>

#include "StVertexRootIO/TDecayVertex.h"
#include "StVertexRootIO/TPrimaryVertex.h"


class TDecayVertexVec : public TObject
{
public:

   TDecayVertexVec() : TObject(), mVertices() {}

   void Add(const TDecayVertex& vertex)
   {
      auto insert_it = std::inserter(mVertices, mVertices.end());
      insert_it = vertex;
   }

   std::vector<TDecayVertex>  mVertices;

   TPrimaryVertex mPrimaryVertex;


   ClassDef(TDecayVertexVec, 1)
};


#endif
