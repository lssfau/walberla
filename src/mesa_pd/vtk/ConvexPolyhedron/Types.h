#include <OpenMesh/Core/Utils/Property.hh>
#include <OpenMesh/Core/Utils/PropertyManager.hh>

namespace walberla {
namespace mesa_pd {

typedef typename OpenMesh::FPropHandleT<size_t> ParticleIdxFacePropertyHandle;
typedef typename OpenMesh::VPropHandleT<size_t> ParticleIdxVertexPropertyHandle;
template <typename MeshType>
using ParticleIdxFacePropertyManager = OpenMesh::PropertyManager<ParticleIdxFacePropertyHandle, MeshType>;
template <typename MeshType>
using ParticleIdxVertexPropertyManager = typename OpenMesh::PropertyManager<ParticleIdxVertexPropertyHandle, MeshType>;


}
}