#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"
#include "field/GhostLayerField.h"
#include "stencil/Directions.h"


namespace walberla {

template<typename LM>
class StreamPullCollideGeneric
{
public:
    StreamPullCollideGeneric(BlockDataID src, BlockDataID dst, real_t omega)
        : src_(src), dst_(dst), omega_(omega)
    {}

    void operator()(IBlock * block)
    {
        using PdfField_T = GhostLayerField<real_t, LM::Stencil::Q>;
        auto src = block->getData<PdfField_T>( src_ );
        auto dst = block->getData<PdfField_T>( dst_ );

        WALBERLA_FOR_ALL_CELLS_XYZ(src, {
                    // Compute conserved moments
                    Vector3<real_t> u(0);
                    real_t rho = 0;
                    using stencil::cx;
                    using stencil::cy;
                    using stencil::cz;
                    for( auto d = LM::Stencil::begin(); d != LM::Stencil::end(); ++d )
                    {
                        const real_t pdf = src->get(x - d.cx(), y - d.cy(), z - d.cz(), d.toIdx());
                        u[0] += pdf * d.cx();
                        u[1] += pdf * d.cy();
                        u[2] += pdf * d.cz();
                        rho += pdf;
                    }

                    // collide
                    const real_t vel_sqr = u.sqrLength() * 1.5;
                    for( auto d = LM::Stencil::begin(); d != LM::Stencil::end(); ++d )
                    {
                        const real_t pdf = src->get(x - d.cx(), y - d.cy(), z - d.cz(), d.toIdx());
                        const real_t vel = d.cx()*u[0] + d.cy()*u[1] + d.cz()*u[2];
                        dst->get(x, y, z, d.toIdx()) = ( 1.0 - omega_ ) * pdf +
                                                       omega_   * LM::w[ d.toIdx() ] * ( rho - vel_sqr + 3.0*vel + 4.5*vel*vel );
                    }
                });
        src->swapDataPointers(*dst);
    }

private:
    BlockDataID src_;
    BlockDataID dst_;
    real_t omega_;
};



} // namespace walberla