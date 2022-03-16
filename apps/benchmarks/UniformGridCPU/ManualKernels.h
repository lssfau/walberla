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
                        u[0] += pdf * real_c(d.cx());
                        u[1] += pdf * real_c(d.cy());
                        u[2] += pdf * real_c(d.cz());
                        rho += pdf;
                    }

                    // collide
                    const real_t vel_sqr = u.sqrLength() * real_t(1.5);
                    for( auto d = LM::Stencil::begin(); d != LM::Stencil::end(); ++d )
                    {
                        const real_t pdf = src->get(x - d.cx(), y - d.cy(), z - d.cz(), d.toIdx());
                        const real_t vel = real_c(d.cx())*u[0] + real_c(d.cy())*u[1] + real_c(d.cz())*u[2];
                        dst->get(x, y, z, d.toIdx()) = ( real_t(1.0) - omega_ ) * pdf +
                                                       omega_ * LM::w[ d.toIdx() ] * ( rho - vel_sqr + real_t(3.0)*vel + real_t(4.5)*vel*vel );
                    }
                });
        src->swapDataPointers(*dst);
    }

private:
    BlockDataID src_;
    BlockDataID dst_;
    real_t omega_;
};




template<typename LM>
class StreamPullCollideD3Q19
{
public:
    StreamPullCollideD3Q19(BlockDataID src, BlockDataID dst, real_t omega)
            : src_(src), dst_(dst), omega_(omega)
    {}

    void operator()(IBlock * block)
    {
        using PdfField_T = GhostLayerField<real_t, LM::Stencil::Q>;
        auto src = block->getData<PdfField_T>( src_ );
        auto dst = block->getData<PdfField_T>( dst_ );
        using namespace stencil;
        using Stencil = typename LM::Stencil;

        real_t omega_trm_ = real_t(1) - omega_;
        real_t omega_w0_  = real_t(3) * ( real_t(1) / real_t( 3) ) * omega_;
        real_t omega_w1_  = real_t(3) * ( real_t(1) / real_t(18) ) * omega_;
        real_t omega_w2_  = real_t(3) * ( real_t(1) / real_t(36) ) * omega_;

        WALBERLA_FOR_ALL_CELLS_XYZ(src, {

            const real_t dd_tmp_NE = src->get(x-1, y-1, z  , Stencil::idx[NE]);
            const real_t dd_tmp_N  = src->get(x  , y-1, z  , Stencil::idx[N]);
            const real_t dd_tmp_NW = src->get(x+1, y-1, z  , Stencil::idx[NW]);
            const real_t dd_tmp_W  = src->get(x+1, y  , z  , Stencil::idx[W]);
            const real_t dd_tmp_SW = src->get(x+1, y+1, z  , Stencil::idx[SW]);
            const real_t dd_tmp_S  = src->get(x  , y+1, z  , Stencil::idx[S]);
            const real_t dd_tmp_SE = src->get(x-1, y+1, z  , Stencil::idx[SE]);
            const real_t dd_tmp_E  = src->get(x-1, y  , z  , Stencil::idx[E]);
            const real_t dd_tmp_T  = src->get(x  , y  , z-1, Stencil::idx[T]);
            const real_t dd_tmp_TE = src->get(x-1, y  , z-1, Stencil::idx[TE]);
            const real_t dd_tmp_TN = src->get(x  , y-1, z-1, Stencil::idx[TN]);
            const real_t dd_tmp_TW = src->get(x+1, y  , z-1, Stencil::idx[TW]);
            const real_t dd_tmp_TS = src->get(x  , y+1, z-1, Stencil::idx[TS]);
            const real_t dd_tmp_B  = src->get(x  , y  , z+1, Stencil::idx[B]);
            const real_t dd_tmp_BE = src->get(x-1, y  , z+1, Stencil::idx[BE]);
            const real_t dd_tmp_BN = src->get(x  , y-1, z+1, Stencil::idx[BN]);
            const real_t dd_tmp_BW = src->get(x+1, y  , z+1, Stencil::idx[BW]);
            const real_t dd_tmp_BS = src->get(x  , y+1, z+1, Stencil::idx[BS]);
            const real_t dd_tmp_C  = src->get(x  , y  , z  , Stencil::idx[C]);

            const real_t velX_trm = dd_tmp_E + dd_tmp_NE + dd_tmp_SE + dd_tmp_TE + dd_tmp_BE;
            const real_t velY_trm = dd_tmp_N + dd_tmp_NW + dd_tmp_TN + dd_tmp_BN;
            const real_t velZ_trm = dd_tmp_T + dd_tmp_TS + dd_tmp_TW;

            const real_t rho = dd_tmp_C + dd_tmp_S + dd_tmp_W + dd_tmp_B + dd_tmp_SW + dd_tmp_BS + dd_tmp_BW + velX_trm + velY_trm + velZ_trm;

            const real_t velX = velX_trm - dd_tmp_W  - dd_tmp_NW - dd_tmp_SW - dd_tmp_TW - dd_tmp_BW;
            const real_t velY = velY_trm + dd_tmp_NE - dd_tmp_S  - dd_tmp_SW - dd_tmp_SE - dd_tmp_TS - dd_tmp_BS;
            const real_t velZ = velZ_trm + dd_tmp_TN + dd_tmp_TE - dd_tmp_B  - dd_tmp_BN - dd_tmp_BS - dd_tmp_BW - dd_tmp_BE;

            const real_t velXX = velX * velX;
            const real_t velYY = velY * velY;
            const real_t velZZ = velZ * velZ;

            const real_t dir_indep_trm = ( real_t(1) / real_t(3) ) * rho - real_t(0.5) * ( velXX + velYY + velZZ );

            dst->get(x,y,z,Stencil::idx[C]) = omega_trm_ * dd_tmp_C + omega_w0_ * dir_indep_trm;

            const real_t vel_trm_E_W = dir_indep_trm + real_t(1.5) * velXX;
            const real_t vel_trm_N_S = dir_indep_trm + real_t(1.5) * velYY;
            const real_t vel_trm_T_B = dir_indep_trm + real_t(1.5) * velZZ;

            dst->get(x,y,z,Stencil::idx[E]) = omega_trm_ * dd_tmp_E + omega_w1_ * ( vel_trm_E_W + velX );
            dst->get(x,y,z,Stencil::idx[W]) = omega_trm_ * dd_tmp_W + omega_w1_ * ( vel_trm_E_W - velX );
            dst->get(x,y,z,Stencil::idx[N]) = omega_trm_ * dd_tmp_N + omega_w1_ * ( vel_trm_N_S + velY );
            dst->get(x,y,z,Stencil::idx[S]) = omega_trm_ * dd_tmp_S + omega_w1_ * ( vel_trm_N_S - velY );
            dst->get(x,y,z,Stencil::idx[T]) = omega_trm_ * dd_tmp_T + omega_w1_ * ( vel_trm_T_B + velZ );
            dst->get(x,y,z,Stencil::idx[B]) = omega_trm_ * dd_tmp_B + omega_w1_ * ( vel_trm_T_B - velZ );

            const real_t velXmY = velX - velY;
            const real_t vel_trm_NW_SE = dir_indep_trm + real_t(1.5) * velXmY * velXmY;

            dst->get(x,y,z,Stencil::idx[NW]) = omega_trm_ * dd_tmp_NW + omega_w2_ * ( vel_trm_NW_SE - velXmY );
            dst->get(x,y,z,Stencil::idx[SE]) = omega_trm_ * dd_tmp_SE + omega_w2_ * ( vel_trm_NW_SE + velXmY );

            const real_t velXpY = velX + velY;
            const real_t vel_trm_NE_SW = dir_indep_trm + real_t(1.5) * velXpY * velXpY;

            dst->get(x,y,z,Stencil::idx[NE]) = omega_trm_ * dd_tmp_NE + omega_w2_ * ( vel_trm_NE_SW + velXpY );
            dst->get(x,y,z,Stencil::idx[SW]) = omega_trm_ * dd_tmp_SW + omega_w2_ * ( vel_trm_NE_SW - velXpY );

            const real_t velXmZ = velX - velZ;
            const real_t vel_trm_TW_BE = dir_indep_trm + real_t(1.5) * velXmZ * velXmZ;

            dst->get(x,y,z,Stencil::idx[TW]) = omega_trm_ * dd_tmp_TW + omega_w2_ * ( vel_trm_TW_BE - velXmZ );
            dst->get(x,y,z,Stencil::idx[BE]) = omega_trm_ * dd_tmp_BE + omega_w2_ * ( vel_trm_TW_BE + velXmZ );

            const real_t velXpZ = velX + velZ;
            const real_t vel_trm_TE_BW = dir_indep_trm + real_t(1.5) * velXpZ * velXpZ;

            dst->get(x,y,z,Stencil::idx[TE]) = omega_trm_ * dd_tmp_TE + omega_w2_ * ( vel_trm_TE_BW + velXpZ );
            dst->get(x,y,z,Stencil::idx[BW]) = omega_trm_ * dd_tmp_BW + omega_w2_ * ( vel_trm_TE_BW - velXpZ );

            const real_t velYmZ = velY - velZ;
            const real_t vel_trm_TS_BN = dir_indep_trm + real_t(1.5) * velYmZ * velYmZ;

            dst->get(x,y,z,Stencil::idx[TS]) = omega_trm_ * dd_tmp_TS + omega_w2_ * ( vel_trm_TS_BN - velYmZ );
            dst->get(x,y,z,Stencil::idx[BN]) = omega_trm_ * dd_tmp_BN + omega_w2_ * ( vel_trm_TS_BN + velYmZ );

            const real_t velYpZ = velY + velZ;
            const real_t vel_trm_TN_BS = dir_indep_trm + real_t(1.5) * velYpZ * velYpZ;

            dst->get(x,y,z,Stencil::idx[TN]) = omega_trm_ * dd_tmp_TN + omega_w2_ * ( vel_trm_TN_BS + velYpZ );
            dst->get(x,y,z,Stencil::idx[BS]) = omega_trm_ * dd_tmp_BS + omega_w2_ * ( vel_trm_TN_BS - velYpZ );

        })
        src->swapDataPointers(*dst);
    }

private:
    BlockDataID src_;
    BlockDataID dst_;
    real_t omega_;
};





} // namespace walberla