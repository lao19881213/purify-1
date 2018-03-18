#ifndef PURIFY_PROJECTION_KERNELS_H
#define PURIFY_PROJECTION_KERNELS_H

#include "purify/config.h"
#include "purify/logging.h"
#include "purify/types.h"

namespace purify {
namespace projection_kernels {
//! parabloic approximation for curvature
inline t_complex w_proj_approx(const t_real &u, const t_real &v, const t_real &w, const t_real &dl,
                               const t_real &dm);
//! returns lambda for w projection kernel
std::function<t_complex(t_real, t_real, t_real)> const
w_projection_kernel(const t_real &cellx, const t_real &celly, const t_uint &imsizex,
                    const t_uint &imsizey);
} // namespace projection_kernels
} // namespace purify

#endif
