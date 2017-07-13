#ifndef FURRY_TINY_MVS_OPT_H
#define FURRY_TINY_MVS_OPT_H


#include "depth_map.h"
#include "cost_volume.h"

#include "mrf_settings.pb.h"

namespace furry {
namespace tiny {

void WtaDepth(const CostVolume &cost_volume, DepthMap *depthmap);

void SolveContinuousDepth(const ContinuousMrfSettings &settings,
                          const CostVolume &cost_volume,
                          DepthMap *depthmap);

void SolveDiscreteDepth(const DiscreteMrfSettings &settings,
                        const CostVolume &cost_volume,
                        DepthMap *depthmap);

void SolveByDenseCrf(const DenseCrfSettings &settings,
                     const CostVolume &cost_volume,
                     DepthMap *depthmap);

} // tiny
} // furry

#endif // FURRY_TINY_MVS_OPT_H
