// aether_sims: the identity of this engine build (the fork tag its source
// carries). Checked by the harness before any row is run or a cache is hit.
#ifndef JOSIM_ENGINEID_HPP
#define JOSIM_ENGINEID_HPP

namespace JoSIM {
// v2 (T35 W0.2): v1 (the Johnson source fix) + the run-level coupling time
// constant .ctrltau; without .ctrltau every run is bit-identical to v1.
constexpr const char* ENGINE_ID = "aether-josim-v2-ctrl-tau";
}  // namespace JoSIM

#endif  // JOSIM_ENGINEID_HPP
