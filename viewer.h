#ifndef _DRX_VIEWER_H_
#define _DRX_VIEWER_H_

// Standalone SDL2-based grid viewer. Intentionally does NOT include any of the
// project's other headers, because helpers.h defines short single-letter
// macros (b, alpha, T, ...) that collide with SDL's symbol names. Keep this
// translation unit isolated.

namespace viewer {
    // Open the window. grid_size is the number of cells along one edge.
    // Returns false on failure.
    bool init(int grid_size);

    // Render one frame from a row-major flat array of grain IDs.
    // Returns false if the user requested to close the window.
    bool render(const int* grain_ids, int grid_size, int iteration, double eps, double p_avg);

    // Tear down the window and SDL.
    void shutdown();
}

#endif
