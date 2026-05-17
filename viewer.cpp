// SDL2 live-view backend for DRX_CA. Compiled as its own translation unit so
// the project's helpers.h macros never reach SDL's headers.

#include "viewer.h"
#include <SDL2/SDL.h>
#include <cstdint>
#include <cstdio>

namespace viewer {

static SDL_Window*   g_window   = nullptr;
static SDL_Renderer* g_renderer = nullptr;
static SDL_Texture*  g_texture  = nullptr;

static int  g_grid_size           = 0;
static int  g_scale               = 4;
static bool g_quit_requested      = false;
static bool g_screenshot_pending  = false;
static int  g_screenshot_counter  = 0;

// matplotlib 'tab20' palette in ARGB8888.
static const uint32_t k_palette[20] = {
    0xff1f77b4u, 0xffaec7e8u, 0xffff7f0eu, 0xffffbb78u,
    0xff2ca02cu, 0xff98df8au, 0xffd62728u, 0xffff9896u,
    0xff9467bdu, 0xffc5b0d5u, 0xff8c564bu, 0xffc49c94u,
    0xffe377c2u, 0xfff7b6d2u, 0xff7f7f7fu, 0xffc7c7c7u,
    0xffbcbd22u, 0xffdbdb8du, 0xff17becfu, 0xff9edae5u,
};

static inline uint32_t color_for(int grain_id) {
    int idx = ((grain_id - 1) % 20 + 20) % 20;
    return k_palette[idx];
}

bool init(int grid_size) {
    g_grid_size = grid_size;
    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        std::fprintf(stderr, "viewer: SDL_Init failed: %s\n", SDL_GetError());
        return false;
    }
    // Pick a scale that gives roughly an 800px window without going huge.
    g_scale = (grid_size > 0 && grid_size < 800) ? (800 / grid_size) : 1;
    if (g_scale < 1) g_scale = 1;

    g_window = SDL_CreateWindow(
        "DRX_CA — q/Esc quit, s screenshot",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
        grid_size * g_scale, grid_size * g_scale,
        SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE
    );
    if (!g_window) {
        std::fprintf(stderr, "viewer: SDL_CreateWindow failed: %s\n", SDL_GetError());
        SDL_Quit();
        return false;
    }
    g_renderer = SDL_CreateRenderer(g_window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if (!g_renderer) {
        g_renderer = SDL_CreateRenderer(g_window, -1, SDL_RENDERER_SOFTWARE);
    }
    if (!g_renderer) {
        std::fprintf(stderr, "viewer: SDL_CreateRenderer failed: %s\n", SDL_GetError());
        SDL_DestroyWindow(g_window);
        SDL_Quit();
        return false;
    }
    g_texture = SDL_CreateTexture(
        g_renderer, SDL_PIXELFORMAT_ARGB8888,
        SDL_TEXTUREACCESS_STREAMING, grid_size, grid_size
    );
    if (!g_texture) {
        std::fprintf(stderr, "viewer: SDL_CreateTexture failed: %s\n", SDL_GetError());
        SDL_DestroyRenderer(g_renderer);
        SDL_DestroyWindow(g_window);
        SDL_Quit();
        return false;
    }
    SDL_SetTextureScaleMode(g_texture, SDL_ScaleModeNearest);
    return true;
}

static void pump_events() {
    SDL_Event ev;
    while (SDL_PollEvent(&ev)) {
        switch (ev.type) {
            case SDL_QUIT:
                g_quit_requested = true;
                break;
            case SDL_KEYDOWN:
                if (ev.key.keysym.sym == SDLK_q || ev.key.keysym.sym == SDLK_ESCAPE) {
                    g_quit_requested = true;
                } else if (ev.key.keysym.sym == SDLK_s) {
                    g_screenshot_pending = true;
                }
                break;
            case SDL_WINDOWEVENT:
                if (ev.window.event == SDL_WINDOWEVENT_CLOSE) {
                    g_quit_requested = true;
                }
                break;
            default:
                break;
        }
    }
}

static void save_screenshot() {
    int w = 0, h = 0;
    SDL_GetRendererOutputSize(g_renderer, &w, &h);
    SDL_Surface* surface = SDL_CreateRGBSurfaceWithFormat(0, w, h, 32, SDL_PIXELFORMAT_ARGB8888);
    if (!surface) {
        std::fprintf(stderr, "viewer: screenshot surface alloc failed: %s\n", SDL_GetError());
        return;
    }
    if (SDL_RenderReadPixels(g_renderer, nullptr, SDL_PIXELFORMAT_ARGB8888, surface->pixels, surface->pitch) != 0) {
        std::fprintf(stderr, "viewer: SDL_RenderReadPixels failed: %s\n", SDL_GetError());
        SDL_FreeSurface(surface);
        return;
    }
    char fname[64];
    std::snprintf(fname, sizeof(fname), "drx_frame_%05d.bmp", g_screenshot_counter++);
    if (SDL_SaveBMP(surface, fname) != 0) {
        std::fprintf(stderr, "viewer: SDL_SaveBMP failed: %s\n", SDL_GetError());
    } else {
        std::printf("viewer: saved %s\n", fname);
    }
    SDL_FreeSurface(surface);
}

bool render(const int* grain_ids, int grid_size, int iteration, double eps, double p_avg) {
    pump_events();
    if (g_quit_requested) return false;

    uint32_t* pixels = nullptr;
    int pitch = 0;
    if (SDL_LockTexture(g_texture, nullptr, reinterpret_cast<void**>(&pixels), &pitch) != 0) {
        std::fprintf(stderr, "viewer: SDL_LockTexture failed: %s\n", SDL_GetError());
        return true;
    }
    for (int y = 0; y < grid_size; ++y) {
        uint32_t* row = reinterpret_cast<uint32_t*>(reinterpret_cast<uint8_t*>(pixels) + y * pitch);
        // Flip y so the image matches the previous gnuplot orientation (origin bottom-left).
        int src_y = grid_size - 1 - y;
        const int* src_row = grain_ids + src_y * grid_size;
        for (int x = 0; x < grid_size; ++x) {
            row[x] = color_for(src_row[x]);
        }
    }
    SDL_UnlockTexture(g_texture);

    SDL_SetRenderDrawColor(g_renderer, 0, 0, 0, 255);
    SDL_RenderClear(g_renderer);
    SDL_RenderCopy(g_renderer, g_texture, nullptr, nullptr);
    SDL_RenderPresent(g_renderer);

    if (g_screenshot_pending) {
        save_screenshot();
        g_screenshot_pending = false;
    }

    char title[160];
    std::snprintf(title, sizeof(title),
        "DRX_CA — iter %d  eps %.3f  p_avg %.3e   (q/Esc quit, s screenshot)",
        iteration, eps, p_avg);
    SDL_SetWindowTitle(g_window, title);
    return true;
}

void shutdown() {
    if (g_texture)  { SDL_DestroyTexture(g_texture);   g_texture  = nullptr; }
    if (g_renderer) { SDL_DestroyRenderer(g_renderer); g_renderer = nullptr; }
    if (g_window)   { SDL_DestroyWindow(g_window);     g_window   = nullptr; }
    SDL_Quit();
}

} // namespace viewer
