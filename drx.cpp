#include <iostream>
#include <random>
#include "helpers.cpp"
#include "drx_grid.cpp"

int main() {
    // Invoking GnuPlot in the background.
    system("gnuplot cadrx.gp &");

    // seeding random number initialization with current time
    srand((unsigned)time(0));

    _grid_ grid;
    _grid_ u_grid;
    // print_grid(grid.cell);
    // cout << "\n\n";
    // print_array(grid.grain_num);

    bool on_border = 0;
    int x = 0, y = 0;
    int nx = 0, ny = 0;
    int x_start = 0, x_end = 3, y_start = 0, y_end = 3;

    float cell_gamma = 0;
    float misorientation = 0;
    
    float eps = 0;
    float delta_t = 0;
    float delta_eps = 0;
    float delta_n = 0;

    vector<int> potential_nucleus_x {-1};
    vector<int> potential_nucleus_y {-1};
    float p_nucleation = 0;
    float p_random = 0;
    int pn_i = 0;
    int nucleus_counter = 0;

    // while (eps < EPS_FINAL) {
    //     delta_t = CELL_SIZE / grid.v_max;
    //     delta_eps = EPS_RATE * delta_t;
    //     if (eps > EPS_CR) {
    //         pn_i = 0;
    //         nucleus_counter = 0;
    //         delta_n = (int)(nucleation_rate * delta_t);
    //         ff(i, 0, GRID_SIZE) {
    //             if (nucleus_counter == delta_n + 1) break;
    //             ff(j, 0, GRID_SIZE) {
    //                 on_border = 0;
    //                 x = y = 0;
    //                 nx = ny = 0;
    //                 x_start = y_start = 0;
    //                 x_end = y_end = 3;
    //                 if (i == 0) x_start++;
    //                 if (j == 0) y_start++;
    //                 if (i == GRID_SIZE - 1) x_end--;
    //                 if (j == GRID_SIZE - 1) y_end--;
    //                 for (x = x_start; x < x_end; x++) {
    //                     for (y = y_start; y < y_end; y++) {
    //                         if (i == x && j == y) continue;
    //                         nx = i + x - 1;
    //                         ny = j + y - 1;
    //                         if (grid.cell[i][j].grain_number != grid.cell[nx][ny].grain_number) {
    //                             on_border = 1;
    //                             break;
    //                         }
    //                     }
    //                     if (on_border == 1) {
    //                         break;
    //                     }
    //                 }
    //                 if (on_border) {
    //                     misorientation = abs(grid.cell[i][j].orientation - grid.cell[nx][ny].orientation) * M_PI / 180;
    //                     cell_gamma = gamma_m * sin(2 * misorientation) * (1 - r_gamma * log(sin(2 * misorientation)));
    //                     p_ciritcal = pow((20 * cell_gamma * EPS_RATE) / (3 * b * sqrt(grid.cell[i][j].dislocation_density) * M * tau * tau), 1.0 / 3.0);
    //                     if (grid.cell[i][j].dislocation_density > p_ciritcal) {
    //                         potential_nucleus_x.push_back(i);
    //                         potential_nucleus_y.push_back(j);
    //                         nucleus_counter++;
    //                     }
    //                 }
    //             }
    //         }
    //         while (nucleus_counter != 0) {
    //             pn_i = rand() % potential_nucleus_x.size();
    //             if (potential_nucleus_y[pn_i] == -1) {
    //                 continue;
    //             }
    //             p_random = rand() / RAND_MAX;
    //             p_nucleation = grid.cell[potential_nucleus_x[pn_i]][potential_nucleus_y[pn_i]].dislocation_density / grid.p_max;
    //             if (p_random < p_nucleation) {
    //                 grid.cell[potential_nucleus_x[pn_i]][potential_nucleus_y[pn_i]].nucleate();
    //                 nucleus_counter--;
    //                 potential_nucleus_x[pn_i] = potential_nucleus_y[pn_i] = -1;
    //             }
    //         }
    //     }


    //     potential_nucleus_x.clear();
    //     potential_nucleus_y.clear();
    // }
    write_to_file(grid.cell);

    return 0;
}