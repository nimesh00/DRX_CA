#include <iostream>
#include <random>
#include <chrono>
#include "helpers.cpp"
#include "drx_grid.cpp"

int main() {
    // Invoking GnuPlot in the background.
    system("gnuplot cadrx.gp &");

    // seeding random number initialization with current time
    srand((unsigned)time(0));

    _grid_ grid;
    _grid_ u_grid;
    deep_copy_grid(grid, u_grid);
    // print_grid(grid.cell);
    // cout << "\n\n";
    // print_array(grid.grain_num);

    // timing variables
    auto start = chrono::steady_clock::now();
    auto end = chrono::steady_clock::now();;

    bool on_border = 0;
    int x = 0, y = 0;
    int nx = 0, ny = 0;
    int x_start = 0, x_end = 3, y_start = 0, y_end = 3;

    float cell_gamma = 0;
    float delta_p = 0;
    float misorientation = 0;
    float grain_size = 0;
    
    float eps = 0;
    float delta_t = 0;
    float delta_eps = 0;
    float delta_n = 0;

    vector<int> potential_nucleus_x {-1};
    vector<int> potential_nucleus_y {-1};
    int pot_n[GRID_SIZE][GRID_SIZE];
    float p_nucleation = 0;
    float p_growth = 0;
    float p_random = 0;
    int pn_i = 0;
    int nucleus_counter = 0;
    int neighbour_info = 0;

    int iteration = 0;

    int i = 0, j = 0;

    write_to_file(grid.cell);

    while (eps < EPS_FINAL) {
        cout << "Iteration: " << ++iteration << endl;
        delta_t = CELL_SIZE / grid.v_max;
        delta_eps = EPS_RATE * delta_t;
        eps += delta_eps;
        if (eps > EPS_CR) {
            pn_i = 0;
            nucleus_counter = 0;

            // Calculate number of new nucleus to be formed
            delta_n = (int)(nucleation_rate * delta_t);
            // cout << "total nucleus: " << delta_n << endl;

            // Probablistic approach to find possible nucleus
            start = chrono::steady_clock::now();

            while (nucleus_counter != delta_n + 1) {
                i = rand() % GRID_SIZE;
                j = rand() % GRID_SIZE;
                if (pot_n[i][j] == 1) {
                    continue;
                }
                neighbour_info = 0;
                neighbour_info = grid.check_neighbours(i, j);
                if (neighbour_info != 0) {
                    nx = (int)(neighbour_info / pow(10, encoder + 1));
                    neighbour_info = neighbour_info % (int)pow(10, encoder + 1);
                    ny = (int)neighbour_info / 10;
                    misorientation = abs(grid.cell[i][j].orientation - grid.cell[nx][ny].orientation) * M_PI / 180;
                    cell_gamma = gamma_m * sin(2 * misorientation) * (1 - r_gamma * log(sin(2 * misorientation)));
                    p_ciritcal = pow((20 * cell_gamma * EPS_RATE) / (3 * b * sqrt(grid.cell[i][j].dislocation_density) * M * tau * tau), 1.0 / 3.0);
                    if (grid.cell[i][j].dislocation_density > p_ciritcal) {
                        potential_nucleus_x.push_back(i);
                        potential_nucleus_y.push_back(j);
                        pot_n[i][j] = 1;
                        nucleus_counter++;
                    }
                }
            }
            end = chrono::steady_clock::now();
            // cout << "probabilistic approach time: " << chrono::duration_cast<chrono::microseconds>(end - start).count() << "us" << endl;


            // // Identify possible nucleation sites (Checking all possible cells) (Deterministic approach) (VERY SLOW!!)

            // start = chrono::steady_clock::now();
            // ff(i, 0, GRID_SIZE) {
            //     if (nucleus_counter == delta_n + 1) break;
            //     ff(j, 0, GRID_SIZE) {
            //         neighbour_info = 0;
            //         neighbour_info = grid.check_neighbours(i, j);
            //         if (neighbour_info != 0) {
            //             nx = (int)(neighbour_info / pow(10, encoder + 1));
            //             neighbour_info = neighbour_info % (int)pow(10, encoder + 1);
            //             ny = (int)neighbour_info / 10;
            //             misorientation = abs(grid.cell[i][j].orientation - grid.cell[nx][ny].orientation) * M_PI / 180;
            //             cell_gamma = gamma_m * sin(2 * misorientation) * (1 - r_gamma * log(sin(2 * misorientation)));
            //             p_ciritcal = pow((20 * cell_gamma * EPS_RATE) / (3 * b * sqrt(grid.cell[i][j].dislocation_density) * M * tau * tau), 1.0 / 3.0);
            //             if (grid.cell[i][j].dislocation_density > p_ciritcal) {
            //                 potential_nucleus_x.push_back(i);
            //                 potential_nucleus_y.push_back(j);
            //                 nucleus_counter++;
            //             }
            //         }
            //     }
            // }
            // end = chrono::steady_clock::now();
            // cout << "Deterministic approach time: " << chrono::duration_cast<chrono::microseconds>(end - start).count() << "us" << endl;

            // perform nucleation
            // nucleus_counter--;
            // cout << nucleus_counter << endl;
            // cout << "vector size: " << potential_nucleus_y[0] << endl;

            ff(pn_i, 0, nucleus_counter) {
                pn_i = rand() % potential_nucleus_x.size();
                if (potential_nucleus_y[pn_i] == -1) {
                    continue;
                }
                if (u_grid.cell[potential_nucleus_x[pn_i]][potential_nucleus_y[pn_i]].N_recrystallized == 1) {
                    continue;
                }
                p_random = rand() / RAND_MAX;
                p_nucleation = grid.cell[potential_nucleus_x[pn_i]][potential_nucleus_y[pn_i]].dislocation_density / grid.p_max;
                if (p_random < p_nucleation && grid.cell[potential_nucleus_x[pn_i]][potential_nucleus_y[pn_i]].N_recrystallized != 1) {
                    u_grid.cell[potential_nucleus_x[pn_i]][potential_nucleus_y[pn_i]].nucleate();
                    nucleus_counter--;
                    potential_nucleus_x[pn_i] = potential_nucleus_y[pn_i] = -1;
                }
            }

            // while (nucleus_counter != 0) {
            //     // cout << "Nucleus Counter: " << nucleus_counter << endl;
            //     pn_i = rand() % potential_nucleus_x.size();
            //     if (potential_nucleus_y[pn_i] == -1) {
            //         continue;
            //     }
            //     if (u_grid.cell[potential_nucleus_x[pn_i]][potential_nucleus_y[pn_i]].N_recrystallized == 1) {
            //         continue;
            //     }
            //     p_random = rand() / RAND_MAX;
            //     p_nucleation = grid.cell[potential_nucleus_x[pn_i]][potential_nucleus_y[pn_i]].dislocation_density / grid.p_max;
            //     if (p_random < p_nucleation && grid.cell[potential_nucleus_x[pn_i]][potential_nucleus_y[pn_i]].N_recrystallized != 1) {
            //         u_grid.cell[potential_nucleus_x[pn_i]][potential_nucleus_y[pn_i]].nucleate();
            //         nucleus_counter--;
            //         potential_nucleus_x[pn_i] = potential_nucleus_y[pn_i] = -1;
            //     }
            // }

            // cout << "Hii\n";
            
            ff (i, 0, GRID_SIZE) {
                ff(j, 0, GRID_SIZE) {
                    // re-initializing here to save doing iterations later
                    pot_n[i][j] = 0;

                    // cell routine
                    on_border = 0;
                    x = y = 0;
                    nx = ny = 0;
                    x_start = y_start = 0;
                    x_end = y_end = 3;
                    if (i == 0) x_start++;
                    if (j == 0) y_start++;
                    if (i == GRID_SIZE - 1) x_end--;
                    if (j == GRID_SIZE - 1) y_end--;

                    ff (x, x_start, x_end) {
                        ff (y, y_start, y_end) {
                            nx = i + x - 1;
                            ny = j + y - 1;
                            if (i == nx && j == ny) continue;
                            if (grid.cell[i][j].grain_number != grid.cell[nx][ny].grain_number && u_grid.cell[nx][ny].N_recrystallized == 1) {
                                on_border = 1;
                                break;
                            }
                        }
                        if (on_border == 1) {
                            break;
                        }
                    }

                    if (on_border == 1) {
                        misorientation = abs(grid.cell[i][j].orientation - grid.cell[nx][ny].orientation);
                        delta_p = abs(grid.cell[i][j].dislocation_density - grid.cell[nx][ny].dislocation_density);
                        grain_size = sqrt(grid.grain_size_nc[grid.grain_num[i][j]] * CELL_SIZE * CELL_SIZE / M_PI);
                        u_grid.cell[i][j].gb_velocity = grid.calculate_cell_velocity(delta_p, misorientation, grain_size);
                        if (u_grid.cell[i][j].gb_velocity > u_grid.v_max) {
                            u_grid.v_max = grid.cell[i][j].gb_velocity;
                        }
                        p_growth = u_grid.cell[i][j].gb_velocity / grid.v_max;
                        p_random = rand() / RAND_MAX;
                        if (p_random < p_growth) {
                            u_grid.cell[i][j].grain_number = u_grid.cell[nx][ny].grain_number;
                            u_grid.cell[i][j].dislocation_density = u_grid.cell[nx][ny].dislocation_density;
                            u_grid.cell[i][j].N_recrystallized = 1;
                            u_grid.cell[i][j].orientation = u_grid.cell[nx][ny].orientation;
                        } else {
                            grid.cell[i][j].update_dislocation_density(delta_eps);
                        }
                    } else {
                        grid.cell[i][j].update_dislocation_density(delta_eps);
                    }
                }
            }
        }
        potential_nucleus_x.clear();
        potential_nucleus_y.clear();
        deep_copy_grid(u_grid, grid);
        // write_to_file(grid.cell);
    }
    write_to_file(u_grid.cell);

    return 0;
}