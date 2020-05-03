#ifndef _DRX_HELPERS_CPP_
#define _DRX_HELPERS_CPP_

#include "helpers.h"

// inline bool cell_on_border(int i, int j) {
//     int x = 0, y = 0;
//     int x_start = 0, x_end = 3, y_start = 0, y_end = 3;
//     if (i == 0) x_start++;
//     if (j == 0) y_start++;
//     if (i == GRID_SIZE - 1) x_end--;
//     if (j == GRID_SIZE - 1) y_end--;

//     for (x = x_start; x < x_end; x++) {
//         for (y = y_start; y < y_end; y++) {
//             if (x == y) continue;

//             if (grid.cell[i][j].grain_number != grid.cell[i + x - 1][j + y - 1].grain_number) {
//                 return 1;
//             }
//         }
//     }

//     // ff(x, x_start, x_end) { 
//     //     ff(y, y_start, y_end) {
//     //         if (x == y) continue;

//     //         if (grid[i][j].grain_number != grid[i + x - 1][j + y - 1].grain_number) {
//     //             return 1;
//     //         }
//     //     }
//     // }
//     return 0;
// }

inline float gamma_l(float misorientation) {
    if (misorientation == 0) {
        return 0;
    }
    else if (misorientation >= CRITICAL_MISORIENTATION) {
        return gamma_o;
    } else {
        return gamma_o * (misorientation / CRITICAL_MISORIENTATION) * (1 - log(misorientation / CRITICAL_MISORIENTATION));
    }
}

inline float mobility(float misorientation) {
    if (misorientation == 0) { 
        return 0;
    } else {
        return Mo * (1 - exp(-5 * pow((misorientation / CRITICAL_MISORIENTATION), 4)));
    }
}

void print_array(int array[GRID_SIZE][GRID_SIZE]) {
    int i = 0, j = 0;
    ff(i, 0, GRID_SIZE) {
        ff(j, 0, GRID_SIZE) {
            cout << array[i][j] << " ";
        }
        cout << "\n";
    }
}

#endif