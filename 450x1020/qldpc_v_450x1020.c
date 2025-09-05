#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include "runtime.h"
#include "LP_1020_450_10.h"
#include "LP_1020_450_10_logic.h"

#ifdef SPIKE
#include <stdio.h>
#elif defined ARA_LINUX
#include <stdio.h>
#else
#include "printf.h"
#endif

#define MIN 1000.0f 

uint8_t qldcp_minsum(int32_t n_iterations, float* L, float alpha, int32_t S[M], int32_t est_error[N], int32_t* it) {

    int32_t i, j, ip, jp, hard_decision, sign, parity;
    float sum, sumv[N], temp[N], min, min1, min2, R_abs, x, one = 1, minus_one = -1;
    int32_t check[M] = {0};
    

    float R[450 * 8] = {0};
    float R_aux[450 * 8] = {0};

    float o[1020 * 5] = {0};
    float o_aux[1020 * 5] = {0};

    int32_t ttt[450 * 8] = {0};

    uint8_t corrected = 0;
    
    int64_t runtime = 0;

    float maxval = 0x7FFFFFFF;

    /************************************************************************************************************** */

    while(*it < n_iterations && !corrected) {
    
        start_timer();

        //mu[i,j] = min(R[i.j']) x mul HD(R[i,j']) xor S[i]
        asm volatile("vsetvli t0, %0, e32, m1, ta, ma" :: "r"(8));
        asm volatile("vfmv.v.f v31, %0" :: "f"(maxval));
        asm volatile("vmv.v.i v30, 1"); 
        asm volatile("vmv.v.i v29, 0");
        asm volatile("vfmv.v.f v28, %0" :: "f"(one));  
        asm volatile("vfmv.v.f v27, %0" :: "f"(minus_one)); 
        asm volatile("vmv.v.i v26, -1");

        for(i = 0; i < M; i++) {

            asm volatile("vle32.v v1, (%0);" :: "r"(&R[i * 8]));

            /************************ MIN ************************/

            //Valor Absoluto
            asm volatile("vmflt.vf v0, v1, %0":: "f"(0.0));
            asm volatile("vfsgnjx.vv v2, v1, v1");

            //Hard decision
            asm volatile("vmerge.vvm v10, v29, v30, v0");            

            //Buscar minimo por fila
            asm volatile("vfredmin.vs v3, v2, v31");
            asm volatile("vfmv.f.s %0, v3" : "=f"(min1));
            asm volatile("vfmv.v.f v3, %0" :: "f"(min1));

            asm volatile("vmfeq.vf v4, v2, %0":: "f"(min1));
            asm volatile("vmsof.m v5, v4");
            asm volatile("vmnot.m v0, v5");

            asm volatile("vfredmin.vs v6, v2, v31, v0.t");
            asm volatile("vfmv.f.s %0, v6" : "=f"(min2));
            asm volatile("vmnot.m v0, v0");

            asm volatile("vfmerge.vfm v7, v3, %0, v0":: "f"(min2));

            /************************ SIGN ************************/

            //Xor de fila
            asm volatile("vredxor.vs v11, v10, v29");
            asm volatile("vmv.x.s %0, v11" : "=r"(parity));

            if(parity){
                asm volatile("vxor.vv v10, v10, v30");
            }

            //Xor con S[i]
            asm volatile("vmv.v.x v12, %0" :: "r"(S[i]));
            asm volatile("vxor.vv v10, v10, v12");

            asm volatile("vmul.vv v22, v10, v26"); 

            /********************* SIGN * MIN *********************/

            asm volatile("vfsgnj.vv v8, v7, v22");

            asm volatile("vse32.v v8, (%0);" ::"r"(&R_aux[i * 8]));

        }
        
        stop_timer();
        runtime = get_timer();
        printf("Loop 1.0: %ld ciclos\n", runtime);
        start_timer();
        
        for (i = 0; i < M; i++) {
            for (j = 0; j < 8; j++) {
                o[h_indices[i][j] * 5 + h_to_ht[i][j]] = R_aux[i * 8 + j];
            }
        }        
        
        stop_timer();
        runtime = get_timer();
        printf("Loop 1.1: %ld ciclos\n", runtime);
        
        /************** LOOP 2 **************/
        start_timer();
        
        //R[i,j] = L[j] + a * sum (mu[i'][j])
        for(j = 0; j < N; j++) {
            sum = 0;
            for(i = 0; i < 5; i++) {
                sum += o[j * 5 + i];
            }
            sumv[j] = sum;
        }
        
        stop_timer();
        runtime = get_timer();
        printf("Loop 2.0: %ld ciclos\n", runtime);
        start_timer();

        asm volatile("vsetvli t0, %0, e32, m1, ta, ma" :: "r"(100));

        asm volatile("vid.v v1");  // v1 = [0,1,2,...,99]
        asm volatile("vdivu.vx v1, v1, %0":: "r"(5));  // v1 = floor(i / 5)
        asm volatile("vsll.vi v1, v1, 2");  // v1 = v1 << 2
        
        asm volatile("vfmv.v.f v2, %0" :: "f"(alpha));

        printf("1OK\n");

        for (j = 0; j < N * 5; j+=100) {

            asm volatile("vle32.v v3, (%0);" :: "r"(&o[j]));           
            asm volatile("vluxei32.v v4, (%0), v1" :: "r"(&sumv[j/5]));
            asm volatile("vfsub.vv v5, v4, v3"); 
        
            asm volatile("vfmul.vv v6, v2, v5");  

            asm volatile("vluxei32.v v7, (%0), v1" :: "r"(&L[j/5]));
            asm volatile("vfadd.vv v8, v7, v6"); 

            asm volatile("vse32.v v8, (%0);" ::"r"(&o_aux[j]));

        }
        
        stop_timer();
        runtime = get_timer();
        printf("Loop 2.1: %ld ciclos\n", runtime);
        
        start_timer();

        for (j = 0; j < N; j++) {
            for (i = 0; i < 5; i++) {
                R[ht_indices[j][i] * 8 + ht_to_h[j][i]] = o_aux[j * 5 + i];
            }
        }

        
        stop_timer();
        runtime = get_timer();
        printf("Loop 2.2: %ld ciclos\n", runtime);
        
        /************** LOOP 3 **************/
        start_timer();

        //e[j] = HD(L[j] + a * sum (mu[i][j]))
        asm volatile("vsetvli t0, %0, e32, m1, ta, ma" :: "r"(102));
        
        asm volatile("vfmv.v.f v1, %0" :: "f"(alpha));                // v1 = alpha

        for(j = 0; j < N; j+=102) {
            
            asm volatile("vle32.v v2, (%0);" :: "r"(&sumv[j]));           // v2 = sumv[j]
            asm volatile("vfmul.vv v3, v2, v1");                          // v3 = v2 * v1

            asm volatile("vle32.v v4, (%0);" :: "r"(&L[j]));              // v4 = L[j]
            asm volatile("vfadd.vv v5, v3, v4");                          // v5 = v4 + v3
            
            asm volatile("vse32.v v5, (%0)" ::"r"(&temp[j])); //??
            
            asm volatile("vmv.v.i v31, 1");
            asm volatile("vmv.v.i v6, 0");
            asm volatile("vmflt.vf v0, v5, %0" :: "f"(0.0f));             // v0 = (v5 < 0.0)
            asm volatile("vmerge.vvm v7, v6, v31, v0");            // v6[i] = v31[i] si v0[i] == 1 else v6[i]
            
            asm volatile("vse32.v v7, (%0);" ::"r"(&est_error[j])); // est_error[j] = 1 where v5 < 0

        }
        
        stop_timer();
        runtime = get_timer();
        printf("Loop 3: %ld ciclos\n", runtime);
        
        /************** LOOP 4 **************/
        start_timer();

        // Cálculo del síndrome esperado
        corrected = 1;
        for (i = 0; i < M; i++) {
            check[i] = 0;
            for (j = 0; j < 8; j++) {   //h_lengths en h son siempre 8
                check[i] ^= est_error[h_indices[i][j]];
            }

            if (check[i] != S[i]) {
                corrected = 0;
            }
        }
        
        stop_timer();
        runtime = get_timer();
        printf("Loop 4: %ld ciclos\n", runtime);

        (*it)++;
    }

    return corrected;
}

uint8_t logic_permutation(int32_t error[N], int32_t logical_error[N_logic]) {

    int32_t i, j;

    for (i = 0; i < N_logic; ++i) {
        logical_error[i] = 0;
        for (j = 0; j < ht_logic_lengths[i]; ++j) {
            logical_error[i] ^= error[ht_logic_indices[i][j]];
        }
    }

    return 0; 
}

int main() {

    //Genreación de semilla    
    srand(0);
    //printf("SEED: %d\n", seed);
    //srand(time(NULL));

    float alpha = 0.75, p_error = 0.02f;
    int32_t error[N] = {0}, est_error[N] = {0}, S[M] = {0}, logical_error[N_logic], error_X[N] = {0}, error_Z[N] = {0};
    float noise[N];
    int32_t i, j, l, n_iterations = 50, error_rate = 10000, it;
    uint8_t corrected = 0;
    int64_t runtime = 0;

    // Inicializar todos los LLRs a 1.0 
    float L[N];
    for (int i = 0; i < N; ++i) {
        L[i] = 1.0f;
    }
    
    for(l=0; l < error_rate; l++) {
    
      // ----------------------------
      // Generación de errores aleatorios (Pauli X, Z, Y)
      // ----------------------------
      for (i = 0; i < N; ++i) {
          
          est_error[i] = 0;
          error_X[i] = 0;
          error_Z[i] = 0;
          noise[i] = (float)rand() / RAND_MAX;

          if (noise[i] < p_error / 3.0f) {
              error_X[i] = 1;
          } else if (noise[i] < 2.0f * p_error / 3.0f) {
              error_Z[i] = 1;
          } else if (noise[i] < p_error) {
              error_X[i] = 1;
              error_Z[i] = 1;
          }
      }

      // Seleccionar qué tipo de error pasar al decodificador (Z)
      for (i = 0; i < N; ++i) {
          error[i] = error_Z[i];
          //error[i] = error_X[i];
          //error[i] = error_X[i] ^ error_Z[i];
      }

      //Calcular el Sindrome generado por el error
      for (i = 0; i < h_M; i++) {
          S[i] = 0;
          for (j = 0; j < h_lengths[i]; j++) {
              S[i] ^= error[h_indices[i][j]];
          }
      }
      
      it = 0;

      // Ejecutar decodificador
      //start_timer();
      corrected = qldcp_minsum(n_iterations, L, alpha, S, est_error, &it);
      //stop_timer();
      
      //runtime = get_timer();
      
      if(corrected){
        printf("%d, %d, %ld\n", l, it, runtime);
      }else{
        printf("%d, %d, 0\n", l, it);
      }
      
    }
    
    return 0;
}