#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include "runtime.h"
#include "Hmatrixn144.h"

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
    float sum, sumv[N], temp[N], min, R_abs, x;
    int32_t check[M] = {0};

    float R[M * h_MAX_NNZ] = {0};
    float R_aux[M * h_MAX_NNZ] = {0};

    float o[N * ht_MAX_NNZ] = {0};
    float o_aux[N * ht_MAX_NNZ] = {0};

    uint8_t corrected = 0;
    
    int64_t runtime = 0;

    /************************************************************************************************************** */

    while(*it < n_iterations && !corrected) {
      
        //mu[i,j] = min(R[i.j']) x mul HD(R[i,j']) xor S[i]
        for(i = 0; i < M; i++) {
            start_timer();
            for(j = 0; j < h_lengths[i]; j++) {     //h_lengths en h son siempre 8
                
                min = MIN;
                parity = S[i];

                //En el mismo bucle buscamos el minimo y el signo
                for (jp = 0; jp < h_lengths[i]; jp++) {    //h_lengths en h son siempre 8
                    if(j != jp) {
                        R_abs = fabs(R[i * h_MAX_NNZ + jp]);
                        if(R_abs < min){
                            min = R_abs;
                        }
                         
                        //Paridad, 1 si es negativo y 0 si es positivo
                        hard_decision = (R[i * h_MAX_NNZ + jp] < 0) ? 1 : 0;
                        parity ^= hard_decision;
                    }
                }
                //Le añadimos el signo segun la paridad
                sign = (parity == 0) ? 1 : -1;
                R_aux[i * h_MAX_NNZ + j] = min * sign;
                               
            }
            stop_timer();
            runtime = get_timer();
            printf("%i: %ld ciclos\n", i, runtime);
        }
        
        start_timer();        
        for (i = 0; i < M; i++) {
            for (j = 0; j < h_lengths[i]; j++) {
                o[h_indices[i][j] * ht_MAX_NNZ + h_to_ht[i][j]] = R_aux[i * h_MAX_NNZ + j];
            }
        }       
        
        stop_timer();
        runtime = get_timer();
        printf("Loop 1.1: %ld ciclos\n", runtime);
        
        start_timer(); 
        
        for(j = 0; j < N; j++) {
            sum = 0;
            for(i = 0; i < 6; i++) {
                sum += o[j * ht_MAX_NNZ + i];
            }
            sumv[j] = sum;
        }
        
        stop_timer();
        runtime = get_timer();
        printf("Loop 2.0: %ld ciclos\n", runtime);
        start_timer();

        //R[i,j] = L[j] + a * sum (mu[i'][j])
        for(j = 0; j < N; j++) {
            for(i = 0; i < ht_lengths[j]; i++) {
        
                sum = sumv[j] - o[j * ht_MAX_NNZ + i];
                o_aux[j * ht_MAX_NNZ + i] = L[j] + alpha * sum;
                
            }
        }
        
        stop_timer();
        runtime = get_timer();
        printf("Loop 2.1: %ld ciclos\n", runtime);
        start_timer();
        
        for (j = 0; j < N; j++) {
            for (i = 0; i < ht_lengths[j]; i++) {
                R[ht_indices[j][i] * h_MAX_NNZ + ht_to_h[j][i]] = o_aux[j * ht_MAX_NNZ + i];
            }
        }
        
        stop_timer();
        runtime = get_timer();
        printf("Loop 2.2: %ld ciclos\n", runtime);
        
        start_timer();

        //e[j] = HD(L[j] + a * sum (mu[i][j]))
        for(j = 0; j < N; j++) {

            x = L[j] + alpha * sumv[j];
            est_error[j] = (x < 0) ? 1 : 0;

        }
        
        stop_timer();
        runtime = get_timer();
        printf("Loop 3.0: %ld ciclos\n", runtime);

        start_timer();

        // Cálculo del síndrome esperado
        corrected = 1;
        for (i = 0; i < M; i++) {
            check[i] = 0;
            for (j = 0; j < h_lengths[i]; j++) {  
                check[i] ^= est_error[h_indices[i][j]];
            }

            if (check[i] != S[i]) {
                corrected = 0;
            }
        }
        
        stop_timer();
        runtime = get_timer();
        printf("Loop 4.0: %ld ciclos\n", runtime);

        (*it)++;
    }

    return corrected;
}


int main() {

    //Genreación de semilla    
    srand(0);
    //printf("SEED: %d\n", seed);
    //srand(time(NULL));

    float alpha = 0.75, p_error = 0.001f;
    int32_t error[N] = {0}, est_error[N] = {0}, S[M] = {0}, error_X[N] = {0}, error_Z[N] = {0};
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