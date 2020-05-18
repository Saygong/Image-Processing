#include "ip_lib.h"
#include "bmp.h"
#include <stdlib.h>
#include <stdio.h>

int main()
{   

    Bitmap *bm;
    
    ip_mat *mat=ip_mat_create(100,100,3,0.0);

    ip_mat_init_random(mat, 150.0, 75.0);

    ip_mat_show(mat);

    bm = ip_mat_to_bitmap(mat);
    bm_save(bm, "random.bmp");
    ip_mat_free(mat);
    


/*

    Bitmap *bm = bm_load("fullmoon.bmp");


    ip_mat *mat = bitmap_to_ip_mat(bm);
    
    mat=ip_mat_corrupt(mat, 100);

    bm = ip_mat_to_bitmap(mat);
    bm_save(bm, "corone.bmp");

    bm_free(bm);
    
    ip_mat_free(mat);

    
    
    

    ip_mat * matrix1, *risultato;
    
    
    Bitmap *bm1 = bm_load("mongolfiere.bmp");
    Bitmap *bmbri;
    
        
    
    ip_mat *mat1 = bitmap_to_ip_mat(bm1);
    
    ip_mat *bright = ip_mat_corrupt(mat, 50);
    
  
    
    bmbri = ip_mat_to_bitmap(bright);
    
    
    
    ip_mat_free(matrix1);
    ip_mat_free(risultato);
    */
    printf("AHOI questa è la fine\n");
    return 0;
}
