#include "ip_lib.h"
#include "bmp.h"
#include <stdlib.h>
#include <stdio.h>

int main()
{   

    Bitmap *bm = bm_load("flower.bmp");
    /*Bitmap *bm1 = bm_load("mongolfiere.bmp");*/

    ip_mat *mat = bitmap_to_ip_mat(bm);
    /*ip_mat *mat1 = bitmap_to_ip_mat(bm1);*/
    ip_mat *filt;
    ip_mat *res;
    ip_mat *res1;
     
    filt = create_average_filter(9,9,1);
    /*res1 = ip_mat_corrupt(mat, 0);*/
    res = ip_mat_convolve(mat, filt);
    printf("sono uscito\n");

    /*res = ip_mat_blend(mat, mat1, 0.5);*/

    bm = ip_mat_to_bitmap(res);
    bm_save(bm, "oi_gesubambino.bmp");

    /*bm = ip_mat_to_bitmap(res1);
    bm_save(bm, "oi_giacomo.bmp");*/

    bm_free(bm);
    
    ip_mat_free(mat);
    ip_mat_free(filt);
    ip_mat_free(res);
    /*ip_mat_free(res1);*/

  /*
    Bitmap *bm;
    
    ip_mat *mat=ip_mat_create(100,100,3,0.0);

    ip_mat_init_random(mat, 150.0, 75.0);

    ip_mat_show(mat);

    bm = ip_mat_to_bitmap(mat);
    bm_save(bm, "random.bmp");
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
