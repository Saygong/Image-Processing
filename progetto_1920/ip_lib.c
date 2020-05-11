/*
 Created by Sebastiano Vascon on 23/03/20.
*/

#include <stdio.h>
#include <math.h>
#include "ip_lib.h"
#include "bmp.h"

float clamp_f(float f, float lo, float hi)
{
    if (f > hi)
    {
        return 255.0;
    }
    else if (f < lo)
    {
        return lo;
    }
    else
    {
        return f;
    }
}

void ip_mat_show(ip_mat * t){
    unsigned int i,l,j;
    printf("Matrix of size %d x %d x %d (hxwxk)\n",t->h,t->w,t->k);
    for (l = 0; l < t->k; l++) {
        printf("Slice %d\n", l);
        for(i=0;i<t->h;i++) {
            for (j = 0; j < t->w; j++) {
                printf("%f ", get_val(t,i,j,l));
            }
            printf("\n");
        }
        printf("\n");
    }
}

void ip_mat_show_stats(ip_mat * t){
    unsigned int k;

    compute_stats(t);

    for(k=0;k<t->k;k++){
        printf("Channel %d:\n", k);
        printf("\t Min: %f\n", t->stat[k].min);
        printf("\t Max: %f\n", t->stat[k].max);
        printf("\t Mean: %f\n", t->stat[k].mean);
    }
}


void print_matrix(ip_mat *t){
    int i,l,j;
    for (l = 0; l < t->k; l++) {
        printf("Slice %d\n", l);
        for(i=0;i< t->h ;i++) {
            for (j = 0; j < t->w; j++) {
                printf("%f ", t->data[l][i][j]);
            }
            printf("\n");
        }
        printf("\n");
    }
}


ip_mat * bitmap_to_ip_mat(Bitmap * img){
    unsigned int i=0,j=0;

    unsigned char R,G,B;

    unsigned int h = img->h;
    unsigned int w = img->w;

    ip_mat * out = ip_mat_create(h, w,3,0);

    for (i = 0; i < h; i++)              /* rows */
    {
        for (j = 0; j < w; j++)          /* columns */
        {
            bm_get_pixel(img, j,i,&R, &G, &B);
            set_val(out,i,j,0,(float) R);
            set_val(out,i,j,1,(float) G);
            set_val(out,i,j,2,(float) B);
        }
    }

    return out;
}

Bitmap * ip_mat_to_bitmap(ip_mat * t){

    Bitmap *b = bm_create(t->w,t->h);

    unsigned int i, j;
    for (i = 0; i < t->h; i++)              /* rows */
    {
        for (j = 0; j < t->w; j++)          /* columns */
        {
            float c1 = clamp_f(get_val(t, i, j, 0), 0.0, 255.0);
            float c2 = clamp_f(get_val(t, i, j, 1), 0.0, 255.0);
            float c3 = clamp_f(get_val(t, i, j, 2), 0.0, 255.0);

            bm_set_pixel(b, j,i, (unsigned char)c1,
                    (unsigned char)c2,
                    (unsigned char)c3);
        }
    }
    return b;
}

float get_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k){
    if(i<a->h && j<a->w &&k<a->k){  /* j>=0 and k>=0 and i>=0 is non sense*/
        return a->data[k][i][j];
    }else{
        printf("Errore get_val!!!");
        exit(1);
    }
}

void set_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k, float v){
    if(i<a->h && j<a->w &&k<a->k){
        a->data[k][i][j]=v;
    }else{
        printf("Errore set_val!!!");
        exit(1);
    }
}

float get_normal_random(){
    float y1 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float y2 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    return cos(2*PI*y2)*sqrt(-2.*log(y1));

}

/*------------------------------------------prima parte------------------------------------------------------*/

/* Inizializza una ip_mat con dimensioni h w e k. Ogni elemento è inizializzato a v.
 * Inoltre crea un vettore di stats per contenere le statische sui singoli canali.
 * */
ip_mat * ip_mat_create(unsigned int h, unsigned int w,unsigned  int k, float v){ 
    int i, j, l;
    float ***p;
    ip_mat *matt;
    stats *statistiche;

    matt=(ip_mat *)malloc(sizeof(ip_mat));  /*inizio array 3d*/
    matt->h = h;
    matt->w = w;
    matt->k = k;

    p=(float***)malloc(k*sizeof(float**));

    for(i=0;i<k;i++){
        p[i]=(float**)malloc(h*sizeof(float*));
    }

    for (i = 0; i < k; i++){
        for (j = 0; j < h; j++){
            p[i][j] = (float*)malloc(w*sizeof(float));
        }
    }

    for (i = 0; i < k; i++){
       for (j = 0; j < h; j++){
            for (l = 0; l < w; l++){
                p[i][j][l] = v;
            }
        }
    }

    matt->data = p;  /*fine array 3d*/

    statistiche = (stats *)malloc(k*sizeof(stats)); /*inizio array stats*/

    for (i = 0; i < k; i++)
    {
        statistiche[i].max = v;
        statistiche[i].min = v;
        statistiche[i].mean = v;
    }

    matt->stat = statistiche; /*fine array stats*/
    
    return matt;
}


/* Libera la memoria (data, stat e la struttura) */
void ip_mat_free(ip_mat *a){
    int i,j;
    for (i = 0; i < a->k; i++){
       for (j = 0; j < a->h; j++){
            free(a->data[i][j]);
        }
    }
    for (i = 0; i < a->k; i++)  
    {
       free(a->data[i]);
    }
    free(a->data);
    free(a->stat);
}


/* Calcola il valore minimo, il massimo e la media per ogni canale
 * e li salva dentro la struttura ip_mat stats
 * */
void compute_stats(ip_mat * t)
{
	int i, j, l;
	int coeff=(t->w)*(t->h);
 
	for (i = 0; i < t->k; i++) 
    {
		float m;    
		float g;	
		float sum=0;
		m = t->data[i][0][0];
		g = t->data[i][0][0];
 
        for (j = 0; j < t->h; j++)
        {
            for (l = 0; l < t->w; l++)
            {
                int va;
                va = t->data[i][j][l];
                if (va > g)
					g = va;
				if (va < m)
					m = va;
				sum = sum + va;
            }
        }
		t->stat[i].min = m;
		t->stat[i].max = g;
		t->stat[i].mean = sum/coeff;
    }
}


/* Inizializza una ip_mat con dimensioni w h e k. w h k ci sono gia
 * Ogni elemento è generato da una gaussiana con media mean e varianza var */
void ip_mat_init_random(ip_mat * t, float mean, float var)
{
    int i, j, l;

    for (i = 0; i < t->k; i++){
        for (j = 0; j < t->h; j++){
            for (l = 0; l < t->w; l++){
                /*Applico formula della gaussiana con media mean e varianza var*/ 
                t->data[i][j][l] = 
                    (1.0 / (var * sqrt(2.0 * 3.1415))) 
                    * exp(-(pow((get_normal_random() - mean), 2) 
                        / (2.0 * pow(var, 2))));
            }
        }
    }
}

/* Crea una copia di una ip_mat e lo restituisce in output */
ip_mat * ip_mat_copy(ip_mat * in){
	int i, j, l;
    ip_mat *copy;
 
    copy = ip_mat_create(in->h, in->w, in->k, 0.0) ;/*inizio array 3d*/
 
    for (i = 0; i < in->k; i++) 
    {
       for (j = 0; j < in->h; j++) 
       {
            for (l = 0; l < in->w; l++) 
            {
               copy->data[i][j][l] = in->data[i][j][l];
            }
       }
    }/*fine array 3d*/
 
    for (i = 0; i < in->k; i++) 
    {   /*inizio array stats*/
       copy->stat[i].max = in->stat[i].max;
       copy->stat[i].min = in->stat[i].min;
       copy->stat[i].mean = in->stat[i].mean;
    }/*fine array stats*/     
  
    return copy;
}


/* Restituisce una sotto-matrice, ovvero la porzione individuata da:
 * t->data[row_start...row_end][col_start...col_end][0...k]
 * La terza dimensione la riportiamo per intero, stiamo in sostanza prendendo un sottoinsieme
 * delle righe e delle colonne.
 * */
ip_mat * ip_mat_subset(ip_mat * t, unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end){
    ip_mat *matt; 
    int i, j, l;

    matt=ip_mat_create(row_end - row_start, col_end - col_start, t->k, 0.0);

    for (i = 0; i < matt->k; i++){  /*copia la sotto matrice*/
        for (j = 0; j < matt->h; j++){
            for (l = 0; l < matt->w; l++){
                matt->data[i][j][l] =  t->data[i][row_start+j][col_start+l];
            }
        }
    }

    /*aggiorna stats in base alla nuova matrice*/
    compute_stats(matt); 
    
    return matt;
}

/* Concatena due ip_mat su una certa dimensione.
 * Ad esempio:
 * ip_mat_concat(ip_mat * a, ip_mat * b, 0);
 *      produrrà un nuovo ip_mat di dimensioni:
 *      out.h = a.h + b.h
 *      out.w = a.w = b.w
 *      out.k = a.k = b.k
 *
 * ip_mat_concat(ip_mat * a, ip_mat * b, 1);
 *      produrrà un nuovo ip_mat di dimensioni:
 *      out.h = a.h = b.h
 *      out.w = a.w + b.w
 *      out.k = a.k = b.k
 *
 * ip_mat_concat(ip_mat * a, ip_mat * b, 2);
 *      produrrà un nuovo ip_mat di dimensioni:
 *      out.h = a.h = b.h
 *      out.w = a.w = b.w
 *      out.k = a.k + b.k
 * */
ip_mat * ip_mat_concat(ip_mat * a, ip_mat * b, int dimensione){
    ip_mat *matt; 
    int i, j, l;

    if(dimensione == 0 && a->w == b->w && a->k == b->k){
        matt = ip_mat_create(a->h + b->h, a->w, a->k, 0.0);

        for (i = 0; i < matt->k; i++){
            for (j = 0; j < matt->h; j++){
                for (l = 0; l < matt->w; l++){
                    if(j < a->h)
                        matt->data[i][j][l] = a->data[i][j][l];
                    else
                        matt->data[i][j][l] = b->data[i][j][l];
                }
            }
        }
    
    }else if(dimensione == 1 && a->h == b->h && a->k == b->k){
        matt = ip_mat_create(a->h, a->w + b->w, a->k, 0.0);
    
        for (i = 0; i < matt->k; i++){
            for (j = 0; j < matt->h; j++){
                for (l = 0; l < matt->w; l++){
                    if(l < a->w)
                        matt->data[i][j][l] = a->data[i][j][l];
                    else
                        matt->data[i][j][l] = b->data[i][j][l];
                }
            }
        }
    
    }else if(dimensione == 2 && a->w == b->w && a->h == b->h){
        matt = ip_mat_create(a->h, a->w, a->k + b->k, 0.0);

        for (i = 0; i < matt->k; i++){
            for (j = 0; j < matt->h; j++){
                for (l = 0; l < matt->w; l++){
                    if(i < a->k)
                        matt->data[i][j][l] = a->data[i][j][l];
                    else
                        matt->data[i][j][l] = b->data[i][j][l];
                }
            }
        }
    }else{
        printf("Errore ip_mat_concat!!!");
        exit(1);
    }
    compute_stats(matt);
    return matt;
}

//I controlli per i valori dei pixel (clamp tra 0 e 255) vengono fatti su ip_mat_to_bitmap perciò non c'è bisogno di farli qui

/* Esegue la somma di due ip_mat (tutte le dimensioni devono essere identiche)
* e la restituisce in output. */
ip_mat * ip_mat_sum(ip_mat * a, ip_mat * b){
    if(a->h == b->h && a->w == b->w && a->k == b->k){
        int i, j, l;
        ip_mat *sum;
        sum = ip_mat_create(a->h, a->w, a->k, 0.0); 
    
        for (i = 0; i < a->k; i++) 
        {
            for (j = 0; j < a->h; j++) 
            {
                for (l = 0; l < a->w; l++) 
                {
                    sum->data[i][j][l] = a->data[i][j][l] + b->data[i][j][l];
                }
            }
        }/*fine array 3d*/
        compute_stats(sum);
        return sum;
    }else{
       printf("Errore ip_mat_sum!!!");
       exit(1);
   }
}



/* Esegue la sottrazione di due ip_mat (tutte le dimensioni devono essere identiche)
* e la restituisce in output.
* */
ip_mat * ip_mat_sub(ip_mat * a, ip_mat * b){
    if(a->h == b->h && a->w == b->w && a->k == b->k){
        int i, j, l;
        ip_mat *sub;    
       sub = ip_mat_create(a->h, a->w, a->k, 0.0); 
    
        for (i = 0; i < a->k; i++) {
            for (j = 0; j < a->h; j++) {
                for (l = 0; l < a->w; l++) {
                    sub->data[i][j][l] = a->data[i][j][l] + b->data[i][j][l];   
                }
            }
        }/*fine array 3d*/
        compute_stats(sub);
        return sub;
    }else{
       printf("Errore ip_mat_sub!!!");
       exit(1);
   }
}



/* Moltiplica un ip_mat per uno scalare c. Si moltiplica c per tutti gli elementi di "a"
 * e si salva il risultato in un nuovo tensore in output. */
ip_mat * ip_mat_mul_scalar(ip_mat *a, float c){
	int i, j, l;
	ip_mat *mulscalar;	
	mulscalar = ip_mat_create(a->h, a->w, a->k, 0.0);
 
	for (i = 0; i < a->k; i++) 
    {
		for (j = 0; j < a->h; j++) 
        {
			for (l = 0; l < a->w; l++) 
            {
                mulscalar->data[i][j][l] = a->data[i][j][l] * c;

			}
		}
	}
	compute_stats(mulscalar);
	return mulscalar;
}


/* Aggiunge ad un ip_mat uno scalare c e lo restituisce in un nuovo tensore in output. */
ip_mat *  ip_mat_add_scalar(ip_mat *a, float c){
   int i, j, l;
   ip_mat *addscalar; 
   addscalar = ip_mat_create(a->h, a->w, a->k, 0.0);
  
   for (i = 0; i < a->k; i++) 
   {
       for (j = 0; j < a->h; j++) 
       {
           for (l = 0; l < a->w; l++) 
           {
               addscalar->data[i][j][l] = a->data[i][j][l] + c;
           }
       }
   }
   compute_stats(addscalar);
   return addscalar;
}




/* Calcola la media di due ip_mat a e b e la restituisce in output.*/
ip_mat * ip_mat_mean(ip_mat * a, ip_mat * b){
    if(a->h == b->h && a->w == b->w && a->k == b->k){
        int i, j, l;
        ip_mat *mean;  
        mean= ip_mat_create(a->h, a->w, a->k, 0.0);
    
        for (i = 0; i < a->k; i++) 
        {
            for (j = 0; j < a->h; j++) 
            {
                for (l = 0; l < a->w; l++) 
                {
                    mean->data[i][j][l] = (a->data[i][j][l] + b->data[i][j][l])/2;
                }
            }
        }
        compute_stats(mean);
        return mean;
    }else{
       printf("Errore ip_mat_mean!!!");
       exit(1);
   }
}

/*------------------------------------------secondo parte------------------------------------------------------*/

/* Converte un'immagine RGB ad una immagine a scala di grigio.
 * Quest'operazione viene fatta calcolando la media per ogni pixel sui 3 canali
 * e creando una nuova immagine avente per valore di un pixel su ogni canale la media appena calcolata.
 * Avremo quindi che tutti i canali saranno uguali.
 * */
ip_mat * ip_mat_to_gray_scale(ip_mat * in){
    int i, j, l;
    ip_mat *grigio;  
    grigio = ip_mat_create(in->h, in->w, in->k, 0.0);

    for (l = 0; l < in->w; l++) {    /*calcola media di ogni pixel sui 3 canali e lo mette dentro ip_mat nuovo*/
        for (j = 0; j < in->h; j++) {
            float sum = 0;

            for (i = 0; i < in->k; i++) {
                sum += in->data[i][j][l]; /*sommo i valori del pixel sui 3 canali*/
            }

            for (i = 0; i < in->k; i++) {
                grigio->data[i][j][l] = sum / 3.0; /*metto sul pixel sui 3 canali il valore medio*/
            }
        }
    }
    compute_stats(grigio);
    return grigio;
}

/* Effettua la fusione (combinazione convessa) di due immagini */
ip_mat * ip_mat_blend(ip_mat * a, ip_mat * b, float alpha){
    if(a->h == b->h && a->w == b->w && a->k == b->k && alpha<=1.0 && alpha>=0.0)
    { 
        /*controllo dimensioni e alpha[0,1]*/
        ip_mat* alpha1 = ip_mat_mul_scalar(a, alpha);
        ip_mat* alpha2 = ip_mat_mul_scalar(b, 1 - alpha);
        ip_mat* blends = ip_mat_sum(alpha1, alpha2);

        compute_stats(blends);
        return blends;
    }
    else
    {
        printf("Errore ip_mat_blend!!!");
        exit(1);
    }
}

/* Operazione di brightening: aumenta la luminosità dell'immagine
 * aggiunge ad ogni pixel un certo valore*/
ip_mat * ip_mat_brighten(ip_mat * a, float bright){
	ip_mat *br;	
 
	br=ip_mat_copy(a);
	br=ip_mat_add_scalar(br,bright);
 
	compute_stats(br);
	return br;
}

/* Operazione di corruzione con rumore gaussiano:
 * Aggiunge del rumore gaussiano all'immagine, il rumore viene enfatizzato
 * per mezzo della variabile amount.
 * out = a + gauss_noise*amount
 * */
ip_mat * ip_mat_corrupt(ip_mat * a, float amount){

    if (amount >= 0.0 && amount <= 255){

        ip_mat* gaussNoise = ip_mat_create(a->h, a->w, a->k, 0);
        ip_mat_init_random(gaussNoise, get_normal_random(), get_normal_random());
        gaussNoise = ip_mat_mul_scalar(gaussNoise, amount);

        ip_mat* corr = ip_mat_sum(a, gaussNoise);

        compute_stats(corr);
        return corr;

    }else{
        printf("Errore ip_mat_corrupt!!!");
        exit(1);
    
    }

}

