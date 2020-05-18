/*
 Created by Sebastiano Vascon on 23/03/20.
*/

#include <stdio.h>
#include <math.h>
#include "ip_lib.h"
#include "bmp.h"


/* Comprime un valore all'interno del range [lo,hi]*/
float clamp_f(float f, float lo, float hi);

/*Ritorna una matrice ad una dimensione i cui valori sono la media delle k dimensioni della matrice originale*/
ip_mat* ip_mat_1D_average(ip_mat* t);


/*----------------------------------------parte base----------------------------------------*/

void ip_mat_show(ip_mat* t) {
    unsigned int i, l, j;
    printf("Matrix of size %d x %d x %d (hxwxk)\n", t->h, t->w, t->k);
    for (l = 0; l < t->k; l++) {
        printf("Slice %d\n", l);
        for (i = 0; i < t->h; i++) {
            for (j = 0; j < t->w; j++) {
                printf("%f ", get_val(t, i, j, l));
            }
            printf("\n");
        }
        printf("\n");
    }
}

void ip_mat_show_stats(ip_mat* t) {
    unsigned int k;

    compute_stats(t);

    for (k = 0; k < t->k; k++) {
        printf("Channel %d:\n", k);
        printf("\t Min: %f\n", t->stat[k].min);
        printf("\t Max: %f\n", t->stat[k].max);
        printf("\t Mean: %f\n", t->stat[k].mean);
    }
}

ip_mat* bitmap_to_ip_mat(Bitmap* img) {
    unsigned int i = 0, j = 0;

    unsigned char R, G, B;

    unsigned int h = img->h;
    unsigned int w = img->w;

    ip_mat* out = ip_mat_create(h, w, 3, 0);

    for (i = 0; i < h; i++)              /* rows */
    {
        for (j = 0; j < w; j++)          /* columns */
        {
            bm_get_pixel(img, j, i, &R, &G, &B);
            set_val(out, i, j, 0, (float)R);
            set_val(out, i, j, 1, (float)G);
            set_val(out, i, j, 2, (float)B);
        }
    }

    return out;
}

Bitmap* ip_mat_to_bitmap(ip_mat* t) {
    unsigned int i, j;

    Bitmap* b = bm_create(t->w, t->h);
    clamp(t, 0.0, 255.0);

    for (i = 0; i < t->h; i++)              /* rows */
    {
        for (j = 0; j < t->w; j++)          /* columns */
        {
            /*Non ho capito se questo controllo è al 100% cura dell'utente oppure va implementato
            qualora andassimo fuori range*/

            bm_set_pixel(b, j, i, (unsigned char)(get_val(t, i, j, 0)),
                (unsigned char)(get_val(t, i, j, 1)),
                (unsigned char)(get_val(t, i, j, 2)));
        }
    }
    return b;
}

float get_val(ip_mat* a, unsigned int i, unsigned int j, unsigned int k) {
    if (i < a->h && j < a->w && k < a->k) {  /* j>=0 and k>=0 and i>=0 is non sense*/
        return a->data[k][i][j];
    }
    else {
        printf("Errore get_val!!!\n");
        exit(1);
    }
}

void set_val(ip_mat* a, unsigned int i, unsigned int j, unsigned int k, float v) {
    if (i < a->h && j < a->w && k < a->k) {
        a->data[k][i][j] = v;
    }
    else {
        printf("Errore set_val!!!\n");
        exit(1);
    }
}

float get_normal_random(float media, float std) {

    float y1 = ((float)(rand()) + 1.) / ((float)(RAND_MAX)+1.);
    float y2 = ((float)(rand()) + 1.) / ((float)(RAND_MAX)+1.);
    float num = cos(2 * PI * y2) * sqrt(-2. * log(y1));

    return media + num * std;
}


/*----------------------------------------prima parte----------------------------------------*/

 /* Inizializza una ip_mat con dimensioni h w e k. Ogni elemento è inizializzato a v.
 * Inoltre crea un vettore di stats per contenere le statische sui singoli canali.
 * */
ip_mat* ip_mat_create(unsigned int h, unsigned int w, unsigned  int k, float v) {
    unsigned int i, j, l;
    float*** p;
    ip_mat* matt;
    stats* statistiche;

    /*Controllo tramite assert che le allocazioni siano andate a buon fine*/
    matt = (ip_mat*)malloc(sizeof(ip_mat));  /*inizio allocazione array 3d*/
    assert(matt);

    matt->h = h;
    matt->w = w;
    matt->k = k;

    p = (float***)malloc(k * sizeof(float**));
    assert(p);

    for (i = 0; i < k; i++) {       
        p[i] = (float**)malloc(h * sizeof(float*));
        assert(p[i]);
    }

    for (i = 0; i < k; i++) {
        for (j = 0; j < h; j++) {
            p[i][j] = (float*)malloc(w * sizeof(float));
            assert(p[i][j]);
        }
    }

    for (i = 0; i < k; i++) {
        for (j = 0; j < h; j++) {
            for (l = 0; l < w; l++) {
                p[i][j][l] = v;
            }
        }
    }

    matt->data = p;  /*fine allocazione array 3d*/

    statistiche = (stats*)malloc(k * sizeof(stats)); /*inizio allcoazione array stats*/
    assert(statistiche);

    for (i = 0; i < k; i++)
    {
        statistiche[i].max = v;
        statistiche[i].min = v;
        statistiche[i].mean = v;
    }

    matt->stat = statistiche; /*fine allocazione array stats*/

    return matt;
}

/* Libera la memoria (data, stat e la struttura) */
void ip_mat_free(ip_mat* a) {
    if(a!=NULL){
        unsigned int i, j;
        for (i = 0; i < a->k; i++) {
            for (j = 0; j < a->h; j++) {
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
}

/* Calcola il valore minimo, il massimo e la media per ogni canale
 * e li salva dentro la struttura ip_mat stats
 * */
void compute_stats(ip_mat* t)
{
    unsigned int i, j, l;
    int coeff = (t->w) * (t->h);

    for (i = 0; i < t->k; i++)
    {
        float m;    /*minimo*/
        float g;    /*massimo*/
        float sum = 0;
        m = t->data[i][0][0];
        g = t->data[i][0][0];

        for (j = 0; j < t->h; j++)
        {
            for (l = 0; l < t->w; l++)
            {
                int va;
                va = t->data[i][j][l]; /*confronto dei valori*/
                if (va > g)
                    g = va;
                if (va < m)
                    m = va;
                sum = sum + va;
            }
        }
        t->stat[i].min = m;
        t->stat[i].max = g;
        t->stat[i].mean = sum / coeff;
    }
}

/* Inizializza una ip_mat con dimensioni w h e k. w h k ci sono gia
 * Ogni elemento è generato da una gaussiana con media mean e varianza var */
void ip_mat_init_random(ip_mat* t, float mean, float std)
{
    unsigned int i, j, l;

    for (i = 0; i < t->k; i++) {
        for (j = 0; j < t->h; j++) {
            for (l = 0; l < t->w; l++) {
                t->data[i][j][l] = get_normal_random_param(mean, var);
            }
        }
    }
}

/* Crea una copia di una ip_mat e lo restituisce in output */
ip_mat* ip_mat_copy(ip_mat* in) {
    unsigned int i, j, l;
    ip_mat* copy;

    copy = ip_mat_create(in->h, in->w, in->k, 0.0);/*inizio array 3d*/

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
ip_mat* ip_mat_subset(ip_mat* t, unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end) {
    ip_mat* matt;
    unsigned int i, j, l;

    matt = ip_mat_create(row_end - row_start, col_end - col_start, t->k, 0.0);

    for (i = 0; i < matt->k; i++) {  /*copia la sotto matrice*/
        for (j = 0; j < matt->h; j++) {
            for (l = 0; l < matt->w; l++) {
                matt->data[i][j][l] = t->data[i][row_start + j][col_start + l];
            }
        }
    }

    /*aggiorna stats in base alla nuova matrice*/
    compute_stats(matt);

    return matt;
}

/* Concatena due ip_mat su una certa dimensione.*/
ip_mat* ip_mat_concat(ip_mat* a, ip_mat* b, int dimensione) {
    ip_mat* matt;
    unsigned int i, j, l;

    if (dimensione == 0 && a->w == b->w && a->k == b->k) { /*concatena su h*/
        matt = ip_mat_create(a->h + b->h, a->w, a->k, 0.0);

        for (i = 0; i < matt->k; i++) {
            for (j = 0; j < matt->h; j++) {
                for (l = 0; l < matt->w; l++) {
                    if (j < a->h)
                        matt->data[i][j][l] = a->data[i][j][l];
                    else
                        matt->data[i][j][l] = b->data[i][j-a->h][l];
                }
            }
        }
    }
    else if (dimensione == 1 && a->h == b->h && a->k == b->k) {/*concatena su w*/
        matt = ip_mat_create(a->h, a->w + b->w, a->k, 0.0);

        for (i = 0; i < matt->k; i++) {
            for (j = 0; j < matt->h; j++) {
                for (l = 0; l < matt->w; l++) {
                    if (l < a->w)
                        matt->data[i][j][l] = a->data[i][j][l];
                    else
                        matt->data[i][j][l] = b->data[i][j][l-a->w];
                }
            }
        }
    }
    else if (dimensione == 2 && a->w == b->w && a->h == b->h) {/*concatena su k*/
        matt = ip_mat_create(a->h, a->w, a->k + b->k, 0.0);

        for (i = 0; i < matt->k; i++) {
            for (j = 0; j < matt->h; j++) {
                for (l = 0; l < matt->w; l++) {
                    if (i < a->k)
                        matt->data[i][j][l] = a->data[i][j][l];
                    else
                        matt->data[i][j][l] = b->data[i-a->k][j][l];
                }
            }
        }
    }
    else {
        printf("Errore ip_mat_concat!!!\n");
        exit(1);
    }

    compute_stats(matt);
    return matt;
}

/* Esegue la somma di due ip_mat (tutte le dimensioni devono essere identiche)
* e la restituisce in output. */
ip_mat* ip_mat_sum(ip_mat* a, ip_mat* b) {
    if (a->h == b->h && a->w == b->w && a->k == b->k) {
        unsigned int i, j, l;
        ip_mat* sum;
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
        }
        compute_stats(sum);
        return sum;
    }
    else {
        printf("Errore ip_mat_sum!!!\n");
        exit(1);
    }
}

/* Esegue la sottrazione di due ip_mat (tutte le dimensioni devono essere identiche)
* e la restituisce in output.
* */
ip_mat* ip_mat_sub(ip_mat* a, ip_mat* b) {
    if (a->h == b->h && a->w == b->w && a->k == b->k) {
        unsigned int i, j, l;
        ip_mat* sub;
        sub = ip_mat_create(a->h, a->w, a->k, 0.0);

        for (i = 0; i < a->k; i++)
        {
            for (j = 0; j < a->h; j++)  
            {
                for (l = 0; l < a->w; l++) 
                {
                    sub->data[i][j][l] = a->data[i][j][l] - b->data[i][j][l];
                }
            }
        }
        compute_stats(sub);
        return sub;
    }
    else {
        printf("Errore ip_mat_sub!!!\n");
        exit(1);
    }
}

/* Moltiplica un ip_mat per uno scalare c. Si moltiplica c per tutti gli elementi di "a"
 * e si salva il risultato in un nuovo tensore in output. */
ip_mat* ip_mat_mul_scalar(ip_mat* a, float c) {
    unsigned int i, j, l;
    ip_mat* mulscalar;
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
ip_mat* ip_mat_add_scalar(ip_mat* a, float c) {
    unsigned int i, j, l;
    ip_mat* addscalar;
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
ip_mat* ip_mat_mean(ip_mat* a, ip_mat* b) {
    if (a->h == b->h && a->w == b->w && a->k == b->k) {
        unsigned int i, j, l;
        ip_mat* mean;
        mean = ip_mat_create(a->h, a->w, a->k, 0.0);

        for (i = 0; i < a->k; i++)
        {
            for (j = 0; j < a->h; j++)
            {
                for (l = 0; l < a->w; l++)
                {
                    mean->data[i][j][l] = (a->data[i][j][l] + b->data[i][j][l]) / 2.0;
                }
            }
        }
        compute_stats(mean);
        return mean;
    }
    else {
        printf("Errore ip_mat_mean!!!\n");
        exit(1);
    }
}



/*----------------------------------------seconda parte----------------------------------------*/

/* Converte un'immagine RGB ad una immagine a scala di grigio.
 * Quest'operazione viene fatta calcolando la media per ogni pixel sui 3 canali
 * e creando una nuova immagine avente per valore di un pixel su ogni canale la media appena calcolata.
 * Avremo quindi che tutti i canali saranno uguali.
 * */
ip_mat* ip_mat_to_gray_scale(ip_mat* in) {
    unsigned int i, j, l;
    ip_mat* grigio;
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
ip_mat* ip_mat_blend(ip_mat* a, ip_mat* b, float alpha) {
    if (a->h == b->h && a->w == b->w && a->k == b->k && alpha <= 1.0 && alpha >= 0.0)
    {
        ip_mat* alpha1;
        ip_mat* alpha2;
        ip_mat* blends;

        /*controllo dimensioni e alpha[0,1]*/
        alpha1 = ip_mat_mul_scalar(a, alpha);
        alpha2 = ip_mat_mul_scalar(b, 1 - alpha);
        blends = ip_mat_sum(alpha1, alpha2);

        ip_mat_free(alpha1);
        ip_mat_free(alpha2);

        compute_stats(blends);
        return blends;
    }
    else
    {
        printf("Errore ip_mat_blend!!!\n");
        exit(1);
    }
}

/* Operazione di brightening: aumenta la luminosità dell'immagine
 * aggiunge ad ogni pixel un certo valore*/
ip_mat* ip_mat_brighten(ip_mat* a, float bright)
{
    ip_mat* tmp;
    ip_mat* br;

    tmp = ip_mat_copy(a);
    br = ip_mat_add_scalar(tmp, bright);

    compute_stats(br);
    return br;
}

/* Operazione di corruzione con rumore gaussiano:
 * Aggiunge del rumore gaussiano all'immagine, il rumore viene enfatizzato
 * per mezzo della variabile amount.
 * out = a + gauss_noise*amount
 * */
ip_mat* ip_mat_corrupt(ip_mat* a, float amount) {
    if (amount >= 0.0 && amount <= 255) {
        ip_mat* corrupted;
        corrupted = ip_mat_create(a->h, a->w, a->k, 0);  

        /*Con std = 1/3 amount, abbiamo il 99.73% di generare numeri compresi tra [-amount, +amount],
            in quanto amount = 3 std*/
        ip_mat_init_random(corrupted, 0, 0.33 * amount);

        compute_stats(corrupted);
        return corrupted;
    }
    else {
        printf("Errore ip_mat_corrupt!!!\n");
        exit(1);

    }

}


/*----------------------------------------terza parte----------------------------------------*/

/*Ritorna una matrice ad una dimensione i cui valori sono la media delle k dimensioni della matrice originale*/
/*Utilizzato nella prima versione del progetto per risolvere il problema di avere filtro ed immagine con un numero
    di canali diversi all'interno di ip_mat_convolve*/
ip_mat* ip_mat_1D_average(ip_mat* t)
{
    unsigned int i, j, k;
    ip_mat* avg;

    avg = ip_mat_create(t->h, t->w, 1, 0);

    for (i = 0; i < t->h; i++)
    {
        for (j = 0; j < t->w; j++)
        {
            float valSum = 0;

            for (k = 0; k < t->k; k++)
            {
                valSum += t->data[k][i][j];
            }

            /*Media delle k dimensioni*/
            avg->data[0][i][j] = valSum / t->k;
        }
    }

    compute_stats(avg);
    return avg;
}

/* Effettua la convoluzione di un ip_mat "a" con un ip_mat "f".
 * La funzione restituisce un ip_mat delle stesse dimensioni di "a".
 * */
ip_mat* ip_mat_convolve(ip_mat* a, ip_mat* f)
{
    /*Filtro ed immagine devono avere lo stesso numero di canali per poter correttamente applicare la convoluzione*/
    assert(a->k == f->k);

    unsigned int i, j, k, ii, jj;
    float pad_h, pad_w;
    ip_mat* padded;
    ip_mat* res;

    pad_h = (f->h - 1) / 2;
    pad_w = (f->w - 1) / 2;

    padded = ip_mat_padding(a, pad_h, pad_w);
    res = ip_mat_create(a->h, a->w, a->k, 0);

    for (k = 0; k < a->k; k++)
    {
        /*Bisogna far "scorrere" il filtro sulla matrice*/
        for (i = 0; i < padded->h - (pad_h * 2); i++)
        {
            for (j = 0; j < padded->w - (pad_w * 2); j++)
            {
                float result = 0.0;
                ip_mat* tmp;
                tmp = ip_mat_subset(padded, i, i + avgFilter->h, j, j + avgFilter->w);

                for (ii = 0; ii < avgFilter->h; ii++)
                {
                    for (jj = 0; jj < avgFilter->w; jj++)
                    {
                        /*Di default utilizzo sempre il primo canale del filtro fin quando non viene implementata
                        una soluzione per gestire canali differenti*/
                        result += (tmp->data[k][ii][jj] * avgFilter->data[0][ii][jj]);
                    }
                }

                res->data[k][i][j] = result;
                ip_mat_free(tmp);
            }
        }
    }

    ip_mat_free(padded);

    compute_stats(res);
    return res;
}

/* Aggiunge un padding all'immagine. Il padding verticale è pad_h mentre quello
 * orizzontale è pad_w.
 * L'output sarà un'immagine di dimensioni:
 *      out.h = a.h + 2*pad_h;
 *      out.w = a.w + 2*pad_w;
 *      out.k = a.k
 * con valori nulli sui bordi corrispondenti al padding e l'immagine "a" riportata
 * nel centro
 * */
ip_mat* ip_mat_padding(ip_mat* a, int pad_h, int pad_w)
{
    unsigned int i, j, k, l;
    ip_mat* padded;
    padded = ip_mat_create(a->h + pad_h * 2, a->w + pad_w * 2, a->k, 0);

    /*Creo prima una matrice di 0 delle dimensioni paddate (a.h + pad_h, a.w + pad_w),
    poi la "riempio" inserendo la matrice originale al centro
    padded = ip_mat_create(a->h + pad_h * 2, a->w + pad_w * 2, a->k, 0);*/

    for (k = 0; k < a->k; k++)
    {
        /*Bisogna far "scorrere" il filtro sulla matrice*/
        for (i = pad_h; i < (padded->h - pad_h); i++)
        {
            for (j = pad_w; j < (padded->w - pad_w); j++)
            {
                padded->data[k][i][j] = a->data[k][i - pad_h][j - pad_w];
            }
        }
    }

    compute_stats(padded);
    return padded;
}

/* Crea un filtro di sharpening */
ip_mat* create_sharpen_filter()
{
    unsigned int i;
    ip_mat* f;

    /*Filtro sharpen:
        0  -1  0
       -1   5 -1
        0  -1  0 */
    f = ip_mat_create(3, 3, 1, 0);

    for (i = 0; i < f->k; i++)
    {
        /*Prima e ultima riga basta aggiungere il -1 al centro perché è già piena di 0*/
        f->data[i][0][1] = -1;
        f->data[i][2][1] = -1;

        /*Riga centrale*/
        f->data[i][1][0] = -1;
        f->data[i][1][1] = 5;
        f->data[i][1][2] = -1;
    }

    return f;
}

/* Crea un filtro per rilevare i bordi */
ip_mat* create_edge_filter()
{
    unsigned int i;
    ip_mat* f;

    /*Filtro edge:
    -1  -1  -1
    -1   8  -1
    -1  -1  -1 */
    f = ip_mat_create(3, 3, 1, -1);

    for (i = 0; i < f->k; i++)
    {
        /*Matrice già riempita di -1, basta inserire 8 in centro*/
        f->data[i][1][1] = 8;
    }

    return f;
}

/* Crea un filtro per aggiungere profondità */
ip_mat* create_emboss_filter()
{
    unsigned int i;
    ip_mat* f;

    /*Filtro emboss:
    -2  -1   0
    -1   1   1
     0   1   2 */

     /*Matrice riempita di 1*/
    f = ip_mat_create(3, 3, 1, 1);

    for (i = 0; i < f->k; i++)
    {
        /*Prima riga*/
        f->data[i][0][0] = -2;
        f->data[i][0][1] = -1;
        f->data[i][0][2] = 0;

        /*Seconda riga*/
        f->data[i][1][0] = -1;

        /*Terza riga*/
        f->data[i][2][0] = 0;
        f->data[i][2][1] = 1;
        f->data[i][2][2] = 2;
    }

    return f;
}

/* Crea un filtro medio per la rimozione del rumore */
ip_mat* create_average_filter(unsigned int h, unsigned int w, unsigned int k)
{
    float c = 1.0f / (w * h);

    /*Ritorna una matrice "riempita" di c*/
    return ip_mat_create(h, w, k, c);
}

/* Crea un filtro gaussiano per la rimozione del rumore */
ip_mat* create_gaussian_filter(unsigned int h, unsigned int w, unsigned int k, float sigma)
{
    if(sigma>0 && (w%2==1) && (h%2==1)){
        unsigned int i, j, kk;
        ip_mat* f;
        ip_mat* normalized;

        /*Controllo che dimensione kernel sia dispari altrimenti non c'è il centro*/

        /*Coordinate centro*/
        unsigned int cx = (w - 1) / 2;
        unsigned int cy = (h - 1) / 2;
        float sum = 0.0;

        f = ip_mat_create(w, h, k, 0);

        for (i = 0; i < h; i++)
        {
            for (j = 0; j < w; j++)
            {
                /*Calcolo il valore gaussiano e aggiorno la somma,
                    che utilizzerò dopo per normalizzare i valori*/
                int xDist = i - cx;
                int yDist = j - cy;
                float gaussVal;
                gaussVal = (1 / (2 * PI * pow(sigma, 2))) * exp(-((pow(xDist, 2) + pow(yDist, 2)) / (2 * pow(sigma, 2))));

                sum += gaussVal;

                for (kk = 0; kk < k; kk++)
                {
                    f->data[kk][i][j] = gaussVal;
                }
            }
        }

        /*Normalizzo dividendo per la somma, in modo che, successivamente,
            la somma di tutti i valori del filtro sia unitaria*/
        normalized = ip_mat_mul_scalar(f, 1.0f / sum);

        ip_mat_free(f);

        compute_stats(normalized);
        return normalized;
    }
    else {
        printf("Errore parametri Gauss!!!\n");
        exit(1);
    }
}


/*-----------------------parte extra-----------------------*/

/* Comprime un valore all'interno del range [lo,hi]*/
float clamp_f(float f, float lo, float hi)
{
    if (f > hi)
    {
        return hi;
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

/* Nell'operazione di clamping i valori <low si convertono in low e i valori >high in high.*/
void clamp(ip_mat* t, float low, float high)
{
    unsigned int i, j, k;

    for (k = 0; k < t->k; k++)
    {
        for (i = 0; i < t->h; i++)
        {
            for (j = 0; j < t->w; j++)
            {
                t->data[k][i][j] = clamp_f(t->data[k][i][j], low, high);
            }
        }
    }
}

/* Effettua una riscalatura dei dati tale che i valori siano in [0,new_max].
 * Utilizzate il metodo compute_stat per ricavarvi il min, max per ogni canale.
 *
 * I valori sono scalati tramite la formula (valore - min) / (max - min)
 *
 * Si considera ogni indice della terza dimensione indipendente, quindi l'operazione
 * di scalatura va ripetuta per ogni "fetta" della matrice 3D.
 * Successivamente moltiplichiamo per new_max gli elementi della matrice in modo da ottenere un range
 * di valori in [0,new_max].
 * */
void rescale(ip_mat* t, float new_max)
{
    unsigned int i, j, k;

    /*Necessario per avere min e max*/
    compute_stats(t);

    for (k = 0; k < t->k; k++)
    {
        float min = t->stat[k].min;
        float max = t->stat[k].max;

        for (i = 0; i < t->h; i++)
        {
            for (j = 0; j < t->w; j++)
            {
                float val = t->data[k][i][j];
                t->data[k][i][j] = (val - min) / (max - min) * new_max;
            }
        }
    }
    compute_stats(t);
}

