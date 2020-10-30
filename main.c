#include <stdio.h>
#include <string.h>
#include <gmp.h>

mpz_t * doubled (mpz_t * p1, mpz_t * p2, mpz_t p){
    
    /*
    The "dbl" doubling formulas
     
    X3 = Y1*(Z1^3-X1^3)
    Y3 = X1*(Y1^3-Z1^3)
    Z3 = Z1*(X1^3-Y1^3)
     
    */
    
    mpz_t a, b, subzx, c, subyz, subxy;
    mpz_init(a);
    mpz_init(b);
    mpz_init(subzx);
    mpz_init(c);
    mpz_init(subyz);
    mpz_init(subxy);
    
    //X3
    
    mpz_mul(a, *(p1 + 2), *(p1 + 2));
    mpz_mul(a, a, *(p1 + 2));                   //Z1^3
    mpz_mul(b, *p1, *p1);
    mpz_mul(b, b, *p1);                         //X1^3
    mpz_sub(subzx, a, b);                       //Z1^3-X1^3
    mpz_mul(*p2, *(p1 + 1), subzx);             //X3
    
    //Y3
    
    mpz_mul(c, *(p1 + 1), *(p1 + 1));
    mpz_mul(c, c, *(p1 + 1));                   //Y1^3
    mpz_sub(subyz, c, a);                       //Y1^3-Z1^3
    mpz_mul(*(p2 + 1), *p1, subyz);             //Y3
    
    //Z3
    
    mpz_sub(subxy, b, c);                       //X1^3-Y1^3
    mpz_mul(*(p2 + 2), *(p1 + 2), subxy);       //Z3
    
    //mod p
    
    mpz_fdiv_r(*p2, *p2, p);
    mpz_fdiv_r(*(p2 + 1), *(p2 + 1), p);
    mpz_fdiv_r(*(p2 + 2), *(p2 + 2), p);
    
    return p2;
}

mpz_t * addition (mpz_t * p1, mpz_t * p2, mpz_t * p3, mpz_t p){
    
    mpz_t y1x2, y1y2, z1y2, z1z2, x1z2, x1x2, x31, x32, y31, y32, z31, z32;
    mpz_init(y1x2);
    mpz_init(y1y2);
    mpz_init(z1y2);
    mpz_init(z1z2);
    mpz_init(x1z2);
    mpz_init(x1x2);
    mpz_init(x31);
    mpz_init(x32);
    mpz_init(y31);
    mpz_init(y32);
    mpz_init(z31);
    mpz_init(z32);
        
    /*
    The "add-2009-bkl" addition formulas
     
    Y1X2 = Y1*X2
    Y1Y2 = Y1*Y2
    Z1Y2 = Z1*Y2
    Z1Z2 = Z1*Z2
    X1Z2 = X1*Z2
    X1X2 = X1*X2
    X3 = Z1Z2*Z1Y2-X1X2*Y1X2
    Y3 = Y1Y2*Y1X2-Z1Z2*X1Z2
    Z3 = X1X2*X1Z2-Y1Y2*Z1Y2
    */
    
    mpz_mul(y1x2, *(p1 + 1), *p2);          //Y1X2
    mpz_mul(y1y2, *(p1 + 1), *(p2 + 1));    //Y1Y2
    mpz_mul(z1y2, *(p1 + 2), *(p2 + 1));    //Z1Y2
    mpz_mul(z1z2, *(p1 + 2), *(p2 + 2));    //Z1Z2
    mpz_mul(x1z2, *p1, *(p2 + 2));          //X1Z2
    mpz_mul(x1x2, *p1, *p2);                //X1X2
    
    //X3
    
    mpz_mul(x31, z1z2, z1y2);
    mpz_mul(x32, x1x2, y1x2);
    mpz_sub(*p3, x31, x32);
    
    //Y3
    
    mpz_mul(y31, y1y2, y1x2);
    mpz_mul(y32, z1z2, x1z2);
    mpz_sub(*(p3 + 1), y31, y32);
    
    //Z3
    
    mpz_mul(z31, x1x2, x1z2);
    mpz_mul(z32, y1y2, z1y2);
    mpz_sub(*(p3 + 2), z31, z32);
    
    //mod p
    
    mpz_fdiv_r(*p3, *p3, p);
    mpz_fdiv_r(*(p3 + 1), *(p3 + 1), p);
    mpz_fdiv_r(*(p3 + 2), *(p3 + 2), p);
    
    return p3;
}

mpz_t * ladder_M (mpz_t * P, mpz_t * O, mpz_t * Q, mpz_t k, mpz_t p){
    
    char str[10000];
    char * string = str;
    
    string = mpz_get_str(string, 2, k);
    
    mpz_t r[3];
    mpz_init(r[0]);
    mpz_init(r[1]);
    mpz_init(r[2]);
    
    mpz_t * R = r;
    
    //Q = O
    
    mpz_set(*Q, *O);
    mpz_set(*(Q + 1), *(O + 1));
    mpz_set(*(Q + 2), *(O + 2));
    
    //R = P
    
    mpz_set(*R, *P);
    mpz_set(*(R + 1), *(P + 1));
    mpz_set(*(R + 2), *(P + 2));
    
    mpz_t r1[3];
    mpz_init(r1[0]);
    mpz_init(r1[1]);
    mpz_init(r1[2]);
    
    mpz_t * R1 = r1;
    
    mpz_t q1[3];
    mpz_init(q1[0]);
    mpz_init(q1[1]);
    mpz_init(q1[2]);
    
    mpz_t * Q1 = q1;
    
    long int i = 0;
    
    while (i < strlen(string)){
        
        if (*(string + i) == '0'){
            
            //R = R + Q и Q = [2]Q
            
            R1 = addition(R, Q, R1, p);
            
            mpz_set(*R, *R1);
            mpz_set(*(R + 1), *(R1 + 1));
            mpz_set(*(R + 2), *(R1 + 2));
            
            Q1 = doubled(Q, Q1, p);
            
            mpz_set(*Q, *Q1);
            mpz_set(*(Q + 1), *(Q1 + 1));
            mpz_set(*(Q + 2), *(Q1 + 2));
            
        }else{
            
            //Q = Q + R и R = [2]R
            
            Q1 = addition(Q, R, Q1, p);
            
            mpz_set(*Q, *Q1);
            mpz_set(*(Q + 1), *(Q1 + 1));
            mpz_set(*(Q + 2), *(Q1 + 2));
            
            R1 = doubled(R, R1, p);
            
            mpz_set(*R, *R1);
            mpz_set(*(R + 1), *(R1 + 1));
            mpz_set(*(R + 2), *(R1 + 2));
            
        }
        i = i + 1;
    }
    
    return Q;
    
}

//Test 1. Лежит ли точка на кривой.

void on_curve(mpz_t * Q, mpz_t p){
    
    mpz_t d, three;
    mpz_init_set_ui(d, 3);
    mpz_init_set_ui(three, 3);
    
    mpz_t x3, y3, z3, xy, full_left_side, dx, yz, dxyz, full_right_side;
    mpz_init(x3);
    mpz_init(y3);
    mpz_init(z3);
    mpz_init(xy);
    mpz_init(full_left_side);
    mpz_init(dx);
    mpz_init(yz);
    mpz_init(dxyz);
    mpz_init(full_right_side);
    
    mpz_mul(x3, *Q, *Q);
    mpz_mul(x3, x3, *Q);                      //x^3
    
    mpz_mul(y3, *(Q + 1), *(Q + 1));
    mpz_mul(y3, y3, *(Q + 1));                //y^3
    
    mpz_mul(z3, *(Q + 2), *(Q + 2));
    mpz_mul(z3, z3, *(Q + 2));                //z^3
    
    mpz_add(xy, x3, y3);
    mpz_add(full_left_side, xy, z3);          //x^3 + y^3 + z^3
    
    mpz_fdiv_r(full_left_side, full_left_side, p);
    
    mpz_mul(dx, d, *Q);
    mpz_mul(yz, *(Q + 1), *(Q + 2));
    mpz_mul(dxyz, dx, yz);
    mpz_mul(full_right_side, dxyz, three);    //3 * d * x * y * z
    
    mpz_fdiv_r(full_right_side, full_right_side, p);
    
    if (mpz_cmp(full_left_side, full_right_side) == 0) {
        
        printf("%s\n", "point is on the curve");
        
    }else{
        
        printf("%s\n", "point is not on the curve");
        
    }
    
}


int main(void) {
    
    mpz_t x;
    mpz_t y;
    mpz_t z;
    mpz_t p;
    mpz_t q;
    
    mpz_init_set_str(x, "93528601524582384978654134272795468222405403055717890280271688132874849008326", 10);
    mpz_init_set_str(y, "14443324612566128911211262381388707474030458136470034119105598903952521080679", 10);
    mpz_init_set_str(z, "1", 10);
    mpz_init_set_str(p, "115792089237316195423570985008687907853269984665640564039457584007913111864739", 10);
    mpz_init_set_str(q, "115792089237316195423570985008687907852907080286716537199505774922796921406320", 10);
    
    mpz_t P[3];
    
    mpz_init_set(P[0], x);
    mpz_init_set(P[1], y);
    mpz_init_set(P[2], z);
    
    mpz_t * pointerP = P;
    
    //Точка бесконечности
    
    mpz_t O[3];
    mpz_init_set_ui(O[0], 1);
    mpz_init_set_si(O[1], -1);
    mpz_init_set_ui(O[2], 0);
    
    mpz_t * pointerO = O;
    
    mpz_t P1[3];
    
    mpz_init(P1[0]);
    mpz_init(P1[1]);
    mpz_init(P1[2]);
    
    mpz_t * pointer1 = P1;
    
    mpz_t k;
    mpz_init_set_ui(k, 9);
    
    //Test 1. Лежит ли точка на кривой. Результат.
    
    printf("%s\n", "Test 1");
    on_curve(ladder_M(pointerP, pointerO, pointer1, k, p), p);
    
    //Test 2. Проверить, что [q]P = O, где q – порядок группы точек.
    
    printf("%s\n", "Test 2");
    printf("%s", "Если Z = 0, тогда это точка бесконечности, Z = ");
    gmp_printf("%Zd \n", *(ladder_M(pointerP, pointerO, pointer1, q, p) + 2));
    
    //Test 3. Проверить, что [q + 1]P = P и [q − 1]P = −P.
    
    mpz_set(k, q);
    mpz_add_ui(k, k, 1);
    
    printf("%s\n", "Проверим, что координаты пропорциональны [q + 1]P = P");
    
    mpz_mul(x, *P, *(ladder_M(pointerP, pointerO, pointer1, k, p) + 1));
    mpz_mul(y, *(P + 1), *(ladder_M(pointerP, pointerO, pointer1, k, p)));
    
    mpz_fdiv_r(x, x, p);
    mpz_fdiv_r(y, y, p);
    
    gmp_printf("%Zd \n""%Zd \n", x, y);
    
    mpz_sub_ui(k, q, 1);
    
    printf("%s\n", "Проверим, что координаты пропорциональны [q − 1]P = −P");
    
    mpz_mul(x, *P, *(ladder_M(pointerP, pointerO, pointer1, k, p)));
    mpz_mul(y, *(P + 1), *(ladder_M(pointerP, pointerO, pointer1, k, p) + 1));
    
    mpz_fdiv_r(x, x, p);
    mpz_fdiv_r(y, y, p);
    
    gmp_printf("%Zd \n""%Zd \n", x, y);
    
    //Test 4. Для двух случайных k1, k2 проверить, что [k1]P + [k2]P = [k1 + k2]P
    
    printf("%s\n", "Test 4");
    printf("%s \n", "Проверим, что координаты пропорциональны");
    
    mpz_t k1, k2;
    mpz_init_set_ui(k1, 5);
    mpz_init_set_ui(k2, 7);
    
    mpz_add(k, k1, k2);
    
    pointer1 = addition(ladder_M(pointerP, pointerO, pointer1, k1, p), ladder_M(pointerP, pointerO, pointer1, k2, p), pointer1, p);
    
    pointerP = ladder_M(pointerP, pointerO, pointer1, k, p);
    
    mpz_mul(x, *pointer1, *(pointerP + 1));
    mpz_mul(y, *(pointer1 + 1), *pointerP);
    
    mpz_fdiv_r(x, x, p);
    mpz_fdiv_r(y, y, p);
    
    printf("%s", "[k1]P + [k2]P = ");
    gmp_printf("%Zd \n", x);
    
    printf("%s", "[k1 + k2]P = ");
    gmp_printf("%Zd \n", y);
    
    return 0;
    
}
